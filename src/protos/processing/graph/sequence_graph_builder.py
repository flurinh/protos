"""
SequenceGraphBuilder: Build graphs from GRN tables and embeddings.

This module provides structure-free graph construction using:
- GRN assignments from sequence alignment
- Per-residue embeddings (ankh_large)
- Binding domain configuration for edges
- Positional encoding mapping from unified_grn_mapping.csv
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from protos.processing.grn.binding_domain import BindingDomainConfig
from protos.processing.grn.grn_utils import normalize_grn_format, sort_grns_str

logger = logging.getLogger(__name__)

# Default paths relative to project root
DEFAULT_BINDING_CONFIG = "data/grn/configs/binding_domain.json"
DEFAULT_POS_ENCODING = "data/grn/configs/unified_grn_mapping.csv"
DEFAULT_GRN_TABLES = {
    "vpod1_3": "data/grn/tables/vpod1_3_generated.csv",
    "bistable": "data/grn/tables/bistable_generated.csv",
    "inoue": "data/grn/tables/inoue_generated.csv",
    "type_ii_opsin_atlas": "data/grn/tables/type_ii_opsin_atlas_generated.csv",
    "type_i_opsin_atlas": "data/grn/tables/type_i_opsin_atlas_generated.csv",
}


class SequenceGraphBuilder:
    """
    Build PyG-compatible graphs from GRN assignments and embeddings.

    This class enables structure-free graph construction using:
    - GRN table: protein → GRN position assignments
    - Embeddings: per-residue transformer embeddings (ankh_large)
    - Binding domain config: defines which GRN pairs form edges
    - Positional encoding mapping: GRN → position index (0-based)

    Output graphs have:
    - x: Node features (embeddings at GRN positions)
    - edge_index: Edges from binding domain configuration
    - pos: Positional encodings (0-based indices, -1 for unknown)
    - grn_str: Comma-separated GRN string
    - protein_id, protein_family: Metadata

    Example usage:
        builder = SequenceGraphBuilder(
            binding_config_path="data/grn/configs/binding_domain.json",
            positional_encoding_path="data/grn/configs/unified_grn_mapping.csv",
        )

        grn_table = pd.read_csv("data/grn/tables/vpod1_3_generated.csv", index_col=0)
        embedding = load_embedding(entity_id)  # shape: (seq_len, 1536)

        graph = builder.build_graph(
            protein_id=entity_id,
            grn_table=grn_table,
            embedding=embedding,
            protein_family="type_ii",
        )
    """

    def __init__(
        self,
        binding_config_path: Optional[str] = None,
        positional_encoding_path: Optional[str] = None,
        protein_family: str = "type_ii",
        binding_config: Optional[BindingDomainConfig] = None,
    ) -> None:
        """
        Initialize SequenceGraphBuilder.

        Args:
            binding_config_path: Path to binding domain JSON config
            positional_encoding_path: Path to unified_grn_mapping.csv
            protein_family: Default protein family (type_i, type_ii, hcrbpii)
            binding_config: Pre-loaded BindingDomainConfig (alternative to path)
        """
        self.protein_family = protein_family

        # Load binding domain configuration
        if binding_config is not None:
            self.binding_config = binding_config
        elif binding_config_path:
            self.binding_config = BindingDomainConfig(config_path=binding_config_path)
        else:
            self.binding_config = None
            logger.warning("No binding config provided - edges will not be built")

        # Load positional encoding mapping
        self._grn_to_pos: Dict[str, Dict[str, int]] = {}  # family -> grn -> pos

        if positional_encoding_path:
            self.load_positional_encoding(positional_encoding_path)

    def load_positional_encoding(self, path: str) -> None:
        """
        Load positional encoding mapping from unified_grn_mapping.csv.

        The CSV format has columns: type_i, type_ii, hcrbpii, cluster_size, paired
        For each GRN in a family column, the positional encoding is the row index (0-based).

        Args:
            path: Path to unified_grn_mapping.csv
        """
        csv_path = Path(path)
        if not csv_path.exists():
            raise FileNotFoundError(f"Positional encoding file not found: {path}")

        # Read CSV - GRN values are strings
        pos_df = pd.read_csv(csv_path)

        # Build GRN -> position mapping for each protein family
        # Position is the 0-based row index
        self._grn_to_pos = {}
        for family in ["type_i", "type_ii", "hcrbpii"]:
            if family not in pos_df.columns:
                continue
            self._grn_to_pos[family] = {}
            for idx, grn in enumerate(pos_df[family]):
                if pd.notna(grn) and str(grn).strip() != "":
                    # Store as string key, 0-based position
                    self._grn_to_pos[family][str(grn)] = idx

        logger.info(
            f"Loaded positional encoding from {path}: "
            + ", ".join(f"{k}: {len(v)} GRNs" for k, v in self._grn_to_pos.items())
        )

    def get_positional_encoding(
        self,
        grns: List[str],
        protein_family: Optional[str] = None,
    ) -> np.ndarray:
        """
        Get positional encoding indices for a list of GRNs.

        Args:
            grns: List of GRN strings
            protein_family: Protein family to use for mapping

        Returns:
            Array of 0-based position indices, -1 for unknown GRNs
        """
        family = protein_family or self.protein_family

        if family not in self._grn_to_pos:
            logger.warning(f"Unknown protein family: {family}, using -1 for all positions")
            return np.full(len(grns), -1, dtype=np.int64)

        grn_map = self._grn_to_pos[family]
        positions = []

        for grn in grns:
            # Try exact match first
            if grn in grn_map:
                positions.append(grn_map[grn])
            else:
                # GRN not in mapping - use -1 (model handles this as padding)
                positions.append(-1)

        return np.array(positions, dtype=np.int64)

    def get_assigned_grns(
        self,
        protein_id: str,
        grn_table: pd.DataFrame,
    ) -> List[str]:
        """
        Extract assigned GRN positions for a protein from the GRN table.

        GRN table format: rows are proteins, columns are GRN positions.
        Values are "AA#" where AA is amino acid and # is sequence position,
        or "-" for gaps.

        Args:
            protein_id: Protein identifier (row index)
            grn_table: DataFrame with GRN assignments

        Returns:
            List of assigned (non-gap) GRN column names
        """
        if protein_id not in grn_table.index:
            raise ValueError(f"Protein {protein_id} not found in GRN table")

        row = grn_table.loc[protein_id]
        assigned = []

        for col in grn_table.columns:
            val = row[col]
            if pd.notna(val) and val != "-" and val != "X":
                assigned.append(col)

        return assigned

    def get_sequence_positions(
        self,
        protein_id: str,
        grns: List[str],
        grn_table: pd.DataFrame,
    ) -> List[int]:
        """
        Get 1-based sequence positions for GRNs.

        Extracts the sequence position number from GRN table values
        (format: "AA#" where # is the sequence position).

        Args:
            protein_id: Protein identifier
            grns: List of GRN strings to get positions for
            grn_table: DataFrame with GRN assignments

        Returns:
            List of 1-based sequence positions
        """
        if protein_id not in grn_table.index:
            raise ValueError(f"Protein {protein_id} not found in GRN table")

        row = grn_table.loc[protein_id]
        positions = []

        for grn in grns:
            if grn in row.index:
                val = row[grn]
                if pd.notna(val) and val != "-" and val != "X":
                    # Extract position number (e.g., "A123" -> 123)
                    pos_str = "".join(c for c in str(val) if c.isdigit())
                    if pos_str:
                        positions.append(int(pos_str))
                    else:
                        positions.append(-1)  # Invalid format
                else:
                    positions.append(-1)  # Gap
            else:
                positions.append(-1)  # GRN not in table

        return positions

    def build_graph(
        self,
        protein_id: str,
        grn_table: pd.DataFrame,
        embedding: np.ndarray,
        protein_family: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Build a graph for a single protein.

        Args:
            protein_id: Protein identifier
            grn_table: DataFrame with GRN assignments (from *_generated.csv)
            embedding: Per-residue embedding array (seq_len, embed_dim)
            protein_family: Protein family (type_i, type_ii, hcrbpii)

        Returns:
            Dictionary with graph components:
            - x: Node features (n_nodes, embed_dim)
            - edge_index: Edge connectivity (2, n_edges)
            - pos: Positional encodings (n_nodes,) - 0-based, -1 for unknown
            - grn_str: Comma-separated GRN string
            - grn_positions: List of GRN strings (same as grn_str split)
            - protein_id: Protein identifier
            - protein_family: Protein family name
            - n_nodes: Number of nodes
            - n_edges: Number of edges
        """
        family = protein_family or self.protein_family

        # 1. Get assigned GRNs from table (non-gap positions)
        assigned_grns = self.get_assigned_grns(protein_id, grn_table)

        # 2. Filter to required GRNs (those in binding domain config)
        if self.binding_config is not None:
            required_grns = self.binding_config.get_required_grns(family)
            # Filter to GRNs that are both assigned and required
            # Use normalized comparison for matching
            required_set = set(required_grns)
            filtered_grns = [g for g in assigned_grns if normalize_grn_format(g) in required_set]
        else:
            filtered_grns = assigned_grns

        # 3. Sort GRNs in canonical order
        sorted_grns = sort_grns_str([normalize_grn_format(g) for g in filtered_grns])

        if len(sorted_grns) == 0:
            logger.warning(f"No GRNs for protein {protein_id} after filtering")
            return self._empty_graph(protein_id, family)

        # 4. Get sequence positions for embedding extraction
        # Map normalized GRNs back to original column names for table lookup
        norm_to_orig = {normalize_grn_format(g): g for g in filtered_grns}
        orig_grns = [norm_to_orig.get(g, g) for g in sorted_grns]
        seq_positions = self.get_sequence_positions(protein_id, orig_grns, grn_table)

        # 5. Extract node features from embedding
        # Embedding shape: (seq_len, embed_dim)
        # Sequence positions are 1-indexed, convert to 0-based for array indexing
        embed_dim = embedding.shape[-1]
        x = np.zeros((len(sorted_grns), embed_dim), dtype=np.float32)

        for i, pos in enumerate(seq_positions):
            emb_idx = pos - 1  # Convert 1-indexed to 0-indexed
            if 0 <= emb_idx < embedding.shape[0]:
                x[i] = embedding[emb_idx]
            else:
                logger.debug(f"Invalid position {pos} for GRN {sorted_grns[i]} in {protein_id}")

        # 6. Get positional encodings (0-based)
        if self._grn_to_pos:
            pos_enc = self.get_positional_encoding(sorted_grns, family)
        else:
            # Fallback: sequential 0-based indices
            pos_enc = np.arange(len(sorted_grns), dtype=np.int64)

        # 7. Build edges from binding domain config
        if self.binding_config is not None:
            src, dst = self.binding_config.build_edge_index(sorted_grns, family)
            edge_index = np.array([src, dst], dtype=np.int64)
        else:
            edge_index = np.zeros((2, 0), dtype=np.int64)

        return {
            "x": x,
            "edge_index": edge_index,
            "pos": pos_enc,
            "grn_str": ",".join(sorted_grns),
            "grn_positions": sorted_grns,
            "protein_id": protein_id,
            "protein_family": family,
            "n_nodes": len(sorted_grns),
            "n_edges": edge_index.shape[1],
        }

    def _empty_graph(self, protein_id: str, family: str) -> Dict[str, Any]:
        """Create an empty graph placeholder."""
        return {
            "x": np.zeros((0, 1), dtype=np.float32),
            "edge_index": np.zeros((2, 0), dtype=np.int64),
            "pos": np.zeros((0,), dtype=np.int64),
            "grn_str": "",
            "grn_positions": [],
            "protein_id": protein_id,
            "protein_family": family,
            "n_nodes": 0,
            "n_edges": 0,
        }

    def build_graphs(
        self,
        grn_table: pd.DataFrame,
        embeddings: Dict[str, np.ndarray],
        protein_ids: Optional[List[str]] = None,
        protein_family: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Build graphs for multiple proteins.

        Args:
            grn_table: DataFrame with GRN assignments
            embeddings: Dict mapping protein_id -> per-residue embedding
            protein_ids: List of proteins to process (defaults to all in table)
            protein_family: Protein family for all graphs

        Returns:
            List of graph dictionaries
        """
        if protein_ids is None:
            protein_ids = list(grn_table.index)

        graphs = []
        for pid in protein_ids:
            if pid not in embeddings:
                logger.warning(f"No embedding for {pid}, skipping")
                continue

            try:
                graph = self.build_graph(
                    protein_id=pid,
                    grn_table=grn_table,
                    embedding=embeddings[pid],
                    protein_family=protein_family,
                )
                graphs.append(graph)
            except Exception as e:
                logger.error(f"Error building graph for {pid}: {e}")
                continue

        return graphs

    def to_pyg_data(self, graph: Dict[str, Any]) -> Any:
        """
        Convert graph dictionary to PyTorch Geometric Data object.

        Requires torch_geometric to be installed.

        Args:
            graph: Graph dictionary from build_graph()

        Returns:
            PyTorch Geometric Data object
        """
        try:
            import torch
            from torch_geometric.data import Data
        except ImportError:
            raise ImportError(
                "torch_geometric required for to_pyg_data(). "
                "Install with: pip install torch-geometric"
            )

        return Data(
            x=torch.from_numpy(graph["x"]),
            edge_index=torch.from_numpy(graph["edge_index"]),
            pos=torch.from_numpy(graph["pos"]),
            grn_str=graph["grn_str"],
            grn_positions=graph["grn_positions"],
            protein_id=graph["protein_id"],
            protein_family=graph["protein_family"],
        )

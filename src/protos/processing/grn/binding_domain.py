"""BindingDomainConfig: manages binding domain configurations for graph construction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from protos.processing.grn.grn_utils import normalize_grn_format


class BindingDomainConfig:
    """
    Manages binding domain configurations that define graph edges and required GRNs.

    A binding domain configuration specifies which GRN positions are connected
    (edges) for each protein family. This is used for building graphs where
    nodes are GRN positions and edges represent structural/functional relationships.

    Example config format:
    {
        "gpcr_a": [
            ["7.45", "7.46"],
            ["5.42", "5.43"],
            ...
        ],
        "microbial_opsins": [
            ["67.004", "67.005"],
            ...
        ]
    }
    """

    def __init__(
        self,
        config_path: Optional[str] = None,
        config_dict: Optional[Dict[str, List[List[str]]]] = None,
    ) -> None:
        """
        Initialize BindingDomainConfig.

        Args:
            config_path: Path to JSON configuration file
            config_dict: Direct configuration dictionary (alternative to file)
        """
        self._config: Dict[str, List[List[str]]] = {}
        self._required_grns_cache: Dict[str, Set[str]] = {}

        if config_path:
            self.load_from_file(config_path)
        elif config_dict:
            self._config = config_dict
            self._normalize_config()

    def load_from_file(self, config_path: str) -> None:
        """Load configuration from a JSON file."""
        path = Path(config_path)
        if not path.exists():
            raise FileNotFoundError(f"Binding domain config not found: {config_path}")

        with open(path, "r") as f:
            self._config = json.load(f)

        self._normalize_config()
        self._required_grns_cache.clear()

    def _normalize_config(self) -> None:
        """Normalize all GRN strings in the configuration."""
        normalized: Dict[str, List[List[str]]] = {}
        for family, edges in self._config.items():
            normalized_edges = []
            for edge in edges:
                if len(edge) >= 2:
                    normalized_edges.append([
                        normalize_grn_format(edge[0]),
                        normalize_grn_format(edge[1]),
                    ])
            normalized[family] = normalized_edges
        self._config = normalized

    def get_families(self) -> List[str]:
        """Get list of available protein families."""
        return list(self._config.keys())

    def get_edges(self, protein_family: str) -> List[Tuple[str, str]]:
        """
        Get edge definitions for a protein family.

        Args:
            protein_family: Name of the protein family (e.g., 'gpcr_a')

        Returns:
            List of (grn1, grn2) tuples defining edges
        """
        edges = self._config.get(protein_family, [])
        return [(e[0], e[1]) for e in edges if len(e) >= 2]

    def get_required_grns(self, protein_family: str) -> Set[str]:
        """
        Get the set of all GRN positions required for a protein family.

        This extracts all unique GRNs from the edge definitions.

        Args:
            protein_family: Name of the protein family

        Returns:
            Set of normalized GRN strings
        """
        if protein_family in self._required_grns_cache:
            return self._required_grns_cache[protein_family]

        required = set()
        for edge in self._config.get(protein_family, []):
            if len(edge) >= 2:
                required.add(edge[0])
                required.add(edge[1])

        self._required_grns_cache[protein_family] = required
        return required

    def build_edge_index(
        self,
        grns: List[str],
        protein_family: str,
        bidirectional: bool = True,
    ) -> Tuple[List[int], List[int]]:
        """
        Build edge index arrays from GRN list and binding domain configuration.

        Creates edges based on the binding domain configuration, only including
        edges where both endpoints are present in the provided GRN list.

        Args:
            grns: Ordered list of GRN strings (defines node indices)
            protein_family: Name of the protein family
            bidirectional: If True, add edges in both directions (default: True)

        Returns:
            Tuple of (source_nodes, destination_nodes) lists
        """
        # Build GRN to index mapping
        grn_to_idx = {grn: i for i, grn in enumerate(grns)}

        src_nodes: List[int] = []
        dst_nodes: List[int] = []

        for edge in self._config.get(protein_family, []):
            if len(edge) < 2:
                continue

            src_grn = edge[0]
            dst_grn = edge[1]

            if src_grn in grn_to_idx and dst_grn in grn_to_idx:
                src_idx = grn_to_idx[src_grn]
                dst_idx = grn_to_idx[dst_grn]
                # Add edge in forward direction
                src_nodes.append(src_idx)
                dst_nodes.append(dst_idx)
                # Add edge in reverse direction for bidirectional graphs
                if bidirectional:
                    src_nodes.append(dst_idx)
                    dst_nodes.append(src_idx)

        return src_nodes, dst_nodes

    def filter_grns(
        self,
        grns: List[str],
        protein_family: str,
    ) -> List[str]:
        """
        Filter a list of GRNs to only those required by the binding domain.

        Args:
            grns: List of GRN strings to filter
            protein_family: Name of the protein family

        Returns:
            Filtered list maintaining original order
        """
        required = self.get_required_grns(protein_family)
        return [grn for grn in grns if grn in required]

    def get_edge_count(self, protein_family: str) -> int:
        """Get number of edges defined for a protein family."""
        return len(self._config.get(protein_family, []))

    def to_dict(self) -> Dict[str, List[List[str]]]:
        """Return configuration as dictionary."""
        return self._config.copy()

    def __repr__(self) -> str:
        families = ", ".join(self.get_families())
        return f"BindingDomainConfig(families=[{families}])"

"""Format schemas for models inputs and outputs.

This module defines the exact structure and requirements for each
input and output format used by the models system.
"""

from typing import Dict, List, Optional, Union, Any
from dataclasses import dataclass, field
from enum import Enum
import numpy as np
import pandas as pd
from pathlib import Path


# Format type definitions with detailed schemas

@dataclass
class SequenceFormat:
    """
    Protein sequence input format.
    
    Attributes:
        sequence: Amino acid sequence using standard one-letter codes
        sequence_id: Optional identifier for the sequence
        metadata: Optional metadata (organism, function, etc.)
    
    Examples:
        - Simple: "MKTAYIAKQRQISFVKSHFSRQ..."
        - With metadata: {"sequence": "MKTAY...", "sequence_id": "BACR_HALSA"}
    """
    sequence: str
    sequence_id: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate sequence format."""
        # Check if sequence contains only valid amino acids
        valid_aas = set("ACDEFGHIKLMNPQRSTVWYX")
        return bool(self.sequence) and all(aa in valid_aas for aa in self.sequence.upper())
    
    @classmethod
    def from_fasta(cls, fasta_content: str) -> List['SequenceFormat']:
        """Parse FASTA format."""
        sequences = []
        current_id = None
        current_seq = []
        
        for line in fasta_content.strip().split('\n'):
            if line.startswith('>'):
                if current_id:
                    sequences.append(cls(
                        sequence=''.join(current_seq),
                        sequence_id=current_id
                    ))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_id:
            sequences.append(cls(
                sequence=''.join(current_seq),
                sequence_id=current_id
            ))
        
        return sequences


@dataclass
class StructureFormat:
    """
    Protein structure input format.
    
    Attributes:
        coordinates: DataFrame with atomic coordinates
        pdb_id: PDB identifier
        chains: Available chain IDs
        metadata: Resolution, method, etc.
    
    Required DataFrame columns:
        - atom_name: Atom identifier (CA, CB, etc.)
        - auth_comp_id: Residue name (ALA, GLY, etc.)
        - auth_seq_id: Residue number
        - auth_chain_id: Chain identifier
        - x, y, z: Atomic coordinates
    
    Optional columns:
        - b_factor: Temperature factor
        - occupancy: Atom occupancy
        - element: Element symbol
    """
    coordinates: pd.DataFrame
    pdb_id: str
    chains: List[str]
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate structure format."""
        required_cols = {'atom_name', 'auth_comp_id', 'auth_seq_id', 
                        'auth_chain_id', 'x', 'y', 'z'}
        return bool(len(self.coordinates)) and required_cols.issubset(self.coordinates.columns)
    
    def get_ca_atoms(self, chain: Optional[str] = None) -> pd.DataFrame:
        """Extract CA atoms for a chain."""
        ca_atoms = self.coordinates[self.coordinates['atom_name'] == 'CA']
        if chain:
            ca_atoms = ca_atoms[ca_atoms['auth_chain_id'] == chain]
        return ca_atoms
    
    def to_graph_coordinates(self) -> np.ndarray:
        """Convert to coordinate array for graph construction."""
        ca_atoms = self.get_ca_atoms()
        return ca_atoms[['x', 'y', 'z']].values.astype(np.float32)


@dataclass 
class EmbeddingFormat:
    """
    Protein embedding input/output format.
    
    Attributes:
        embeddings: Embedding matrix [seq_len, embed_dim] or [embed_dim]
        embedding_type: Type of embedding (per_residue, mean, cls)
        model_name: Model that generated embeddings
        layer: Which layer the embeddings are from
        metadata: Additional information
    
    Examples:
        - Per-residue: shape [100, 1280] for 100-residue protein
        - Mean pooled: shape [1280]
        - CLS token: shape [1280]
    """
    embeddings: np.ndarray
    embedding_type: str = "per_residue"
    model_name: Optional[str] = None
    layer: Optional[int] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate embedding format."""
        if self.embedding_type == "per_residue":
            return len(self.embeddings.shape) == 2
        else:  # mean or cls
            return len(self.embeddings.shape) == 1
    
    def get_pooled(self, method: str = "mean") -> np.ndarray:
        """Get pooled representation."""
        if len(self.embeddings.shape) == 1:
            return self.embeddings
        
        if method == "mean":
            return np.mean(self.embeddings, axis=0)
        elif method == "max":
            return np.max(self.embeddings, axis=0)
        else:
            raise ValueError(f"Unknown pooling method: {method}")


@dataclass
class MSAFormat:
    """
    Multiple Sequence Alignment input format.
    
    Attributes:
        alignment: Aligned sequences [n_seqs, seq_len]
        sequence_ids: Identifiers for each sequence
        target_seq_idx: Index of target sequence
        weights: Optional sequence weights
        metadata: Coverage, diversity metrics, etc.
    """
    alignment: List[str]
    sequence_ids: List[str]
    target_seq_idx: int = 0
    weights: Optional[np.ndarray] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate MSA format."""
        if not self.alignment:
            return False
        
        # Check all sequences have same length
        seq_lens = [len(seq) for seq in self.alignment]
        return len(set(seq_lens)) == 1
    
    def get_depth(self) -> int:
        """Get MSA depth."""
        return len(self.alignment)
    
    def get_conservation(self) -> np.ndarray:
        """Calculate per-position conservation."""
        # Simple conservation score
        n_seqs = len(self.alignment)
        seq_len = len(self.alignment[0])
        conservation = np.zeros(seq_len)
        
        for i in range(seq_len):
            col = [seq[i] for seq in self.alignment]
            # Most common residue frequency
            from collections import Counter
            counts = Counter(col)
            conservation[i] = counts.most_common(1)[0][1] / n_seqs
        
        return conservation


@dataclass
class GraphFormat:
    """
    Graph representation format (for PyTorch Geometric).
    
    Attributes:
        node_features: Node feature matrix [n_nodes, n_features]
        edge_index: Edge connectivity [2, n_edges]
        edge_features: Optional edge features [n_edges, n_edge_features]
        pos: Optional node positions [n_nodes, 3]
        batch: Optional batch assignment
        metadata: Graph construction parameters
    """
    node_features: np.ndarray
    edge_index: np.ndarray
    edge_features: Optional[np.ndarray] = None
    pos: Optional[np.ndarray] = None
    batch: Optional[np.ndarray] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate graph format."""
        n_nodes = self.node_features.shape[0]
        
        # Check edge index
        if self.edge_index.shape[0] != 2:
            return False
        if self.edge_index.max() >= n_nodes:
            return False
        
        # Check optional arrays
        if self.pos is not None and self.pos.shape[0] != n_nodes:
            return False
        
        return True
    
    def to_data(self):
        """Convert to PyTorch Geometric Data object."""
        try:
            import torch
            from torch_geometric.data import Data
            
            data = Data(
                x=torch.tensor(self.node_features, dtype=torch.float32),
                edge_index=torch.tensor(self.edge_index, dtype=torch.long)
            )
            
            if self.edge_features is not None:
                data.edge_attr = torch.tensor(self.edge_features, dtype=torch.float32)
            if self.pos is not None:
                data.pos = torch.tensor(self.pos, dtype=torch.float32)
            if self.batch is not None:
                data.batch = torch.tensor(self.batch, dtype=torch.long)
            
            return data
        except ImportError:
            raise ImportError("PyTorch Geometric required for graph conversion")


@dataclass
class GRNFormat:
    """
    Generic Residue Numbering format.
    
    Attributes:
        grn_positions: GRN position assignments
        sequence: Original sequence
        grn_system: Numbering system used
        metadata: Family-specific information
    
    Examples:
        - Series: positions indexed by GRN (1.50, 2.50, etc.)
        - DataFrame: multiple proteins with GRN columns
    """
    grn_positions: Union[pd.Series, pd.DataFrame]
    sequence: Optional[str] = None
    grn_system: str = "ballesteros_weinstein"
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate GRN format."""
        if isinstance(self.grn_positions, pd.Series):
            # Check index format (e.g., "1.50", "2.50")
            return all('.' in str(idx) for idx in self.grn_positions.index)
        else:
            # DataFrame should have GRN columns
            return any('.' in col for col in self.grn_positions.columns)
    
    def get_conserved_positions(self, threshold: float = 0.9) -> List[str]:
        """Get highly conserved GRN positions."""
        if isinstance(self.grn_positions, pd.DataFrame):
            # Count non-null values per position
            conservation = self.grn_positions.notna().mean()
            return list(conservation[conservation > threshold].index)
        else:
            # Single sequence - return all positions
            return list(self.grn_positions.index)


@dataclass
class PropertyFormat:
    """
    Property/label format for predictions.
    
    Attributes:
        properties: Property values (dict or single value)
        property_names: Names of properties
        property_types: Types (regression, classification, etc.)
        confidence: Optional confidence scores
        metadata: Source, experimental conditions, etc.
    
    Examples:
        - Single: {"lambda_max": 470}
        - Multiple: {"lambda_max": 470, "photocycle": "fast", "pump_type": "proton"}
    """
    properties: Union[Dict[str, Any], Any]
    property_names: Optional[List[str]] = None
    property_types: Optional[Dict[str, str]] = None
    confidence: Optional[Dict[str, float]] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def validate(self) -> bool:
        """Validate property format."""
        if isinstance(self.properties, dict):
            return len(self.properties) > 0
        else:
            return self.properties is not None
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame format."""
        if isinstance(self.properties, dict):
            return pd.DataFrame([self.properties])
        else:
            return pd.DataFrame([{"property": self.properties}])


# Output format schemas

@dataclass
class StructureOutput:
    """
    Structure prediction output format.
    
    Attributes:
        coordinates: Predicted atomic coordinates
        plddt: Per-residue confidence scores
        pae: Predicted aligned error (optional)
        ranking_confidence: Overall structure confidence
        metadata: Prediction parameters
    """
    coordinates: pd.DataFrame
    plddt: np.ndarray
    pae: Optional[np.ndarray] = None
    ranking_confidence: Optional[float] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def to_pdb(self, output_path: Path):
        """Save as PDB file."""
        # PDB writing logic
        pass
    
    def to_cif(self, output_path: Path):
        """Save as mmCIF file."""
        # CIF writing logic
        pass


@dataclass
class AttentionOutput:
    """
    Attention weights output format.
    
    Attributes:
        attention_weights: Attention matrices [layers, heads, seq_len, seq_len]
        layer_names: Names of layers
        aggregated: Pre-computed aggregations (mean, max)
    """
    attention_weights: np.ndarray
    layer_names: List[str]
    aggregated: Optional[Dict[str, np.ndarray]] = None
    
    def get_layer(self, layer: Union[int, str]) -> np.ndarray:
        """Get attention for specific layer."""
        if isinstance(layer, int):
            return self.attention_weights[layer]
        else:
            idx = self.layer_names.index(layer)
            return self.attention_weights[idx]


@dataclass
class LogitsOutput:
    """
    Raw models logits output format.
    
    Attributes:
        logits: Raw models outputs
        probabilities: Softmax probabilities (if applicable)
        predictions: Argmax predictions (if applicable)
        class_names: Class labels for classification
    """
    logits: np.ndarray
    probabilities: Optional[np.ndarray] = None
    predictions: Optional[np.ndarray] = None
    class_names: Optional[List[str]] = None


# Format conversion utilities

class FormatConverter:
    """Utilities for converting between formats."""
    
    @staticmethod
    def sequence_to_embedding_input(sequence: SequenceFormat) -> Dict[str, Any]:
        """Convert sequence to embedding models input."""
        return {
            "sequence": sequence.sequence,
            "sequence_id": sequence.sequence_id
        }
    
    @staticmethod
    def structure_to_sequence(structure: StructureFormat, chain: str = 'A') -> SequenceFormat:
        """Extract sequence from structure."""
        chain_data = structure.coordinates[
            structure.coordinates['auth_chain_id'] == chain
        ]
        
        # Get CA atoms in order
        ca_atoms = chain_data[chain_data['atom_name'] == 'CA'].sort_values('auth_seq_id')
        
        # Convert 3-letter to 1-letter codes
        aa_3to1 = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        
        sequence = ''.join([
            aa_3to1.get(res, 'X') 
            for res in ca_atoms['auth_comp_id']
        ])
        
        return SequenceFormat(
            sequence=sequence,
            sequence_id=f"{structure.pdb_id}_{chain}"
        )
    
    @staticmethod
    def structure_to_graph(structure: StructureFormat, 
                          method: str = "knn",
                          k: int = 10,
                          threshold: float = 8.0) -> GraphFormat:
        """Convert structure to graph representation."""
        # Get CA coordinates
        coords = structure.to_graph_coordinates()
        n_residues = len(coords)
        
        # Build edges
        if method == "knn":
            from sklearn.neighbors import NearestNeighbors
            nbrs = NearestNeighbors(n_neighbors=k+1)
            nbrs.fit(coords)
            distances, indices = nbrs.kneighbors(coords)
            
            edges = []
            for i in range(n_residues):
                for j in range(1, k+1):  # Skip self
                    if distances[i, j] < threshold:
                        edges.append([i, indices[i, j]])
            
            edge_index = np.array(edges).T
        else:
            # Distance threshold method
            from scipy.spatial.distance import cdist
            dist_matrix = cdist(coords, coords)
            edge_index = np.array(np.where(
                (dist_matrix < threshold) & (dist_matrix > 0)
            ))
        
        # Simple node features (can be extended)
        node_features = np.eye(20)[:n_residues]  # Dummy features
        
        return GraphFormat(
            node_features=node_features,
            edge_index=edge_index,
            pos=coords
        )


# Validation functions

def validate_model_input(input_data: Any, expected_format: str) -> bool:
    """Validate input data matches expected format."""
    format_validators = {
        "sequence": lambda x: isinstance(x, (str, SequenceFormat)),
        "structure": lambda x: isinstance(x, (pd.DataFrame, StructureFormat)),
        "embedding": lambda x: isinstance(x, (np.ndarray, EmbeddingFormat)),
        "msa": lambda x: isinstance(x, (list, MSAFormat)),
        "graph": lambda x: hasattr(x, 'node_features') or isinstance(x, GraphFormat),
        "grn": lambda x: isinstance(x, (pd.Series, pd.DataFrame, GRNFormat)),
        "property": lambda x: isinstance(x, (dict, PropertyFormat))
    }
    
    validator = format_validators.get(expected_format)
    if validator:
        return validator(input_data)
    return False


def validate_model_output(output_data: Any, expected_format: str) -> bool:
    """Validate output data matches expected format."""
    format_validators = {
        "embedding": lambda x: isinstance(x, (np.ndarray, dict)),
        "structure": lambda x: isinstance(x, (pd.DataFrame, StructureOutput)),
        "property": lambda x: isinstance(x, (dict, float, PropertyFormat)),
        "attention": lambda x: isinstance(x, (np.ndarray, AttentionOutput)),
        "logits": lambda x: isinstance(x, (np.ndarray, LogitsOutput)),
        "graph": lambda x: hasattr(x, 'node_features') or isinstance(x, GraphFormat)
    }
    
    validator = format_validators.get(expected_format)
    if validator:
        return validator(output_data)
    return False
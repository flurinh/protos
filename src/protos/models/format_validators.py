"""Format validation utilities for models inputs and outputs.

This module provides validation and conversion utilities to ensure
data compatibility between Protos processors and AI models.
"""

from typing import Any, Dict, List, Optional, Union, Tuple
import numpy as np
import pandas as pd
from pathlib import Path
import logging

from protos.models.format_schemas import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    MSAFormat, GraphFormat, GRNFormat, PropertyFormat,
    StructureOutput, AttentionOutput, LogitsOutput,
    FormatConverter
)

logger = logging.getLogger(__name__)


class FormatValidator:
    """Validates and prepares data formats for models."""
    
    def __init__(self):
        """Initialize validator with format mappings."""
        self.format_classes = {
            "sequence": SequenceFormat,
            "structure": StructureFormat,
            "embedding": EmbeddingFormat,
            "msa": MSAFormat,
            "graph": GraphFormat,
            "grn": GRNFormat,
            "property": PropertyFormat
        }
        
        self.converter = FormatConverter()
    
    def validate_input(self, data: Any, format_type: str) -> Tuple[bool, Optional[str]]:
        """
        Validate input data for a specific format.
        
        Args:
            data: Input data to validate
            format_type: Expected format type
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            if format_type == "sequence":
                return self._validate_sequence(data)
            elif format_type == "structure":
                return self._validate_structure(data)
            elif format_type == "embedding":
                return self._validate_embedding(data)
            elif format_type == "msa":
                return self._validate_msa(data)
            elif format_type == "graph":
                return self._validate_graph(data)
            elif format_type == "grn":
                return self._validate_grn(data)
            elif format_type == "property":
                return self._validate_property(data)
            else:
                return False, f"Unknown format type: {format_type}"
        except Exception as e:
            return False, f"Validation error: {str(e)}"
    
    def _validate_sequence(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate sequence format."""
        if isinstance(data, str):
            # Simple string sequence
            valid_aas = set("ACDEFGHIKLMNPQRSTVWYX")
            if all(aa in valid_aas for aa in data.upper()):
                return True, None
            else:
                return False, "Invalid amino acid characters in sequence"
        
        elif isinstance(data, SequenceFormat) or (hasattr(data, '__class__') and data.__class__.__name__ == 'SequenceFormat'):
            if hasattr(data, 'validate') and data.validate():
                return True, None
            else:
                return False, "Invalid SequenceFormat object"
        
        elif isinstance(data, dict):
            # Dictionary with sequence key
            if 'sequence' in data:
                return self._validate_sequence(data['sequence'])
            else:
                return False, "Dictionary missing 'sequence' key"
        
        else:
            return False, f"Invalid sequence format: {type(data)}"
    
    def _validate_structure(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate structure format."""
        if isinstance(data, pd.DataFrame):
            required_cols = {'atom_name', 'auth_comp_id', 'auth_seq_id', 
                           'auth_chain_id', 'x', 'y', 'z'}
            if required_cols.issubset(data.columns):
                return True, None
            else:
                missing = required_cols - set(data.columns)
                return False, f"Missing required columns: {missing}"
        
        elif isinstance(data, StructureFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid StructureFormat object"
        
        else:
            return False, f"Invalid structure format: {type(data)}"
    
    def _validate_embedding(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate embedding format."""
        if isinstance(data, np.ndarray):
            if len(data.shape) in [1, 2]:
                return True, None
            else:
                return False, f"Invalid embedding shape: {data.shape}"
        
        elif isinstance(data, EmbeddingFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid EmbeddingFormat object"
        
        else:
            return False, f"Invalid embedding format: {type(data)}"
    
    def _validate_msa(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate MSA format."""
        if isinstance(data, list):
            if all(isinstance(seq, str) for seq in data):
                # Check all sequences have same length
                lengths = [len(seq) for seq in data]
                if len(set(lengths)) == 1:
                    return True, None
                else:
                    return False, "MSA sequences have different lengths"
            else:
                return False, "MSA must be list of strings"
        
        elif isinstance(data, MSAFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid MSAFormat object"
        
        else:
            return False, f"Invalid MSA format: {type(data)}"
    
    def _validate_graph(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate graph format."""
        if hasattr(data, 'x') and hasattr(data, 'edge_index'):
            # PyTorch Geometric Data object
            return True, None
        
        elif isinstance(data, GraphFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid GraphFormat object"
        
        elif isinstance(data, dict):
            if 'node_features' in data and 'edge_index' in data:
                return True, None
            else:
                return False, "Graph dict missing required keys"
        
        else:
            return False, f"Invalid graph format: {type(data)}"
    
    def _validate_grn(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate GRN format."""
        if isinstance(data, pd.Series):
            # Check index contains GRN positions
            if any('.' in str(idx) for idx in data.index):
                return True, None
            else:
                return False, "GRN Series missing position format (X.YY)"
        
        elif isinstance(data, pd.DataFrame):
            # Check columns contain GRN positions
            if any('.' in col for col in data.columns):
                return True, None
            else:
                return False, "GRN DataFrame missing position columns"
        
        elif isinstance(data, GRNFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid GRNFormat object"
        
        else:
            return False, f"Invalid GRN format: {type(data)}"
    
    def _validate_property(self, data: Any) -> Tuple[bool, Optional[str]]:
        """Validate property format."""
        if isinstance(data, (dict, float, int, str)):
            return True, None
        
        elif isinstance(data, PropertyFormat):
            if data.validate():
                return True, None
            else:
                return False, "Invalid PropertyFormat object"
        
        else:
            return False, f"Invalid property format: {type(data)}"
    
    def prepare_input(self, data: Dict[str, Any], 
                     required_formats: List[str]) -> Dict[str, Any]:
        """
        Prepare and validate input data for models.
        
        Args:
            data: Raw input data from Protos processors
            required_formats: List of required input formats
            
        Returns:
            Prepared input data
            
        Raises:
            ValueError: If validation fails
        """
        prepared = {}
        
        for format_type in required_formats:
            if format_type not in data:
                raise ValueError(f"Missing required input format: {format_type}")
            
            # Validate
            is_valid, error = self.validate_input(data[format_type], format_type)
            if not is_valid:
                raise ValueError(f"Invalid {format_type}: {error}")
            
            # Convert to standard format if needed
            prepared[format_type] = self._standardize_format(
                data[format_type], 
                format_type
            )
        
        return prepared
    
    def _standardize_format(self, data: Any, format_type: str) -> Any:
        """Convert data to standard format object."""
        format_class = self.format_classes.get(format_type)
        
        if format_class and not isinstance(data, format_class):
            # Convert to standard format
            if format_type == "sequence" and isinstance(data, str):
                return SequenceFormat(sequence=data)
            
            elif format_type == "embedding" and isinstance(data, np.ndarray):
                return EmbeddingFormat(
                    embeddings=data,
                    embedding_type="per_residue" if len(data.shape) == 2 else "pooled"
                )
            
            # Add more conversions as needed
        
        return data
    
    def validate_output(self, data: Any, format_type: str) -> Tuple[bool, Optional[str]]:
        """
        Validate models output data.
        
        Args:
            data: Output data from models
            format_type: Expected output format
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            if format_type == "embedding":
                if isinstance(data, (np.ndarray, dict)):
                    return True, None
                else:
                    return False, f"Invalid embedding output: {type(data)}"
            
            elif format_type == "structure":
                if isinstance(data, (pd.DataFrame, StructureOutput)):
                    return True, None
                else:
                    return False, f"Invalid structure output: {type(data)}"
            
            elif format_type == "property":
                if isinstance(data, (dict, float, int, str, PropertyFormat)):
                    return True, None
                else:
                    return False, f"Invalid property output: {type(data)}"
            
            elif format_type == "attention":
                if isinstance(data, (np.ndarray, AttentionOutput)):
                    return True, None
                else:
                    return False, f"Invalid attention output: {type(data)}"
            
            elif format_type == "logits":
                if isinstance(data, (np.ndarray, LogitsOutput)):
                    return True, None
                else:
                    return False, f"Invalid logits output: {type(data)}"
            
            else:
                return False, f"Unknown output format: {format_type}"
                
        except Exception as e:
            return False, f"Output validation error: {str(e)}"


class FormatAdapter:
    """Adapts between Protos and models-specific formats."""
    
    def __init__(self):
        """Initialize adapter with converters."""
        self.converter = FormatConverter()
        self.validator = FormatValidator()
    
    def adapt_for_model(self, data: Dict[str, Any], model_name: str, 
                       model_config: Dict[str, Any]) -> Any:
        """
        Adapt Protos data for specific models input.
        
        Args:
            data: Validated input data
            model_name: Name of the models
            model_config: Model configuration
            
        Returns:
            Model-specific input format
        """
        if model_name == "esm2":
            return self._adapt_for_esm2(data, model_config)
        elif model_name == "ankh":
            return self._adapt_for_ankh(data, model_config)
        elif model_name == "lambda":
            return self._adapt_for_lambda(data, model_config)
        elif model_name == "boltz1":
            return self._adapt_for_boltz(data, model_config)
        else:
            # Generic adaptation
            return data
    
    def _adapt_for_esm2(self, data: Dict[str, Any], config: Dict[str, Any]) -> Any:
        """Adapt for ESM-2 input format."""
        if 'sequence' in data:
            seq_data = data['sequence']
            if isinstance(seq_data, SequenceFormat):
                return [(seq_data.sequence_id or "protein", seq_data.sequence)]
            elif isinstance(seq_data, str):
                return [("protein", seq_data)]
        
        raise ValueError("ESM-2 requires sequence input")
    
    def _adapt_for_ankh(self, data: Dict[str, Any], config: Dict[str, Any]) -> Any:
        """Adapt for Ankh input format."""
        if 'sequence' in data:
            seq_data = data['sequence']
            if isinstance(seq_data, SequenceFormat):
                return seq_data.sequence
            elif isinstance(seq_data, str):
                return seq_data
        
        raise ValueError("Ankh requires sequence input")
    
    def _adapt_for_lambda(self, data: Dict[str, Any], config: Dict[str, Any]) -> Any:
        """Adapt for Lambda graph models format."""
        # Lambda can use multiple inputs
        if 'graph' in data:
            # Already in graph format
            graph_data = data['graph']
            if isinstance(graph_data, GraphFormat):
                return graph_data.to_data()
            else:
                return graph_data
        
        elif 'structure' in data:
            # Convert structure to graph
            struct_data = data['structure']
            if isinstance(struct_data, StructureFormat):
                graph = self.converter.structure_to_graph(
                    struct_data,
                    **config.get('preprocessing_params', {})
                )
                return graph.to_data()
        
        raise ValueError("Lambda requires graph or structure input")
    
    def _adapt_for_boltz(self, data: Dict[str, Any], config: Dict[str, Any]) -> Any:
        """Adapt for Boltz input format."""
        inputs = {}
        
        if 'sequence' in data:
            seq_data = data['sequence']
            if isinstance(seq_data, SequenceFormat):
                inputs['sequence'] = seq_data.sequence
            else:
                inputs['sequence'] = seq_data
        
        if 'msa' in data:
            msa_data = data['msa']
            if isinstance(msa_data, MSAFormat):
                inputs['msa'] = msa_data.alignment
            else:
                inputs['msa'] = msa_data
        
        return inputs
    
    def adapt_from_model(self, output: Any, output_format: str, 
                        model_name: str) -> Any:
        """
        Adapt models output to Protos format.
        
        Args:
            output: Raw models output
            output_format: Expected output format
            model_name: Name of the models
            
        Returns:
            Protos-compatible output format
        """
        # Validate output
        is_valid, error = self.validator.validate_output(output, output_format)
        if not is_valid:
            logger.warning(f"Invalid models output: {error}")
        
        # Model-specific adaptations
        if model_name == "esm2" and output_format == "embedding":
            # ESM-2 returns dict with multiple outputs
            if isinstance(output, dict) and 'embeddings' in output:
                return EmbeddingFormat(
                    embeddings=output['embeddings'],
                    embedding_type="per_residue",
                    model_name="esm2"
                )
        
        elif model_name == "boltz1" and output_format == "structure":
            # Boltz returns coordinates and confidence
            if isinstance(output, dict):
                return StructureOutput(
                    coordinates=output.get('coordinates'),
                    plddt=output.get('plddt'),
                    pae=output.get('pae'),
                    ranking_confidence=output.get('ranking_confidence')
                )
        
        # Default: return as-is
        return output


# Convenience functions

def validate_model_compatibility(entity_formats: List[str], 
                               model_inputs: List[str]) -> bool:
    """
    Check if entity has required formats for models.
    
    Args:
        entity_formats: Available formats for entity
        model_inputs: Required input formats for models
        
    Returns:
        True if compatible
    """
    return set(model_inputs).issubset(set(entity_formats))


def suggest_format_conversions(available: List[str], 
                             required: List[str]) -> List[Tuple[str, str]]:
    """
    Suggest possible format conversions.
    
    Args:
        available: Available formats
        required: Required formats
        
    Returns:
        List of (from_format, to_format) conversions
    """
    conversions = []
    
    # Define possible conversions
    conversion_map = {
        ('structure', 'sequence'): "Extract sequence from structure",
        ('structure', 'graph'): "Convert structure to graph",
        ('sequence', 'embedding'): "Generate embeddings from sequence",
        ('embedding', 'property'): "Predict properties from embeddings",
    }
    
    for required_fmt in required:
        if required_fmt not in available:
            # Find possible conversions
            for (from_fmt, to_fmt), desc in conversion_map.items():
                if from_fmt in available and to_fmt == required_fmt:
                    conversions.append((from_fmt, to_fmt))
    
    return conversions
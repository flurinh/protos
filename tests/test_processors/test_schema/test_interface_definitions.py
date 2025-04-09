"""
Tests for the interface_definitions module.

These tests verify that the standardized interface definitions work correctly
and provide consistent method signatures for operations on different data types.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch

from protos.processing.schema.interface_definitions import (
    # Interface classes
    StructureInterface,
    GRNInterface,
    SequenceInterface,
    
    # Return value formats
    AlignmentResult,
    DistanceResult,
    GRNMappingResult,
    FeatureResult,
    
    # Validation decorators
    validate_structure_operation,
    validate_grn_operation,
    validate_sequence_operation
)

class TestStructureInterface:
    """Tests for the StructureInterface abstract class."""
    
    def test_structure_interface_abstract(self):
        """Test that StructureInterface is an abstract class."""
        # Attempting to instantiate should raise TypeError
        with pytest.raises(TypeError):
            interface = StructureInterface()
    
    def test_structure_interface_methods(self):
        """Test that StructureInterface defines required methods."""
        # Check that the expected methods are defined
        assert hasattr(StructureInterface, 'load_structure')
        assert hasattr(StructureInterface, 'save_structure')
        assert hasattr(StructureInterface, 'get_atoms')
        assert hasattr(StructureInterface, 'get_residues')
        assert hasattr(StructureInterface, 'get_chains')
        assert hasattr(StructureInterface, 'get_sequence')
        assert hasattr(StructureInterface, 'calculate_distance')
    
    def test_structure_interface_implementation(self):
        """Test implementing the StructureInterface with a concrete class."""
        # Create a concrete implementation
        class ConcreteStructureProcessor(StructureInterface):
            @staticmethod
            def load_structure(file_path, structure_id=None):
                return pd.DataFrame({'test': [1, 2, 3]})
                
            @staticmethod
            def save_structure(structure_df, file_path):
                pass
                
            @staticmethod
            def get_atoms(structure_df, atom_names=None, chain_id=None, residue_ids=None):
                return pd.DataFrame({'test': [1, 2, 3]})
                
            @staticmethod
            def get_residues(structure_df, chain_id=None, residue_ids=None):
                return pd.DataFrame({'test': [1, 2, 3]})
                
            @staticmethod
            def get_chains(structure_df):
                return ['A', 'B', 'C']
                
            @staticmethod
            def get_sequence(structure_df, chain_id):
                return 'ACDEFGHI'
                
            @staticmethod
            def calculate_distance(structure_df, atom1_index, atom2_index):
                return 3.8
        
        # Instantiate the concrete implementation
        processor = ConcreteStructureProcessor()
        
        # Test that the methods can be called
        assert isinstance(processor.load_structure('test.cif'), pd.DataFrame)
        assert processor.get_chains(pd.DataFrame()) == ['A', 'B', 'C']
        assert processor.get_sequence(pd.DataFrame(), 'A') == 'ACDEFGHI'
        assert processor.calculate_distance(pd.DataFrame(), 1, 2) == 3.8


class TestGRNInterface:
    """Tests for the GRNInterface abstract class."""
    
    def test_grn_interface_abstract(self):
        """Test that GRNInterface is an abstract class."""
        # Attempting to instantiate should raise TypeError
        with pytest.raises(TypeError):
            interface = GRNInterface()
    
    def test_grn_interface_methods(self):
        """Test that GRNInterface defines required methods."""
        # Check that the expected methods are defined
        assert hasattr(GRNInterface, 'load_grn_table')
        assert hasattr(GRNInterface, 'save_grn_table')
        assert hasattr(GRNInterface, 'assign_grns')
        assert hasattr(GRNInterface, 'map_grns_to_structure')
        assert hasattr(GRNInterface, 'get_residue_by_grn')
    
    def test_grn_interface_implementation(self):
        """Test implementing the GRNInterface with a concrete class."""
        # Create a concrete implementation
        class ConcreteGRNProcessor(GRNInterface):
            @staticmethod
            def load_grn_table(file_path):
                return pd.DataFrame({'test': [1, 2, 3]})
                
            @staticmethod
            def save_grn_table(grn_table, file_path):
                pass
                
            @staticmethod
            def assign_grns(sequence, reference_id, grn_table):
                return {1: '1x50', 2: '2x50', 3: '3x50'}
                
            @staticmethod
            def map_grns_to_structure(structure_df, grn_mapping, chain_id):
                return pd.DataFrame({'test': [1, 2, 3]})
                
            @staticmethod
            def get_residue_by_grn(structure_df, grn, chain_id):
                return pd.DataFrame({'test': [1]})
        
        # Instantiate the concrete implementation
        processor = ConcreteGRNProcessor()
        
        # Test that the methods can be called
        assert isinstance(processor.load_grn_table('test.csv'), pd.DataFrame)
        assert processor.assign_grns('ACDEFG', 'ref1', pd.DataFrame()) == {1: '1x50', 2: '2x50', 3: '3x50'}
        assert isinstance(processor.get_residue_by_grn(pd.DataFrame(), '1x50', 'A'), pd.DataFrame)


class TestSequenceInterface:
    """Tests for the SequenceInterface abstract class."""
    
    def test_sequence_interface_abstract(self):
        """Test that SequenceInterface is an abstract class."""
        # Attempting to instantiate should raise TypeError
        with pytest.raises(TypeError):
            interface = SequenceInterface()
    
    def test_sequence_interface_methods(self):
        """Test that SequenceInterface defines required methods."""
        # Check that the expected methods are defined
        assert hasattr(SequenceInterface, 'align_sequences')
        assert hasattr(SequenceInterface, 'load_fasta')
        assert hasattr(SequenceInterface, 'save_fasta')
        assert hasattr(SequenceInterface, 'calculate_identity')
        assert hasattr(SequenceInterface, 'validate_sequence')
    
    def test_sequence_interface_implementation(self):
        """Test implementing the SequenceInterface with a concrete class."""
        # Create a concrete implementation
        class ConcreteSequenceProcessor(SequenceInterface):
            @staticmethod
            def align_sequences(query_sequence, target_sequence, gap_open=-10.0, gap_extend=-0.5):
                return ('ACDEF', 'AC-EF', 0.8)
                
            @staticmethod
            def load_fasta(file_path):
                return {'seq1': 'ACDEF', 'seq2': 'ACGHI'}
                
            @staticmethod
            def save_fasta(sequences, file_path, width=80):
                pass
                
            @staticmethod
            def calculate_identity(sequence1, sequence2):
                return 80.0
                
            @staticmethod
            def validate_sequence(sequence):
                return True
        
        # Instantiate the concrete implementation
        processor = ConcreteSequenceProcessor()
        
        # Test that the methods can be called
        assert processor.align_sequences('ACDEF', 'ACGHI') == ('ACDEF', 'AC-EF', 0.8)
        assert processor.load_fasta('test.fasta') == {'seq1': 'ACDEF', 'seq2': 'ACGHI'}
        assert processor.calculate_identity('ACDEF', 'ACGHI') == 80.0
        assert processor.validate_sequence('ACDEF') is True


class TestValidationDecorators:
    """Tests for the validation decorator functions."""
    
    def test_validate_structure_operation(self):
        """Test the structure operation validation decorator."""
        # Create a mock for validate_structure_df
        with patch('protos.processing.schema.schema_definitions.validate_structure_df') as mock_validate:
            # Set up the mock to return True for validation
            mock_validate.return_value = True
            
            # Create a function to decorate
            @validate_structure_operation
            def test_function(df):
                return df
            
            # Create a test DataFrame
            test_df = pd.DataFrame({'test': [1, 2, 3]})
            
            # Call the decorated function
            result = test_function(test_df)
            
            # Check that validation was called
            mock_validate.assert_called_once_with(test_df)
            
            # Check that the function returned the DataFrame
            assert result is test_df
    
    def test_validate_grn_operation(self):
        """Test the GRN operation validation decorator."""
        # Create a function to decorate
        @validate_grn_operation
        def test_function(grn_mapping):
            return grn_mapping
        
        # Create a test GRN mapping
        test_mapping = {1: '1x50', 2: '2x50', 3: '3x50'}
        
        # Call the decorated function
        result = test_function(test_mapping)
        
        # Check that the function returned the mapping
        assert result is test_mapping
        
        # Test with invalid input (not a dictionary)
        with pytest.raises(ValueError):
            test_function("not a dictionary")
    
    def test_validate_sequence_operation(self):
        """Test the sequence operation validation decorator."""
        # Create a function to decorate
        @validate_sequence_operation
        def test_function(sequence):
            return sequence
        
        # Create a test sequence
        test_sequence = "ACDEFGHI"
        
        # Call the decorated function
        result = test_function(test_sequence)
        
        # Check that the function returned the sequence
        assert result is test_sequence
        
        # Test with invalid input (not a string)
        with pytest.raises(ValueError):
            test_function(123)  # Integer instead of string


class TestReturnValueFormats:
    """Tests for the standardized return value formats."""
    
    def test_alignment_result(self):
        """Test the AlignmentResult type annotation."""
        # Check that AlignmentResult is defined as a tuple with 3 elements
        from typing import get_type_hints, get_origin, get_args
        
        # AlignmentResult should be a tuple with 3 elements: str, str, float
        assert get_origin(AlignmentResult) is tuple
        args = get_args(AlignmentResult)
        assert len(args) == 3
        assert args[0] is str  # aligned_query
        assert args[1] is str  # aligned_target
        assert args[2] is float  # score
        
        # Create a conforming result
        result: AlignmentResult = ('ACDEF', 'AC-EF', 0.8)
        assert isinstance(result[0], str)
        assert isinstance(result[1], str)
        assert isinstance(result[2], float)
    
    def test_distance_result(self):
        """Test the DistanceResult type annotation."""
        # Check that DistanceResult is defined as a dictionary
        from typing import get_type_hints, get_origin, get_args
        
        # DistanceResult should be a dict with tuple keys and float values
        assert get_origin(DistanceResult) is dict
        args = get_args(DistanceResult)
        assert len(args) == 2
        
        key_type, value_type = args
        assert get_origin(key_type) is tuple
        assert get_args(key_type) == (int, int)
        assert value_type is float
        
        # Create a conforming result
        result: DistanceResult = {(1, 2): 3.8, (1, 3): 5.2}
        assert isinstance(result, dict)
        for key, value in result.items():
            assert isinstance(key, tuple)
            assert isinstance(key[0], int)
            assert isinstance(key[1], int)
            assert isinstance(value, float)
    
    def test_grn_mapping_result(self):
        """Test the GRNMappingResult type annotation."""
        # Check that GRNMappingResult is defined as a dictionary
        from typing import get_type_hints, get_origin, get_args
        
        # GRNMappingResult should be a dict with int keys and str values
        assert get_origin(GRNMappingResult) is dict
        args = get_args(GRNMappingResult)
        assert len(args) == 2
        assert args[0] is int  # residue position
        assert args[1] is str  # grn
        
        # Create a conforming result
        result: GRNMappingResult = {1: '1x50', 2: '2x50', 3: '3x50'}
        assert isinstance(result, dict)
        for key, value in result.items():
            assert isinstance(key, int)
            assert isinstance(value, str)
    
    def test_feature_result(self):
        """Test the FeatureResult type annotation."""
        # Check that FeatureResult is defined as a dictionary
        from typing import get_type_hints, get_origin, get_args, Any
        
        # FeatureResult should be a dict with str keys and Any values
        assert get_origin(FeatureResult) is dict
        args = get_args(FeatureResult)
        assert len(args) == 2
        assert args[0] is str  # feature name
        
        # Create a conforming result
        result: FeatureResult = {
            'sequence_length': 100,
            'hydrophobicity': 0.45,
            'secondary_structure': ['H', 'H', 'L', 'E', 'E']
        }
        assert isinstance(result, dict)
        for key in result:
            assert isinstance(key, str)
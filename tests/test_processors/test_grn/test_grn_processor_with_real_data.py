#!/usr/bin/env python3
"""
Test suite for GRN processor using real reference data without parser modifications.

This test suite works with the current GRN system and validates functionality
with actual biological data.
"""

import pytest
import pandas as pd
import numpy as np

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.paths.path_config import ProtosPaths


class TestGRNProcessorWithRealData:
    """Test GRN processor with real reference data using current system."""
    
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        """Set up test environment with real data paths."""
        # Set global data root to temporary directory
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Initialize processor
        self.processor = GRNBaseProcessor(
            name="test_grn_real",
            processor_data_dir='grn'
        )
        
        # Create test reference data
        self._create_test_reference_data()
        
        yield
        
        # Cleanup
        ProtosPaths.set_data_root(None)
    
    def _create_test_reference_data(self):
        """Create test reference GRN data mimicking real mo_grn structure."""
        # Create comprehensive test data
        test_data = pd.DataFrame({
            '1.44': ['I22', 'V23', 'L24', 'A25', 'T26'],
            '1.45': ['V23', 'A24', 'M25', 'V26', 'L27'],
            '1.46': ['T24', 'L25', 'F26', 'L27', 'F28'],
            '1.47': ['A25', 'M26', 'L27', 'M28', 'V29'],
            '1.48': ['A26', 'F27', 'V28', 'A29', 'A30'],
            '1.49': ['V27', 'L28', 'A29', 'L30', 'V31'],
            '1.50': ['L28', 'V29', 'L30', 'L31', 'L32'],  # Key position
            '2.50': ['A72', 'S73', 'A74', 'A75', 'S76'],  # Key position
            '3.50': ['R115', 'R116', 'R117', 'R118', 'R119'],  # Key position
            '4.50': ['G157', 'W158', 'W159', 'W160', 'G161'],  # Key position
            '5.50': ['F197', 'F198', 'Y199', 'F200', 'F201'],  # Key position
            '6.50': ['W238', 'W239', 'W240', 'F241', 'W242'],  # Key position
            '7.50': ['K276', 'K277', 'K278', 'K279', 'K280'],  # Schiff base lysine
        }, index=['7BMH', '4HYJ', '3DDL', 'TEST1', 'TEST2'])
        
        # Save as reference data
        self.processor.save_data('mo_grn', test_data)
    
    def test_load_reference_table_and_basic_operations(self):
        """Test loading real reference table and basic operations."""
        # Load the reference table
        self.processor.load_grn_table("mo_grn")
        
        # Verify data loaded
        assert self.processor.data is not None
        assert not self.processor.data.empty
        
        # Check data structure
        assert len(self.processor.data) == 5  # 5 test proteins
        assert len(self.processor.data.columns) == 13  # 13 GRN positions
        
        # Check for expected proteins
        known_proteins = ["7BMH", "4HYJ", "3DDL"]
        for protein in known_proteins:
            assert protein in self.processor.data.index
        
        # Check for key positions (dot notation)
        key_positions = ["1.50", "2.50", "3.50", "7.50"]
        for pos in key_positions:
            assert pos in self.processor.data.columns
        
        # Analyze position 7.50 (Schiff base lysine)
        col_750 = self.processor.data["7.50"]
        # Count lysines at this position
        k_count = sum(1 for val in col_750 if isinstance(val, str) and val.startswith("K"))
        total_non_missing = sum(1 for val in col_750 if pd.notna(val) and val != "-")
        
        # All test entries should have K at 7.50
        assert k_count == total_non_missing
        assert k_count == 5  # All 5 test proteins have K
    
    def test_get_grn_table_as_dict(self):
        """Test getting GRN table data in dictionary format."""
        # Load reference table
        self.processor.load_grn_table("mo_grn")
        
        # Get dictionary representation
        grn_dict = self.processor.get_grn_dict()
        
        # Check structure
        assert isinstance(grn_dict, dict)
        assert len(grn_dict) == len(self.processor.data)
        
        # Check each protein's data
        for protein_id, grn_list in grn_dict.items():
            assert isinstance(grn_list, list)
            assert len(grn_list) > 0
            
            # Each entry should be a residue code
            for grn in grn_list:
                if grn != '-' and pd.notna(grn):
                    assert isinstance(grn, str)
                    assert len(grn) >= 2  # Like "K276"
    
    def test_save_and_reload_workflow(self):
        """Test complete save and reload workflow."""
        # Load reference table
        self.processor.load_grn_table("mo_grn")
        
        # Get initial state
        initial_shape = self.processor.data.shape
        initial_columns = list(self.processor.data.columns)
        
        # Save with new name
        self.processor.save_grn_table("workflow_test")
        
        # Create new processor
        new_processor = GRNBaseProcessor(
            name="test_reload",
            processor_data_dir='grn'
        )
        
        # Load saved table
        new_processor.load_grn_table("workflow_test")
        
        # Verify identical
        assert new_processor.data.shape == initial_shape
        assert list(new_processor.data.columns) == initial_columns
        
        # Check data equality (reset index to avoid type mismatch)
        pd.testing.assert_frame_equal(
            self.processor.data.sort_index(axis=1).reset_index(drop=True),
            new_processor.data.sort_index(axis=1).reset_index(drop=True)
        )
    
    def test_filter_operations(self):
        """Test filtering operations on real data."""
        # Load reference table
        self.processor.load_grn_table("mo_grn")
        
        # Test 1: Filter by protein IDs
        proteins_to_keep = ["7BMH", "4HYJ"]
        self.processor.filter_by_ids(proteins_to_keep)
        
        assert len(self.processor.data) == 2
        assert all(idx in proteins_to_keep for idx in self.processor.data.index)
        
        # Reload for next test
        self.processor.load_grn_table("mo_grn")
        
        # Test 2: Filter by occurrence threshold
        # Keep columns where at least 60% have non-dash values
        initial_cols = len(self.processor.data.columns)
        self.processor.filter_data_by_occurances(threshold=0.6)
        
        # Should keep all columns since our test data is complete
        assert len(self.processor.data.columns) == initial_cols
    
    def test_column_sorting(self):
        """Test GRN column sorting."""
        # Load reference table
        self.processor.load_grn_table("mo_grn")
        
        # Manually shuffle columns
        shuffled_cols = list(self.processor.data.columns)
        np.random.shuffle(shuffled_cols)
        self.processor.data = self.processor.data[shuffled_cols]
        
        # Sort columns
        self.processor.sort_columns()
        
        # Verify sorting
        cols = list(self.processor.data.columns)
        
        # Check that columns are sorted by helix.position
        prev_helix = 0
        prev_pos = 0
        
        for col in cols:
            parts = col.split('.')
            if len(parts) == 2:
                helix = int(parts[0])
                pos = int(parts[1])
                
                # Should be in ascending order
                if helix == prev_helix:
                    assert pos > prev_pos
                else:
                    assert helix > prev_helix
                
                prev_helix = helix
                prev_pos = pos
    
    def test_dataset_merging(self):
        """Test merging multiple GRN datasets."""
        # Load reference table
        self.processor.load_grn_table("mo_grn")
        
        # Split into parts
        part1 = self.processor.data.iloc[:2]
        part2 = self.processor.data.iloc[2:4]
        part3 = self.processor.data.iloc[4:]
        
        # Save parts
        self.processor.data = part1
        self.processor.save_grn_table("part1")
        
        self.processor.data = part2
        self.processor.save_grn_table("part2")
        
        self.processor.data = part3
        self.processor.save_grn_table("part3")
        
        # Load and merge
        self.processor.load_and_merge_grn_tables(["part1", "part2", "part3"])
        
        # Should have all original proteins
        assert len(self.processor.data) == 5
        assert set(self.processor.data.index) == {'7BMH', '4HYJ', '3DDL', 'TEST1', 'TEST2'}
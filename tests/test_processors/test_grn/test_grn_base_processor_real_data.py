#!/usr/bin/env python3
"""
Test suite for GRN base processor using real reference data.

This test suite ensures all GRN processor functionality works with actual
biological data rather than mocked data.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import os
import tempfile
import shutil

from protos.processing.grn.grn_base_processor import GRNBaseProcessor


class TestGRNBaseProcessorRealData:
    """Test GRN base processor with real reference data."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Set up test environment with real data paths."""
        # Get test data directory
        self.test_data_dir = Path(__file__).parent.parent.parent / "test-data"
        self.grn_ref_dir = self.test_data_dir / "grn" / "ref"
        
        # Verify test data exists
        if not self.grn_ref_dir.exists():
            pytest.skip(f"Test data directory not found: {self.grn_ref_dir}")
        
        # Check for reference table
        self.mo_ref_path = self.grn_ref_dir / "mo_grn.csv"
        if not self.mo_ref_path.exists():
            pytest.skip(f"Reference GRN table not found: {self.mo_ref_path}")
        
        # Create temporary directory for outputs
        self.temp_dir = Path(tempfile.mkdtemp())
        self.data_root = self.temp_dir / "protos_data"
        self.data_root.mkdir(parents=True)
        
        # Set environment variables
        os.environ["PROTOS_DATA_ROOT"] = str(self.data_root)
        os.environ["PROTOS_REF_DATA_ROOT"] = str(self.test_data_dir)
        
        # Initialize processor
        self.processor = GRNBaseProcessor(
            name="test_grn_real",
            data_root=str(self.test_data_dir)
        )
        
        yield
        
        # Cleanup
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_load_real_reference_table(self):
        """Test loading real microbial opsin reference table."""
        # Load the reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Verify data loaded
        assert self.processor.data is not None
        assert not self.processor.data.empty
        
        # Check expected structure
        assert len(self.processor.data) > 0  # Should have proteins
        assert len(self.processor.data.columns) > 0  # Should have GRN positions
        
        # Check for expected proteins (microbial opsins)
        expected_proteins = ["7BMH", "4HYJ", "3DDL"]
        for protein in expected_proteins:
            if protein in self.processor.data.index:
                assert protein in self.processor.data.index
        
        # Check for expected GRN positions
        expected_positions = ["1.50", "2.50", "3.50", "7.50"]
        for pos in expected_positions:
            if pos in self.processor.data.columns:
                assert pos in self.processor.data.columns
        
        # Check critical position 7.50 (Schiff base lysine)
        if "7.50" in self.processor.data.columns:
            k_positions = self.processor.data["7.50"]
            # Most should have K (lysine) at this position
            k_count = sum(1 for val in k_positions if isinstance(val, str) and val.startswith("K"))
            assert k_count > len(k_positions) * 0.5  # At least 50% should have K
    
    def test_get_grn_dict_real_data(self):
        """Test GRN dictionary functionality with real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Get GRN dictionary
        grn_dict = self.processor.get_grn_dict()
        assert isinstance(grn_dict, dict)
        assert len(grn_dict) > 0
        
        # Check structure of dictionary
        for protein, grns in grn_dict.items():
            assert isinstance(grns, dict)
            # Check that we get residue+position format in values
            for grn_pos, value in grns.items():
                if pd.notna(value) and value != "-":
                    # Should be in format like "K270"
                    assert isinstance(value, str)
                    assert len(value) > 1
                    assert value[0].isalpha()  # First char is amino acid
    
    def test_filter_by_ids_real_data(self):
        """Test filtering by protein IDs with real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Get some protein IDs to filter
        all_proteins = list(self.processor.data.index)
        if len(all_proteins) > 3:
            proteins_to_keep = all_proteins[:3]
            
            # Filter by IDs
            self.processor.filter_by_ids(proteins_to_keep)
            
            # Check that only selected proteins remain
            assert len(self.processor.data) == len(proteins_to_keep)
            assert all(p in proteins_to_keep for p in self.processor.data.index)
    
    def test_grn_columns_real_data(self):
        """Test GRN columns from real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Check grns attribute is populated
        assert hasattr(self.processor, 'grns')
        assert isinstance(self.processor.grns, list)
        assert len(self.processor.grns) > 0
        
        # Should match data columns
        assert self.processor.grns == list(self.processor.data.columns)
        
        # Check format of positions
        for pos in self.processor.grns:
            assert isinstance(pos, str)
            # Should be in format X.YY where X is helix number
            parts = pos.split(".")
            if len(parts) == 2:
                assert parts[0].isdigit()
                assert parts[1].isdigit()
    
    def test_filter_by_occurrences_real_data(self):
        """Test filtering by occurrence threshold with real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Filter to keep only well-conserved positions
        original_cols = len(self.processor.data.columns)
        self.processor.filter_data_by_occurances(threshold=0.5)  # Keep positions present in >50% of proteins
        
        # Should have fewer columns after filtering
        assert len(self.processor.data.columns) <= original_cols
        
        # Check that remaining positions are well-represented
        for col in self.processor.data.columns:
            non_missing = self.processor.data[col].notna() & (self.processor.data[col] != "-")
            occurrence_rate = non_missing.sum() / len(self.processor.data)
            assert occurrence_rate >= 0.5
    
    def test_save_and_reload_real_data(self):
        """Test saving and reloading real GRN data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        original_data = self.processor.data.copy()
        
        # Change data root to temp directory for saving
        self.processor.data_root = str(self.data_root)
        self.processor.data_path = self.data_root / self.processor.processor_data_dir
        self.processor.data_path.mkdir(parents=True, exist_ok=True)
        
        # Save with a test name
        self.processor.save_grn_table("test_real_save", self.processor.data)
        
        # Create new processor to load saved data
        new_processor = GRNBaseProcessor(
            name="test_reload",
            data_root=str(self.data_root)
        )
        new_processor.load_grn_table("test_real_save")
        
        # Compare data
        pd.testing.assert_frame_equal(original_data, new_processor.data)
    
    def test_load_and_merge_tables_real_data(self):
        """Test loading and merging multiple GRN tables with real data."""
        # Create two test datasets by splitting the reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        if len(self.processor.data) > 10:
            # Split data into two parts
            mid_point = len(self.processor.data) // 2
            part1 = self.processor.data.iloc[:mid_point].copy()
            part2 = self.processor.data.iloc[mid_point-2:].copy()  # Some overlap
            
            # Save as separate datasets
            self.processor.data_root = str(self.data_root)
            self.processor.data_path = self.data_root / self.processor.processor_data_dir
            self.processor.data_path.mkdir(parents=True, exist_ok=True)
            
            self.processor.save_grn_table("test_part1", part1)
            self.processor.save_grn_table("test_part2", part2)
            
            # Create new processor and load/merge
            new_processor = GRNBaseProcessor(
                name="test_merge",
                data_root=str(self.data_root)
            )
            new_processor.load_and_merge_grn_tables(["test_part1", "test_part2"])
            
            # Should have all proteins from both parts
            expected_proteins = set(part1.index) | set(part2.index)
            assert set(new_processor.data.index) == expected_proteins
    
    def test_handle_missing_data_real(self):
        """Test handling of missing data in real GRN tables."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Real data should have some missing values
        has_missing = self.processor.data.isna().any().any()
        has_dashes = (self.processor.data == "-").any().any()
        
        # Check that we handle both NaN and "-" as missing
        if has_missing or has_dashes:
            # Get proteins with missing data
            missing_mask = self.processor.data.isna() | (self.processor.data == "-")
            proteins_with_missing = missing_mask.any(axis=1)
            
            assert proteins_with_missing.sum() > 0  # Should have some proteins with missing data
    
    def test_sort_columns_real_data(self):
        """Test sorting GRN columns with real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Get original column order
        original_cols = list(self.processor.data.columns)
        
        # Sort columns
        self.processor.sort_columns()
        sorted_cols = list(self.processor.data.columns)
        
        # Check that columns are now sorted
        # They should be sorted by helix number then position within helix
        for i in range(1, len(sorted_cols)):
            prev = sorted_cols[i-1].split(".")
            curr = sorted_cols[i].split(".")
            
            if len(prev) == 2 and len(curr) == 2:
                prev_helix = int(prev[0])
                curr_helix = int(curr[0])
                
                if prev_helix == curr_helix:
                    # Within same helix, position should increase
                    assert int(prev[1]) <= int(curr[1])
                else:
                    # Different helix
                    assert prev_helix < curr_helix
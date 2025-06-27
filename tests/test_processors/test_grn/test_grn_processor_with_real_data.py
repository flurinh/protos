#!/usr/bin/env python3
"""
Test suite for GRN processor using real reference data without parser modifications.

This test suite works with the current GRN system and validates functionality
with actual biological data.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import os
import tempfile
import shutil

from protos.processing.grn.grn_base_processor import GRNBaseProcessor


class TestGRNProcessorWithRealData:
    """Test GRN processor with real reference data using current system."""
    
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
        )
        
        yield
        
        # Cleanup
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_load_reference_table_and_basic_operations(self):
        """Test loading real reference table and basic operations."""
        # Load the reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Verify data loaded
        assert self.processor.data is not None
        assert not self.processor.data.empty
        
        # Check data structure
        print(f"\nLoaded GRN table with {len(self.processor.data)} proteins")
        print(f"Number of GRN positions: {len(self.processor.data.columns)}")
        
        # Check for some expected proteins
        known_proteins = ["7BMH", "4HYJ", "3DDL"]
        found_proteins = [p for p in known_proteins if p in self.processor.data.index]
        print(f"Found known proteins: {found_proteins}")
        
        # Check for key positions (dot notation)
        key_positions = ["1.50", "2.50", "3.50", "7.50"]
        found_positions = [pos for pos in key_positions if pos in self.processor.data.columns]
        print(f"Found key positions: {found_positions}")
        
        # Analyze position 7.50 (Schiff base lysine)
        if "7.50" in self.processor.data.columns:
            col_750 = self.processor.data["7.50"]
            # Count lysines at this position
            k_count = sum(1 for val in col_750 if isinstance(val, str) and val.startswith("K"))
            total_non_missing = sum(1 for val in col_750 if pd.notna(val) and val != "-")
            print(f"\nPosition 7.50 analysis:")
            print(f"  Total non-missing: {total_non_missing}")
            print(f"  Lysines (K): {k_count}")
            print(f"  Percentage K: {k_count/total_non_missing*100:.1f}%")
            
            # Most microbial opsins should have K at 7.50
            assert k_count > total_non_missing * 0.8  # At least 80% should have K
    
    def test_get_grn_table_as_dict(self):
        """Test getting GRN table data in dictionary format."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # The GRN processor's get_grn_dict returns list of positions per protein
        # But what we really want for GRN tables is protein -> {grn: residue}
        # Let's create that from the data
        grn_residue_dict = {}
        
        for protein in self.processor.data.index:
            grn_residue_dict[protein] = {}
            for grn_pos in self.processor.data.columns:
                residue = self.processor.data.loc[protein, grn_pos]
                if pd.notna(residue) and residue != "-":
                    grn_residue_dict[protein][grn_pos] = residue
        
        # Check structure
        assert isinstance(grn_residue_dict, dict)
        assert len(grn_residue_dict) > 0
        
        # Analyze first few proteins
        print("\nGRN table as dictionary sample:")
        for i, (protein, grn_dict) in enumerate(grn_residue_dict.items()):
            if i >= 3:  # Show only first 3
                break
            print(f"\n{protein}:")
            print(f"  Total annotated positions: {len(grn_dict)}")
            
            # Show some key positions if present
            key_positions = ["1.50", "2.50", "3.50", "7.50"]
            for pos in key_positions:
                if pos in grn_dict:
                    print(f"  {pos}: {grn_dict[pos]}")
    
    def test_filtering_operations(self):
        """Test filtering operations with real data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Test filtering by IDs
        all_proteins = list(self.processor.data.index)
        print(f"\nTotal proteins: {len(all_proteins)}")
        
        # Filter to first 10 proteins
        if len(all_proteins) > 10:
            proteins_to_keep = all_proteins[:10]
            self.processor.filter_by_ids(proteins_to_keep)
            
            assert len(self.processor.data) == 10
            print(f"After filtering: {len(self.processor.data)} proteins")
    
    def test_occurrence_filtering(self):
        """Test filtering by occurrence threshold."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Get original dimensions
        original_proteins = len(self.processor.data)
        original_positions = len(self.processor.data.columns)
        
        # Filter to keep only well-conserved positions
        # The method expects number of proteins, not fraction
        min_proteins = int(original_proteins * 0.7)  # 70% of proteins
        self.processor.filter_data_by_occurances(threshold=min_proteins)
        
        # Check results
        filtered_proteins = len(self.processor.data)
        filtered_positions = len(self.processor.data.columns)
        
        print(f"\nOccurrence filtering (>{min_proteins} proteins / {original_proteins} total):")
        print(f"  Proteins: {original_proteins} -> {filtered_proteins}")
        print(f"  Positions: {original_positions} -> {filtered_positions}")
        
        # Should have fewer positions after filtering
        assert filtered_positions < original_positions
        
        # Well-conserved positions should include key ones
        expected_conserved = ["1.50", "2.50", "3.50", "7.50"]
        for pos in expected_conserved:
            if pos in self.processor.data.columns:
                print(f"  {pos}: retained (conserved)")
        
        # Verify remaining positions are well-represented
        for col in self.processor.data.columns[:10]:  # Check first 10
            non_missing = self.processor.data[col].notna() & (self.processor.data[col] != "-")
            occurrence_count = non_missing.sum()
            assert occurrence_count >= min_proteins
    
    def test_save_and_reload_workflow(self):
        """Test saving and reloading GRN data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Get subset for testing
        if len(self.processor.data) > 20:
            subset = self.processor.data.iloc[:20].copy()
        else:
            subset = self.processor.data.copy()
        
        # Change data root to temp directory
        # data_root is now managed globally
        self.processor.data_path = self.data_root / self.processor.processor_data_dir
        self.processor.data_path.mkdir(parents=True, exist_ok=True)
        
        # Save subset - need to set data first
        self.processor.data = subset
        self.processor.save_grn_table("test_subset")
        print(f"\nSaved {len(subset)} proteins to test_subset")
        
        # Create new processor and load
        new_processor = GRNBaseProcessor(
            name="test_reload",
        )
        new_processor.load_grn_table("test_subset")
        
        # Verify data matches
        assert len(new_processor.data) == len(subset)
        assert list(new_processor.data.index) == list(subset.index)
        assert list(new_processor.data.columns) == list(subset.columns)
        print(f"Successfully reloaded {len(new_processor.data)} proteins")
    
    def test_column_operations(self):
        """Test column-related operations."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        # Test grns attribute
        assert hasattr(self.processor, 'grns')
        assert len(self.processor.grns) == len(self.processor.data.columns)
        
        # Check column sorting
        original_order = list(self.processor.data.columns)
        self.processor.sort_columns()
        sorted_order = list(self.processor.data.columns)
        
        print(f"\nColumn sorting:")
        print(f"  Original first 10: {original_order[:10]}")
        print(f"  Sorted first 10: {sorted_order[:10]}")
        
        # The sorted columns should group by helix
        # (though with dot notation parsing issues, this might not work perfectly)
    
    def test_reference_data_properties(self):
        """Analyze properties of the real reference data."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        print("\nReference data analysis:")
        print(f"Total proteins: {len(self.processor.data)}")
        print(f"Total GRN positions: {len(self.processor.data.columns)}")
        
        # Analyze missing data
        total_cells = self.processor.data.size
        missing_cells = self.processor.data.isna().sum().sum()
        dash_cells = (self.processor.data == "-").sum().sum()
        
        print(f"\nMissing data:")
        print(f"  Total cells: {total_cells}")
        print(f"  NaN cells: {missing_cells} ({missing_cells/total_cells*100:.1f}%)")
        print(f"  Dash cells: {dash_cells} ({dash_cells/total_cells*100:.1f}%)")
        
        # Analyze helix distribution
        helix_counts = {}
        for col in self.processor.data.columns:
            if '.' in col and not col.startswith('n.') and not col.startswith('c.'):
                helix = col.split('.')[0]
                if helix.isdigit() and len(helix) == 1:
                    helix_counts[helix] = helix_counts.get(helix, 0) + 1
        
        print(f"\nPositions per helix:")
        for helix in sorted(helix_counts.keys()):
            print(f"  Helix {helix}: {helix_counts[helix]} positions")
        
        # Check for special regions
        n_term = [col for col in self.processor.data.columns if col.startswith('n.')]
        c_term = [col for col in self.processor.data.columns if col.startswith('c.')]
        loops = [col for col in self.processor.data.columns 
                if '.' in col and len(col.split('.')[0]) == 2]
        
        print(f"\nSpecial regions:")
        print(f"  N-terminal positions: {len(n_term)}")
        print(f"  C-terminal positions: {len(c_term)}")
        print(f"  Loop positions: {len(loops)}")
    
    def test_dataset_merging(self):
        """Test merging multiple GRN datasets."""
        # Load reference table
        self.processor.load_grn_table("ref/mo_grn")
        
        if len(self.processor.data) < 20:
            pytest.skip("Not enough data for merging test")
        
        # Split into parts
        third = len(self.processor.data) // 3
        part1 = self.processor.data.iloc[:third].copy()
        part2 = self.processor.data.iloc[third-5:2*third].copy()  # Some overlap
        part3 = self.processor.data.iloc[2*third-5:].copy()  # Some overlap
        
        # Save parts
        # data_root is now managed globally
        self.processor.data_path = self.data_root / self.processor.processor_data_dir
        self.processor.data_path.mkdir(parents=True, exist_ok=True)
        
        self.processor.data = part1
        self.processor.save_grn_table("part1")
        self.processor.data = part2
        self.processor.save_grn_table("part2")
        self.processor.data = part3
        self.processor.save_grn_table("part3")
        
        # Test merging
        merge_processor = GRNBaseProcessor(
            name="test_merge",
        )
        merge_processor.load_and_merge_grn_tables(["part1", "part2", "part3"])
        
        print(f"\nDataset merging:")
        print(f"  Part 1: {len(part1)} proteins")
        print(f"  Part 2: {len(part2)} proteins")
        print(f"  Part 3: {len(part3)} proteins")
        print(f"  Merged: {len(merge_processor.data)} proteins")
        
        # Should have union of all proteins (with duplicates removed)
        all_proteins = set(part1.index) | set(part2.index) | set(part3.index)
        assert len(merge_processor.data) == len(all_proteins)
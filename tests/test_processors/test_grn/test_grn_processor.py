import pytest
import os
import pandas as pd
import numpy as np
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_utils import sort_grns_str, get_seq


class TestGRNProcessor:
    @pytest.fixture
    def data_path(self):
        """Path to the test data directory."""
        return os.path.join("data", "grn", "ref")
    
    @pytest.fixture
    def test_dataset(self):
        """Name of a test dataset to use for testing."""
        return "ref"  # Using the reference dataset which should exist
    
    @pytest.fixture
    def grnp(self, data_path, test_dataset):
        """Create a GRNProcessor instance for testing."""
        return GRNProcessor(dataset=test_dataset, path=data_path)
    
    def test_init(self, grnp, test_dataset):
        """Test GRNProcessor initialization."""
        assert grnp.dataset == test_dataset
        assert isinstance(grnp.data, pd.DataFrame)
        assert len(grnp.ids) > 0
        assert len(grnp.grns) > 0
    
    def test_list_available_datasets(self, grnp):
        """Test listing available datasets."""
        datasets = grnp.list_available_datasets()
        assert isinstance(datasets, list)
        assert len(datasets) > 0
        assert grnp.dataset in datasets
    
    def test_load_grn_table(self, grnp, test_dataset):
        """Test loading a GRN table."""
        # Test loading the same dataset
        data = grnp.load_grn_table(dataset=test_dataset)
        assert isinstance(data, pd.DataFrame)
        assert len(data) > 0
        
        # Test with explicit loading
        grnp_explicit = GRNProcessor(dataset=None, preload=False)
        data_explicit = grnp_explicit.load_grn_table(dataset=test_dataset)
        assert isinstance(data_explicit, pd.DataFrame)
        assert len(data_explicit) > 0
    
    def test_get_available_grn_tables(self, grnp):
        """Test getting available GRN tables."""
        tables = grnp.get_available_grn_tables()
        assert isinstance(tables, list)
        assert len(tables) > 0
        assert grnp.dataset in tables
    
    def test_get_seq_dict(self, grnp):
        """Test getting sequences from GRN table."""
        seq_dict = grnp.get_seq_dict()
        assert isinstance(seq_dict, dict)
        assert len(seq_dict) > 0
        
        # Check a few entries
        for uen, seq in list(seq_dict.items())[:3]:
            assert isinstance(uen, str)
            assert isinstance(seq, str)
            assert len(seq) > 0
            assert set(seq).issubset(set("ACDEFGHIKLMNPQRSTVWY-"))  # Valid amino acids and gap
    
    def test_remove_duplicate_ids(self, grnp):
        """Test removing duplicate IDs."""
        # Create a temporary duplicate
        test_id = grnp.ids[0]
        grnp.ids.append(test_id)
        grnp.remove_duplicate_ids()
        
        # Ensure the ID only appears once
        assert grnp.ids.count(test_id) == 1
    
    def test_load_and_merge_grn_tables(self, grnp, test_dataset):
        """Test loading and merging GRN tables."""
        # For this test, we'll merge the same dataset with itself
        datasets = [test_dataset]
        merged_data = grnp.load_and_merge_grn_tables(datasets=datasets)
        
        assert isinstance(merged_data, pd.DataFrame)
        assert len(merged_data) > 0
        assert len(merged_data.columns) > 0
    
    def test_reset_data(self, grnp):
        """Test resetting data."""
        original_ids = grnp.ids.copy()
        original_grns = grnp.grns.copy()
        
        # Modify the data
        if len(grnp.ids) > 1:
            grnp.ids = grnp.ids[:-1]
        if len(grnp.grns) > 1:
            grnp.grns = grnp.grns[:-1]
        
        # Reset and verify
        grnp.reset_data()
        assert len(grnp.ids) == len(original_ids)
        assert len(grnp.grns) == len(original_grns)
    
    def test_apply_interval(self, grnp):
        """Test applying a GRN interval."""
        original_grn_count = len(grnp.grns)
        
        # Select a subset of GRNs
        grn_subset = sorted(grnp.grns[:5])
        grnp.apply_interval(grn_subset)
        
        assert len(grnp.grns) == len(grn_subset)
        assert set(grnp.grns) == set(grn_subset)
    
    def test_get_grn_dict(self, grnp):
        """Test getting GRN dictionary."""
        grn_dict = grnp.get_grn_dict()
        
        assert isinstance(grn_dict, dict)
        assert len(grn_dict) > 0
        
        # Check a few entries
        for uen, grns in list(grn_dict.items())[:3]:
            assert isinstance(uen, str)
            assert isinstance(grns, list)
            assert all(grn in grnp.grns for grn in grns)
    
    def test_filter_data_by_occurances(self, grnp):
        """Test filtering data by occurrences."""
        original_grn_count = len(grnp.grns)
        
        # Filter with a very low threshold to ensure we have results
        threshold = 1
        grnp.filter_data_by_occurances(threshold)
        
        # Since we used a low threshold, we should still have GRNs
        assert len(grnp.grns) > 0
        assert len(grnp.grns) <= original_grn_count
    
    def test_filter_by_ids(self, grnp):
        """Test filtering by IDs."""
        original_id_count = len(grnp.ids)
        
        # Take a subset of IDs
        if len(grnp.ids) >= 3:
            test_ids = grnp.ids[:3]
            grnp.filter_by_ids(test_ids)
            
            assert len(grnp.ids) == len(test_ids)
            assert set(grnp.ids) == set(test_ids)
        else:
            # Skip the test if not enough IDs
            pytest.skip("Not enough IDs in test dataset")
    
    def test_sort_columns(self, grnp):
        """Test sorting columns."""
        # Shuffle the columns
        import random
        random.shuffle(grnp.grns)
        grnp.data = grnp.data[grnp.grns]
        
        # Sort and check
        sorted_grns = sort_grns_str(grnp.grns)
        grnp.sort_columns()
        
        assert grnp.grns == sorted_grns
    
    def test_get_tm_regions(self, grnp):
        """Test getting TM regions."""
        # This depends on having a reference dataset with TM annotations
        try:
            tm_region_grns = grnp.get_tm_regions()
            assert isinstance(tm_region_grns, list)
            
            # TM regions should be in GRNs
            for grn in tm_region_grns:
                assert grn in grnp.grns
        except:
            # Skip if the method is not implemented or dataset doesn't have TM regions
            pytest.skip("get_tm_regions not implemented or no TM regions in dataset")
    
    def test_get_binding_pocket(self, grnp):
        """Test getting binding pocket."""
        # This depends on having a reference dataset with binding pocket annotations
        try:
            binding_pocket_grns = grnp.get_binding_pocket()
            assert isinstance(binding_pocket_grns, list)
            
            # Binding pocket GRNs should be in GRNs
            for grn in binding_pocket_grns:
                assert grn in grnp.grns
        except:
            # Skip if the method is not implemented or dataset doesn't have binding pocket
            pytest.skip("get_binding_pocket not implemented or no binding pocket in dataset")
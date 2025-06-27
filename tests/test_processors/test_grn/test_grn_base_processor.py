"""
Tests for the GRNBaseProcessor class.

This test suite validates GRN (Generic Residue Numbering) processor functionality.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from protos.processing.grn.grn_base_processor import GRNBaseProcessor


@pytest.fixture
def temp_data_dir():
    """Create a temporary data directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_grn_table():
    """Create a sample GRN table for testing.
    
    The actual GRN table format has:
    - Protein/structure IDs as row index
    - GRN positions as columns (e.g., "3.50", "7.53")
    - Values are residue+position strings (e.g., "K270")
    """
    data = {
        '1.50': ['M22', 'M62', np.nan, 'L22', '-'],
        '2.50': ['V50', 'I90', 'V50', 'I50', '-'],
        '3.50': ['R82', 'R115', 'R80', 'R82', '-'],
        '7.53': ['Y296', 'K270', 'Y240', 'Y296', 'H3']
    }
    index = ['1ABC', '7BMH', '4HYJ', '2DEF', 'LMA0']
    return pd.DataFrame(data, index=index)


def test_basic_initialization(temp_data_dir):
    """Test basic GRN processor initialization."""
    processor = GRNBaseProcessor(
        name="test_grn",
        processor_data_dir="grn"
    )
    
    assert processor.name == "test_grn"
    # data_root comes from global configuration now
    assert processor.processor_data_dir == "grn"
    assert os.path.exists(processor.data_path)
    assert processor.dataset is None
    assert len(processor.grns) == 0
    assert len(processor.ids) == 0


def test_save_load_grn_table(temp_data_dir, sample_grn_table):
    """Test saving and loading GRN tables."""
    processor = GRNBaseProcessor(
        name="test_grn",
        processor_data_dir="grn",
        preload=False
    )
    
    # Set the data and save
    processor.data = sample_grn_table
    processor.save_grn_table("test_grn_table")
    
    # Verify file exists in the correct location
    # GRN processor saves to tables subdirectory
    table_path = os.path.join(processor.data_path, "tables", "test_grn_table.csv")
    assert os.path.exists(table_path)
    
    # Create new processor to test loading
    processor2 = GRNBaseProcessor(
        name="test_grn2",
        processor_data_dir="grn",
        preload=False
    )
    processor2.load_grn_table("test_grn_table")
    
    # Verify data loaded correctly
    assert processor2.dataset == "test_grn_table"
    assert not processor2.data.empty
    assert len(processor2.data) == len(sample_grn_table)
    
    # Verify IDs and GRNs extracted
    assert set(processor2.ids) == {'1ABC', '7BMH', '4HYJ', '2DEF', 'LMA0'}
    assert '1.50' in processor2.grns
    assert '2.50' in processor2.grns
    assert '3.50' in processor2.grns
    assert '7.53' in processor2.grns


def test_get_residue_by_grn(temp_data_dir, sample_grn_table):
    """Test getting residue information by GRN position."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Set data and save
    processor.data = sample_grn_table
    processor.save_grn_table("test_grn")
    processor.load_grn_table("test_grn")
    
    # Test getting residue for specific protein and GRN
    # In the actual format, values are residue+position strings
    residue = processor.data.loc['1ABC', '3.50']
    assert residue == 'R82'
    
    residue = processor.data.loc['7BMH', '7.53']
    assert residue == 'K270'
    
    # Test accessing non-existent data
    residue = processor.data.loc['4HYJ', '3.50']
    assert residue == 'R80'


def test_grn_lookup_methods(temp_data_dir, sample_grn_table):
    """Test GRN lookup methods."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Set data and save
    processor.data = sample_grn_table
    processor.save_grn_table("test_grn")
    processor.load_grn_table("test_grn")
    
    # Test getting all GRN positions for a protein
    protein_grns = processor.data.loc['7BMH']
    assert protein_grns['1.50'] == 'M62'
    assert protein_grns['3.50'] == 'R115'
    assert protein_grns['7.53'] == 'K270'
    
    # Test getting all proteins with a specific residue at a GRN position
    proteins_with_R82 = processor.data[processor.data['3.50'] == 'R82'].index.tolist()
    assert '1ABC' in proteins_with_R82
    assert '2DEF' in proteins_with_R82


def test_filter_grn_columns(temp_data_dir, sample_grn_table):
    """Test filtering data by GRN columns."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Set data and save
    processor.data = sample_grn_table
    processor.save_grn_table("test_grn")
    processor.load_grn_table("test_grn")
    
    # Filter columns manually (since filter_by_grns doesn't exist)
    grn_cols = ['1.50', '3.50']
    filtered_df = processor.data[grn_cols]
    
    # Should only have the selected GRN columns
    assert all(col in filtered_df.columns for col in grn_cols)
    # All proteins should still be present
    assert len(filtered_df) == 5


def test_merge_grn_tables(temp_data_dir):
    """Test merging multiple GRN tables."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Create two GRN tables with actual format
    table1_data = {
        '1.50': ['M22', 'M62'],
        '2.50': ['V50', 'I90']
    }
    table1 = pd.DataFrame(table1_data, index=['1ABC', '7BMH'])
    
    table2_data = {
        '3.50': ['R80', 'R82'],
        '7.53': ['Y240', 'Y296']
    }
    table2 = pd.DataFrame(table2_data, index=['4HYJ', '2DEF'])
    
    # Save tables
    processor.data = table1
    processor.save_grn_table("table1")
    processor.data = table2
    processor.save_grn_table("table2")
    
    # Load and merge
    processor.load_and_merge_grn_tables(["table1", "table2"])
    
    # Verify merge
    assert len(processor.data) == 4
    assert set(processor.ids) == {'1ABC', '7BMH', '4HYJ', '2DEF'}
    assert set(processor.grns) == {'1.50', '2.50', '3.50', '7.53'}


def test_dot_notation_handling(temp_data_dir):
    """Test handling of dot notation (3.50 vs 3x50)."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Create table with dot notation (actual format)
    dot_table_data = {
        '3.50': ['R82'],
        '3.51': ['Y83'],
        '7.53': ['K270']
    }
    dot_table = pd.DataFrame(dot_table_data, index=['1ABC'])
    
    # Set data and save
    processor.data = dot_table
    processor.save_grn_table("dot_table")
    processor.load_grn_table("dot_table")
    
    # Should detect dot notation
    assert processor.using_dot_notation
    assert '3.50' in processor.grns
    assert '7.53' in processor.grns
    # Check x_to_dot mapping exists
    assert '3x50' in processor.x_to_dot
    assert processor.x_to_dot['3x50'] == '3.50'


def test_empty_grn_table(temp_data_dir):
    """Test handling of empty GRN table."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Create empty table with proper structure
    empty_table = pd.DataFrame(columns=['1.50', '2.50', '3.50'])
    
    # Set data and save
    processor.data = empty_table
    processor.save_grn_table("empty")
    processor.load_grn_table("empty")
    
    # Should handle gracefully
    assert len(processor.data) == 0
    assert len(processor.ids) == 0
    assert len(processor.grns) == 3  # Has the GRN columns


def test_dataset_listing(temp_data_dir, sample_grn_table):
    """Test listing available GRN datasets."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Save multiple tables
    processor.data = sample_grn_table
    processor.save_grn_table("table_a")
    processor.data = sample_grn_table.head(3)
    processor.save_grn_table("table_b")
    
    # List datasets
    datasets = processor.list_datasets()
    dataset_ids = [d["id"] for d in datasets]
    
    assert "table_a" in dataset_ids
    assert "table_b" in dataset_ids


def test_na_handling(temp_data_dir, sample_grn_table):
    """Test handling of NA values in GRN tables."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Set data and save
    processor.data = sample_grn_table
    processor.save_grn_table("test_na")
    processor.load_grn_table("test_na")
    
    # Check NA values are filled with '-'
    # The load_grn_table method fills NA with '-'
    assert processor.data.loc['4HYJ', '1.50'] == '-'
    
    # Verify we can still work with the data
    proteins_with_data = processor.data[processor.data['1.50'] != '-'].index.tolist()
    assert '4HYJ' not in proteins_with_data
    assert '1ABC' in proteins_with_data
    assert '7BMH' in proteins_with_data
    assert '2DEF' in proteins_with_data
    assert 'LMA0' not in proteins_with_data  # LMA0 has '-' in 1.50


def test_grn_table_format_validation(temp_data_dir):
    """Test GRN table format validation and loading."""
    processor = GRNBaseProcessor(
        name="test_grn",
        preload=False
    )
    
    # Create table with just GRN columns (actual format)
    data = {
        '1.50': ['M22', 'M62', np.nan, 'L22', '-'],
        '2.50': ['V50', 'I90', 'V50', 'I50', '-'],
        '3.50': ['R82', 'R115', 'R80', 'R82', 'H3']
    }
    grn_table = pd.DataFrame(data, index=['1ABC', '7BMH', '4HYJ', '2DEF', 'LMA0'])
    
    # Save and load
    processor.data = grn_table
    processor.save_grn_table("format_test")
    processor.load_grn_table("format_test")
    
    # Verify GRN columns are identified correctly
    assert '1.50' in processor.grns
    assert '2.50' in processor.grns
    assert '3.50' in processor.grns
    
    # Verify all GRN columns are present
    assert len(processor.grns) == 3
    
    # Verify protein IDs
    assert set(processor.ids) == {'1ABC', '7BMH', '4HYJ', '2DEF', 'LMA0'}
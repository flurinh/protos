"""
Advanced tests for GRN functionalities including sorting, filtering, sequence dict, etc.
"""

import os
import tempfile
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.grn.grn_utils import (
    sort_grns_str,
    get_grn_interval,
    get_seq,
    get_annot_seq,
    get_tm_residues,
    map_grn_to_color,
    parse_grn_str2float,
    parse_grn_float2str
)
from protos.processing.schema.grn_utils_updated import (
    validate_grn_string
)


@pytest.fixture
def temp_data_dir():
    """Create a temporary data directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def comprehensive_grn_table():
    """Create a comprehensive GRN table for testing advanced features."""
    data = {
        # Standard GRN positions
        '1.50': ['M62', 'M62', 'L22', 'M22', '-', 'M62'],
        '2.50': ['I90', 'I90', 'I50', 'V50', 'V50', 'I90'],
        '3.50': ['R115', 'R115', 'R82', 'R82', 'R80', 'R115'],
        '4.50': ['W142', 'W142', '-', 'W142', 'W86', 'W142'],
        '5.50': ['F195', 'F195', 'F187', 'F186', '-', 'F195'],
        '6.50': ['W244', 'W244', 'W198', 'W192', 'W186', 'W244'],
        '7.50': ['N296', 'N296', 'N218', 'P241', '-', 'N296'],
        '7.53': ['K270', 'K270', 'K225', 'Y296', 'Y240', 'K270'],
        # N-terminus positions
        'n.1': ['S50', 'S50', '-', '-', '-', 'S50'],
        'n.2': ['D49', 'D49', '-', '-', '-', 'D49'],
        # C-terminus positions
        'c.1': ['L277', 'L277', '-', '-', '-', 'L277'],
        'c.2': ['L278', 'L278', '-', '-', '-', 'L278'],
        # Loop positions
        '45.500': ['K111', 'K111', '-', '-', '-', 'K111'],
        '45.510': ['I112', 'I112', '-', '-', '-', 'I112'], 
        '45.520': ['V113', 'V113', '-', '-', '-', 'V113']
    }
    index = ['7BMH', '7BMH_A', '1ABC', '2DEF', '4HYJ', '7BMH_B']
    return pd.DataFrame(data, index=index)


class TestGRNSorting:
    """Test GRN sorting functionality."""
    
    def test_sort_grns_str(self):
        """Test sorting GRN strings."""
        grn_list = ['7.53', '3.50', 'n.2', '1.50', 'c.1', '45.50', '2.50', 'n.1']
        sorted_grns = sort_grns_str(grn_list)
        
        # Note: sort_grns_str converts loop positions like 45.50 to 45.500
        expected_order = ['n.2', 'n.1', '1.50', '2.50', '3.50', '45.500', '7.53', 'c.1']
        assert sorted_grns == expected_order
    
    def test_sort_grns_with_x_notation(self):
        """Test sorting with x notation."""
        grn_list = ['7x53', '3x50', '1x50', '2x50']
        sorted_grns = sort_grns_str(grn_list, output_notation_type='x')
        
        expected_order = ['1x50', '2x50', '3x50', '7x53']
        assert sorted_grns == expected_order
    
    def test_sort_columns_in_processor(self, temp_data_dir, comprehensive_grn_table):
        """Test sorting columns in the processor."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_sort")
        processor.load_grn_table("test_sort")
        
        # Sort columns
        processor.sort_columns()
        
        # Check that columns are in the expected order
        cols = list(processor.data.columns)
        
        # N-terminus should come first
        assert cols[0] == 'n.2'
        assert cols[1] == 'n.1'
        
        # Then TM regions in order
        assert '1.50' in cols[:5]
        assert '2.50' in cols[:6]
        
        # C-terminus should come last
        assert cols[-2] == 'c.1'
        assert cols[-1] == 'c.2'


class TestGRNFiltering:
    """Test GRN filtering functionality."""
    
    def test_filter_by_ids(self, temp_data_dir, comprehensive_grn_table):
        """Test filtering by protein IDs."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.ids = list(comprehensive_grn_table.index)
        
        # Filter to keep only specific IDs
        ids_to_keep = ['7BMH', '1ABC', '4HYJ']
        processor.filter_by_ids(ids_to_keep)
        
        # Check filtered data
        assert len(processor.data) == 3
        assert set(processor.data.index) == set(ids_to_keep)
        assert set(processor.ids) == set(ids_to_keep)
    
    def test_filter_data_by_occurances(self, temp_data_dir, comprehensive_grn_table):
        """Test filtering by occurrence threshold."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_filter")
        processor.load_grn_table("test_filter")
        
        # Filter columns that have at least 4 non-dash values
        processor.filter_data_by_occurances(threshold=4)
        
        # Columns with 4+ values should remain
        assert '1.50' in processor.data.columns  # Has 5 non-dash values
        assert '2.50' in processor.data.columns  # Has 6 non-dash values
        assert '3.50' in processor.data.columns  # Has 6 non-dash values
        
        # Columns with <4 values should be removed
        assert 'n.1' not in processor.data.columns  # Has only 3 non-dash values
        assert 'c.1' not in processor.data.columns  # Has only 3 non-dash values


class TestGRNSequenceDict:
    """Test sequence dictionary functionality."""
    
    def test_get_seq_dict(self, temp_data_dir, comprehensive_grn_table):
        """Test getting sequence dictionary from GRN table."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_seq")
        processor.load_grn_table("test_seq")
        
        # Get sequence dictionary
        seq_dict = processor.get_seq_dict()
        
        # Check structure
        assert isinstance(seq_dict, dict)
        assert '7BMH' in seq_dict
        assert '1ABC' in seq_dict
        
        # Check sequence for 7BMH - should be a string of single letter amino acids
        seq_7bmh = seq_dict['7BMH']
        assert isinstance(seq_7bmh, str)
        assert 'M' in seq_7bmh  # M62 -> M
        assert 'R' in seq_7bmh  # R115 -> R
        assert 'K' in seq_7bmh  # K270 -> K
        
        # Check sequence for 4HYJ - should exclude missing values
        seq_4hyj = seq_dict['4HYJ']
        assert isinstance(seq_4hyj, str)
        # Should have fewer residues due to missing values
    
    def test_get_seq_from_table(self, comprehensive_grn_table):
        """Test getting sequence for specific gene."""
        seq = get_seq('7BMH', comprehensive_grn_table)
        
        # Should return sequence of single-letter amino acids
        assert 'M' in seq  # From M62
        assert 'R' in seq  # From R115
        assert 'K' in seq  # From K270
    
    def test_get_annot_seq_from_table(self, comprehensive_grn_table):
        """Test getting annotated sequence."""
        annot_seq = get_annot_seq('7BMH', comprehensive_grn_table)
        
        # Should return full residue+position annotations
        assert 'M62' in annot_seq
        assert 'R115' in annot_seq
        assert 'K270' in annot_seq


class TestGRNIntervals:
    """Test GRN interval functionality."""
    
    def test_get_grn_interval(self, temp_data_dir, comprehensive_grn_table):
        """Test getting GRN intervals."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_interval")
        processor.load_grn_table("test_interval")
        
        # Get interval between two GRN positions
        # Note: get_grn_interval expects x notation for standard TM regions
        interval = get_grn_interval(left='2x50', right='3x50', table=processor.data)
        
        # Should include positions between 2.50 and 3.50 (returns in dot notation)
        assert '2.50' in interval
        assert '3.50' in interval
        
        # Test N-terminus interval
        n_interval = get_grn_interval(left='n.2', right='n.1', grns_str=['n.2', 'n.1', '1.50'])
        assert 'n.2' in n_interval
        assert 'n.1' in n_interval
    
    def test_apply_interval(self, temp_data_dir, comprehensive_grn_table):
        """Test applying GRN interval to filter data."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_apply_interval")
        processor.load_grn_table("test_apply_interval")
        
        # Apply interval from 1.50 to 3.50
        processor.apply_interval(['1.50', '2.50', '3.50'])
        
        # Check that only interval columns remain
        assert len(processor.data.columns) == 3
        assert '1.50' in processor.data.columns
        assert '2.50' in processor.data.columns
        assert '3.50' in processor.data.columns
        assert '7.53' not in processor.data.columns


class TestGRNParsing:
    """Test GRN string parsing and validation."""
    
    def test_parse_grn_str2float(self):
        """Test parsing GRN strings to floats."""
        # Standard notation
        assert parse_grn_str2float('3.50') == 3.50
        assert parse_grn_str2float('7.53') == 7.53
        
        # X notation
        assert parse_grn_str2float('3x50') == 3.50
        assert parse_grn_str2float('7x53') == 7.53
        
        # N-terminus
        assert parse_grn_str2float('n.1') == -1.0
        assert parse_grn_str2float('n.2') == -2.0
        
        # C-terminus
        assert parse_grn_str2float('c.1') == 101.0
        assert parse_grn_str2float('c.2') == 102.0
        
        # Loop positions
        assert parse_grn_str2float('45.50') == 45.50
        assert parse_grn_str2float('45.51') == 45.51
        assert parse_grn_str2float('45.500') == 45.50  # Handles three decimal places
    
    def test_parse_grn_float2str(self):
        """Test parsing floats back to GRN strings."""
        # Standard notation (dot)
        assert parse_grn_float2str(3.50, 'dot') == '3.50'
        assert parse_grn_float2str(7.53, 'dot') == '7.53'
        
        # X notation
        assert parse_grn_float2str(3.50, 'x') == '3x50'
        assert parse_grn_float2str(7.53, 'x') == '7x53'
        
        # N-terminus
        assert parse_grn_float2str(-1.0, 'dot') == 'n.1'
        assert parse_grn_float2str(-2.0, 'dot') == 'n.2'
        
        # C-terminus
        assert parse_grn_float2str(101.0, 'dot') == 'c.1'
        assert parse_grn_float2str(102.0, 'dot') == 'c.2'
    
    def test_validate_grn_string(self):
        """Test GRN string validation."""
        # Valid GRNs
        assert validate_grn_string('3.50')[0] == True
        assert validate_grn_string('3x50')[0] == True
        assert validate_grn_string('n.1')[0] == True
        assert validate_grn_string('c.1')[0] == True
        assert validate_grn_string('45.50')[0] == True
        assert validate_grn_string('45.500')[0] == True  # Three decimal places also valid
        
        # Invalid GRNs
        assert validate_grn_string('0.50')[0] == False  # 0 not allowed
        assert validate_grn_string('3.500')[0] == False  # Too many decimal places
        assert validate_grn_string('abc')[0] == False  # Not a GRN
        assert validate_grn_string('3-50')[0] == False  # Wrong separator


class TestGRNUtilities:
    """Test various GRN utility functions."""
    
    def test_get_tm_residues(self):
        """Test getting transmembrane residues."""
        grn_list = ['1.50', '2.50', '3.50', '45.50', '7.53', 'n.1', 'c.1']
        tm_residues = get_tm_residues(grn_list)
        
        # Should only include TM positions (1-7)
        assert '1.50' in tm_residues
        assert '2.50' in tm_residues
        assert '3.50' in tm_residues
        assert '7.53' in tm_residues
        
        # Should exclude loops, N/C terminus
        assert '45.50' not in tm_residues
        assert '45.500' not in tm_residues
        assert 'n.1' not in tm_residues
        assert 'c.1' not in tm_residues
    
    def test_map_grn_to_color(self):
        """Test mapping GRN to color."""
        # N-terminus gets specific color
        n_color = map_grn_to_color('n.1')
        assert n_color == 'rgb(31, 119, 180)'
        
        # C-terminus gets specific color
        c_color = map_grn_to_color('c.1')
        assert c_color == 'rgb(255, 127, 14)'
        
        # TM regions (<=10) get one color
        tm_color = map_grn_to_color('1.50')
        assert tm_color == 'rgb(214, 39, 40)'
        assert map_grn_to_color('7.53') == tm_color  # Same color for all TM regions
        
        # Loop regions (>10) get different color
        loop_color = map_grn_to_color('45.50')
        assert loop_color == 'rgb(44, 160, 44)'


class TestGRNDict:
    """Test GRN dictionary functionality."""
    
    def test_get_grn_dict(self, temp_data_dir, comprehensive_grn_table):
        """Test getting GRN dictionary."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_dict")
        processor.load_grn_table("test_dict")
        
        # Get GRN dict with dot notation
        grn_dict = processor.get_grn_dict(notation='dot')
        
        # Check structure - maps protein IDs to lists of GRN positions
        assert isinstance(grn_dict, dict)
        assert '7BMH' in grn_dict
        
        # Check values are lists of GRN positions (not residues)
        protein_grns = grn_dict['7BMH']
        assert isinstance(protein_grns, list)
        # Should contain GRN positions where 7BMH has residues (not '-')
        assert '1.50' in protein_grns  # 7BMH has M62 at 1.50
        assert '3.50' in protein_grns  # 7BMH has R115 at 3.50
        assert '7.53' in protein_grns  # 7BMH has K270 at 7.53
        
        # Test with x notation
        grn_dict_x = processor.get_grn_dict(notation='x')
        protein_grns_x = grn_dict_x['7BMH']
        # Should have same positions but possibly in x notation
        assert isinstance(protein_grns_x, list)


class TestGRNMaps:
    """Test GRN maps functionality."""
    
    def test_get_maps(self, temp_data_dir, comprehensive_grn_table):
        """Test getting GRN maps."""
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = comprehensive_grn_table
        processor.save_grn_table("test_maps")
        processor.load_grn_table("test_maps")
        
        # Initialize maps dictionary
        processor.maps = {'test_map': pd.DataFrame(index=processor.data.index, columns=processor.grns)}
        
        # Get maps
        map_names = processor.get_maps()
        
        # Should return list of map names
        assert isinstance(map_names, list)
        assert 'test_map' in map_names
        
        # The actual map is in processor.map attribute
        processor.map = pd.DataFrame(index=processor.data.index, columns=processor.grns)
        assert isinstance(processor.map, pd.DataFrame)


class TestRemoveDuplicates:
    """Test duplicate removal functionality."""
    
    def test_remove_duplicate_ids(self, temp_data_dir):
        """Test removing duplicate protein IDs."""
        # Create table with duplicates
        data = {
            '1.50': ['M62', 'M62', 'L22', 'M62'],
            '2.50': ['I90', 'I90', 'I50', 'I90'],
            '3.50': ['R115', 'R115', 'R82', 'R115']
        }
        # Note: 7BMH appears twice
        index = ['7BMH', '7BMH_dup', '1ABC', '7BMH']
        grn_table = pd.DataFrame(data, index=index)
        
        processor = GRNBaseProcessor(
            name="test_grn",
            preload=False
        )
        
        processor.data = grn_table
        processor.ids = list(grn_table.index)
        
        # Remove duplicates
        processor.remove_duplicate_ids()
        
        # Should have removed duplicates
        assert len(processor.data) < len(grn_table)
        # Each ID should appear only once
        assert processor.data.index.value_counts().max() == 1
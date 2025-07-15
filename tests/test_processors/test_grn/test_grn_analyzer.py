"""Test GRN Analyzer functionality."""

import pytest
import pandas as pd
from protos.processing.grn import GRNProcessor
from protos.processing.grn.grn_analyzer import GRNAnalyzer


class TestGRNAnalyzer:
    """Test GRN Analyzer with data management principles."""
    
    @pytest.fixture
    def grn_processor(self):
        """Create a GRN processor with test data."""
        processor = GRNProcessor(name="test_analyzer")
        
        # Create test GRN data
        test_data = pd.DataFrame({
            '1.50': ['M20', 'L20', 'M20'],
            '2.50': ['V50', 'I50', 'V50'],
            '3.50': ['T90', 'S90', 'T90'],
            '7.50': ['K216', 'K216', 'R216']
        }, index=['PROT1', 'PROT2', 'PROT3'])
        
        processor.data = test_data
        return processor
    
    def test_analyzer_initialization_with_processor(self, grn_processor):
        """Test initializing analyzer with processor."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        assert analyzer.data is not None
        assert len(analyzer.data) == 3
        assert analyzer.table_name == 'GRN default table'
        assert len(analyzer.features) == 3
        assert len(analyzer.map.columns) == 4
    
    def test_analyzer_initialization_with_data(self):
        """Test initializing analyzer with direct data."""
        test_data = pd.DataFrame({
            '1.50': ['M20'],
            '2.50': ['V50']
        }, index=['TEST1'])
        
        analyzer = GRNAnalyzer(grn_data=test_data, table_name='Test Analysis')
        
        assert analyzer.data is not None
        assert len(analyzer.data) == 1
        assert analyzer.table_name == 'Test Analysis'
    
    def test_analyzer_requires_data(self):
        """Test that analyzer requires either processor or data."""
        with pytest.raises(ValueError, match="Either grn_data or grn_processor"):
            GRNAnalyzer()
    
    def test_get_seq(self, grn_processor):
        """Test getting sequence from GRN data."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        seq = analyzer.get_seq('PROT1')
        assert seq == 'MVTK'  # M20 -> M, V50 -> V, T90 -> T, K216 -> K
    
    def test_get_seqs(self, grn_processor):
        """Test getting all sequences."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        seqs = analyzer.get_seqs()
        assert len(seqs) == 3
        assert seqs['PROT1'] == 'MVTK'  # M, V, T, K
        assert seqs['PROT2'] == 'LISK'  # L, I, S, K 
        assert seqs['PROT3'] == 'MVTR'  # M, V, T, R
    
    def test_get_lens(self, grn_processor):
        """Test getting sequence lengths."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        lens = analyzer.get_lens()
        assert len(lens) == 3
        assert all(l == 4 for l in lens)  # All sequences should be length 4 (4 positions)
    
    def test_get_interval(self, grn_processor):
        """Test getting specific GRN interval."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        # Get interval for helices 1 and 7
        interval_data = analyzer.get_interval(['1.50', '7.50'])
        
        assert interval_data.shape == (3, 2)
        assert '1.50' in interval_data.columns
        assert '7.50' in interval_data.columns
        assert '2.50' not in interval_data.columns
    
    def test_get_entries(self, grn_processor):
        """Test getting specific entries."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        entries = analyzer.get_entries(['PROT1', 'PROT3'])
        
        assert len(entries) == 2
        assert 'PROT1' in entries.index
        assert 'PROT3' in entries.index
        assert 'PROT2' not in entries.index
    
    def test_populate_map_features(self, grn_processor):
        """Test populating map features for amino acids."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        # Map cysteine occurrences
        analyzer.populate_map_features(mode='aminoacid', selection=['C'])
        
        # Should be all zeros since no cysteines in test data
        assert analyzer.map.sum().sum() == 0
        
        # Map lysine occurrences
        analyzer.populate_map_features(mode='aminoacid', selection=['K'])
        
        # Should have 2 lysines at position 7.50
        assert analyzer.map['7.50'].sum() == 2
        assert analyzer.map['1.50'].sum() == 0
    
    def test_populate_length_features(self, grn_processor):
        """Test populating length features."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        # Get overall lengths
        analyzer.populate_length_features()
        
        assert 'length' in analyzer.features.columns
        assert all(analyzer.features['length'] == 4)  # 4 positions
        
        # Get length for specific interval
        analyzer.populate_length_features(grn_interval=['1.50', '2.50'])
        
        assert 'length_1.50_2.50' in analyzer.features.columns
    
    def test_reset_data(self, grn_processor):
        """Test resetting data from processor."""
        analyzer = GRNAnalyzer(grn_processor=grn_processor)
        
        # Modify data
        analyzer.data = analyzer.data.head(1)
        assert len(analyzer.data) == 1
        
        # Reset
        analyzer.reset_data()
        assert len(analyzer.data) == 3  # Back to original
    
    def test_analyzer_with_real_grn_table(self):
        """Test analyzer with a real GRN table."""
        processor = GRNProcessor(name="test_real_analyzer")
        
        # Load a reference table if available
        try:
            ref_table = processor.load_reference_table('mo_ref')
            processor.data = ref_table
            
            analyzer = GRNAnalyzer(grn_processor=processor, table_name='MO Reference Analysis')
            
            # Basic checks
            assert len(analyzer) > 0
            assert '7.50' in analyzer.data.columns
            
            # Check sequence extraction
            first_id = analyzer.data.index[0]
            seq = analyzer.get_seq(first_id)
            assert len(seq) > 0
            
        except FileNotFoundError:
            pytest.skip("Reference table not available")
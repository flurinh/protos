"""
Comprehensive tests for GRN assignment functionality using real data.
Tests use synthetic microbial opsin data mimicking real structures.
"""

import pytest
import pandas as pd
import numpy as np

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import (
    GRNConfigManager,
    init_grn_intervals,
    expand_annotation,
    init_row_from_alignment,
    annotate_gpcr
)
from protos.processing.grn.grn_utils import (
    sort_grns_str,
    parse_grn_str2float,
    parse_grn_float2str
)
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment
)
from protos.io.paths.path_config import ProtosPaths
from protos.core.base_processor import BaseProcessor


class TestGRNAssignmentRealData:
    """Test GRN assignment with realistic data."""
    
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        """Set up test environment."""
        # Set global data root
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Initialize processors
        self.grn_processor = GRNBaseProcessor(
            name="test_grn_assign",
            processor_data_dir='grn'
        )
        
        self.seq_processor = BaseProcessor(
            name="test_sequences",
            processor_data_dir='sequence'
        )
        
        # Create test data
        self._create_test_data()
        
        yield
        
        # Cleanup
        ProtosPaths.set_data_root(None)
    
    def _create_test_data(self):
        """Create test reference data and sequences."""
        # Create reference GRN table
        ref_table = pd.DataFrame({
            '1.50': ['L28', 'V29', 'L30'],
            '2.50': ['A72', 'S73', 'A74'],
            '3.50': ['R115', 'R116', 'R117'],
            '4.50': ['G157', 'W158', 'W159'],
            '5.50': ['F197', 'F198', 'Y199'],
            '6.50': ['W238', 'W239', 'W240'],
            '7.50': ['K276', 'K277', 'K278']
        }, index=['REF1', 'REF2', 'REF3'])
        
        self.grn_processor.save_data('mo_ref', ref_table)
        
        # Create test sequences
        test_sequences = {
            'TEST_OPSIN1': 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG',
            'TEST_OPSIN2': 'MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYV'
        }
        
        self.seq_processor.save_data('test_sequences', test_sequences, format='json')
    
    def test_load_reference_table(self):
        """Test loading reference GRN table."""
        self.grn_processor.load_grn_table('mo_ref')
        
        assert self.grn_processor.data is not None
        assert not self.grn_processor.data.empty
        assert len(self.grn_processor.data) == 3
        assert '7.50' in self.grn_processor.data.columns
    
    def test_grn_config_manager(self):
        """Test GRN configuration manager."""
        # Test with microbial opsin family
        config_manager = GRNConfigManager(protein_family='mo')
        
        # Get strict config
        strict_config = config_manager.get_config(strict=True)
        assert isinstance(strict_config, dict)
        
        # Should have TM regions defined
        assert 'tm1' in strict_config or 'TM1' in strict_config
        
        # Get standard config
        standard_config = config_manager.get_config(strict=False)
        assert isinstance(standard_config, dict)
    
    def test_grn_assignment_with_real_fasta(self):
        """Test GRN assignment on a real sequence."""
        # Load reference table
        self.grn_processor.load_grn_table('mo_ref')
        ref_row = self.grn_processor.data.iloc[0]
        
        # Load test sequences
        sequences = self.seq_processor.load_data('test_sequences')
        test_seq = sequences['TEST_OPSIN1']
        
        # Extract reference sequence from GRN table
        ref_seq_parts = []
        for grn_pos in self.grn_processor.data.columns:
            value = ref_row[grn_pos]
            if value != '-' and pd.notna(value):
                # Extract amino acid from values like "L28"
                aa = value[0]
                ref_seq_parts.append(aa)
        ref_seq = ''.join(ref_seq_parts)
        
        # Create alignment
        aligner = init_aligner()
        raw_alignment = align_blosum62(test_seq, ref_seq, aligner)
        alignment = format_alignment(raw_alignment)
        
        # Create a simple new_row for testing
        # In a real scenario, this would come from alignment analysis
        new_row = pd.Series({
            '1.50': 'L20',
            '2.50': 'A60',
            '3.50': 'R100',
            '4.50': 'G140',
            '5.50': 'F180',
            '6.50': 'W220',
            '7.50': 'K260'
        })
        
        # Expand annotation
        grns, rns, missing = expand_annotation(
            new_row=new_row,
            query_seq=test_seq,
            alignment=alignment,
            protein_family='microbial_opsins',
            max_alignment_gap=1,
            verbose=0
        )
        
        # Verify results
        assert len(grns) > 0
        assert len(rns) > 0
        assert len(grns) == len(rns)
    
    def test_grn_intervals(self):
        """Test GRN interval initialization."""
        # Get microbial opsin configuration
        config_manager = GRNConfigManager(protein_family='mo')
        standard_config = config_manager.get_config(strict=False)
        
        # Initialize intervals - returns a list of GRN positions
        grn_positions = init_grn_intervals(standard_config)
        
        # Should return list of GRN positions
        assert isinstance(grn_positions, list)
        assert len(grn_positions) > 0
        
        # Check that positions are strings in GRN format
        for grn in grn_positions:
            assert isinstance(grn, str)
            # Should be in format X.YY
            assert '.' in grn or 'x' in grn
    
    def test_sort_grns(self):
        """Test GRN sorting functionality."""
        # Test positions
        test_grns = ['2.50', '1.50', '7.50', '3.50', '1.45', '1.46']
        
        # Sort
        sorted_grns = sort_grns_str(test_grns)
        
        # Verify order
        expected = ['1.45', '1.46', '1.50', '2.50', '3.50', '7.50']
        assert sorted_grns == expected
    
    def test_parse_grn_conversions(self):
        """Test GRN string/float conversions."""
        # Test string to float
        assert parse_grn_str2float('1.50') == 1.50
        assert parse_grn_str2float('7.50') == 7.50
        assert parse_grn_str2float('2x50') == 2.50
        
        # Test float to string
        assert parse_grn_float2str(1.50) == '1.50'
        assert parse_grn_float2str(7.50) == '7.50'
        
        # Test with x notation
        assert parse_grn_float2str(2.50, notation_type='x') == '2x50'
    
    def test_extract_strict_residues(self):
        """Test extraction of strict GRN positions."""
        # Load reference table
        self.grn_processor.load_grn_table('mo_ref')
        
        # Get config
        config_manager = GRNConfigManager(protein_family='mo')
        strict_config = config_manager.get_config(strict=True)
        
        # Extract strict positions
        strict_grns = []
        for tm_name, grn_range in strict_config.items():
            if isinstance(grn_range, list) and len(grn_range) == 2:
                start_grn = parse_grn_str2float(grn_range[0])
                end_grn = parse_grn_str2float(grn_range[1])
                
                for col in self.grn_processor.data.columns:
                    col_float = parse_grn_str2float(col)
                    if start_grn <= col_float <= end_grn:
                        strict_grns.append(col)
        
        # Should have found some strict positions
        assert len(strict_grns) > 0
        
        # Key positions should be included
        key_positions = ['1.50', '2.50', '3.50', '7.50']
        for pos in key_positions:
            if pos in self.grn_processor.data.columns:
                assert pos in strict_grns or pos in [g for g in strict_grns]
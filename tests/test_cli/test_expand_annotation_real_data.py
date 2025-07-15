"""
Tests for expand_annotation using synthetic opsin data.

Uses processor methods and dataset IDs instead of direct file access.
"""

import pytest
import pandas as pd

from protos.cli.grn.assign_grns import (
    expand_annotation,
    init_aligner,
    align_blosum62,
    format_alignment
)
from protos.processing.grn import GRNProcessor
from protos.core.base_processor import BaseProcessor
from protos.io.paths.path_config import ProtosPaths


class TestExpandAnnotationRealData:
    """Test expand_annotation with synthetic opsin data."""
    
    @pytest.fixture
    def setup_processors(self, tmp_path):
        """Set up GRN and sequence processors with test data."""
        # Set global data root to tmp directory
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Initialize processors
        grn_processor = GRNProcessor(
            name='test_expand'
        )
        
        from protos.processing.sequence import SequenceProcessor
        seq_processor = SequenceProcessor(
            name='test_expand'
        )
        
        # Create minimal test GRN reference data
        test_grn_ref = pd.DataFrame({
            '1.50': ['L21', 'L22', 'L23'],
            '2.50': ['A65', 'A66', 'A67'],
            '3.50': ['R107', 'R108', 'R109'],
            '4.50': ['G149', 'G150', 'G151'],
            '5.50': ['F189', 'F190', 'F191'],
            '6.50': ['W230', 'W231', 'W232'],
            '7.50': ['K268', 'K269', 'K270']
        }, index=['REF1', 'REF2', 'REF3'])
        
        grn_processor.save_data('test_ref', test_grn_ref)
        grn_processor.load_grn_table('test_ref')
        
        # Create test sequences
        test_sequences = {
            'TEST_SEQ1': 'MRPSGTAGAALLALLAALCPASRLEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYV',
            'TEST_SEQ2': 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN',
            'TEST_SEQ3': 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'
        }
        
        # Save sequences using the sequence processor
        seq_processor.save_data('test_sequences', test_sequences, format='json')
        
        return grn_processor, seq_processor
    
    def test_expand_annotation_single_sequence(self, setup_processors):
        """Test expand_annotation on a single sequence."""
        grn_processor, seq_processor = setup_processors
        
        # Load test sequences
        sequences = seq_processor.load_data('test_sequences')
        test_seq = sequences['TEST_SEQ1']
        
        # Get reference from GRN processor
        ref_row = grn_processor.data.loc['REF1']
        
        # Extract reference sequence from GRN table
        ref_seq_parts = []
        for grn, value in ref_row.items():
            if value != '-' and isinstance(value, str) and len(value) > 0:
                aa = value[0]  # First character is the amino acid
                ref_seq_parts.append(aa)
        ref_seq = ''.join(ref_seq_parts)
        
        # Create alignment
        aligner = init_aligner(open_gap_score=-22)
        raw_alignment = align_blosum62(test_seq, ref_seq, aligner)
        alignment = format_alignment(raw_alignment)
        
        # Create initial GRN assignments (simulating what would come from alignment)
        new_row = pd.Series({
            '1.50': 'L21',
            '2.50': 'A65',
            '3.50': 'R107',
            '4.50': 'G149',
            '5.50': 'F189',
            '6.50': 'W230',
            '7.50': 'K268'
        })
        
        # Run expand_annotation
        grns, rns, missing = expand_annotation(
            new_row=new_row,
            query_seq=test_seq,
            alignment=alignment,
            protein_family='microbial_opsins',
            max_alignment_gap=1,
            verbose=0
        )
        
        # Verify results
        assert len(grns) > 0, "No GRNs were assigned"
        assert len(rns) > 0, "No residue numbers were assigned"
        assert len(grns) == len(rns), "GRNs and residue numbers don't match"
        
        # Save results using processor - prepare GRN table properly
        result_df = pd.DataFrame({grn: [rn] for grn, rn in zip(grns, rns)}, index=['TEST_SEQ1'])
        grn_processor.data = result_df
        grn_processor.ids = result_df.index.tolist()
        grn_processor.grns = result_df.columns.tolist()
        grn_processor.save_grn_table("test_seq1_expanded")
        
        # Verify the results were saved correctly
        grn_processor.load_grn_table("test_seq1_expanded")
        loaded_result = grn_processor.data
        assert loaded_result is not None
        assert 'TEST_SEQ1' in loaded_result.index
    
    def test_expand_annotation_multiple_sequences(self, setup_processors):
        """Test expand_annotation on multiple sequences."""
        grn_processor, seq_processor = setup_processors
        
        # Load test sequences
        sequences = seq_processor.load_data('test_sequences')
        
        all_results = {}
        
        for seq_name, query_seq in sequences.items():
            # Get reference sequence
            ref_row = grn_processor.data.loc['REF1']
            ref_seq_parts = []
            for grn, value in ref_row.items():
                if value != '-' and isinstance(value, str) and len(value) > 0:
                    aa = value[0]
                    ref_seq_parts.append(aa)
            ref_seq = ''.join(ref_seq_parts)
            
            # Create alignment
            aligner = init_aligner(open_gap_score=-22)
            raw_alignment = align_blosum62(query_seq, ref_seq, aligner)
            alignment = format_alignment(raw_alignment)
            
            # Create initial assignments
            new_row = pd.Series({
                '1.50': f'L{20}',
                '2.50': f'A{60}',
                '3.50': f'R{100}',
                '4.50': f'G{140}',
                '5.50': f'F{180}',
                '6.50': f'W{220}',
                '7.50': f'K{min(260, len(query_seq)-1)}'
            })
            
            # Run expand_annotation
            grns, rns, missing = expand_annotation(
                new_row=new_row,
                query_seq=query_seq,
                alignment=alignment,
                protein_family='microbial_opsins',
                max_alignment_gap=1,
                verbose=0
            )
            
            all_results[seq_name] = {
                'grns': grns,
                'rns': rns,
                'missing': missing
            }
            
            # Basic validation
            assert len(grns) > 0, f"No GRNs assigned for {seq_name}"
            assert len(grns) == len(rns), f"Mismatch for {seq_name}"
        
        # Save combined results
        combined_data = {}
        for seq_name, result in all_results.items():
            grn_rn_dict = dict(zip(result['grns'], result['rns']))
            combined_data[seq_name] = grn_rn_dict
        
        df = pd.DataFrame.from_dict(combined_data, orient='index')
        df.fillna('-', inplace=True)
        grn_processor.save_data("multiple_sequences_expanded", df)
        
        # Verify the combined results were saved
        loaded_df = grn_processor.load_data("multiple_sequences_expanded")
        assert loaded_df is not None
        assert len(loaded_df) == len(sequences)  # Should have all sequences
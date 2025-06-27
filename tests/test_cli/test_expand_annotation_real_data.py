"""
Tests for expand_annotation using real microbial opsin data.
"""

import pytest
import pandas as pd
from pathlib import Path

from protos.cli.grn.assign_grns import (
    expand_annotation,
    init_aligner,
    align_blosum62,
    format_alignment
)
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.paths.path_config import ProtosPaths
from protos.io.fasta_utils import read_fasta


class TestExpandAnnotationRealData:
    """Test expand_annotation with real microbial opsin data."""
    
    @pytest.fixture
    def setup_paths(self, tmp_path):
        """Set up paths for testing."""
        # Set global data root
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Get test data paths
        test_data_dir = Path(__file__).parent.parent / "test-data"
        return {
            'fasta_file': test_data_dir / "sequence" / "fasta" / "opsin_sequences_from_yaml.fasta",
            'grn_ref_file': test_data_dir / "grn" / "ref" / "mo_ref.csv",
            'output_dir': tmp_path / "grn" / "tables"
        }
    
    @pytest.fixture
    def grn_processor(self, setup_paths, tmp_path):
        """Create GRN processor with real mo_ref data."""
        # Create directory structure
        grn_dir = tmp_path / "grn" / "tables"
        grn_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy mo_ref to expected location
        import shutil
        mo_ref_copy = grn_dir / "mo_ref.csv"
        shutil.copy(setup_paths['grn_ref_file'], mo_ref_copy)
        
        # Initialize processor
        processor = GRNBaseProcessor(
            name='mo_test',
            processor_data_dir='grn'
        )
        
        # Load the mo_ref table using dataset ID
        processor.load_grn_table('mo_ref')
        
        return processor
    
    @pytest.fixture
    def opsin_sequences(self, setup_paths):
        """Load opsin sequences from FASTA file."""
        sequences = read_fasta(str(setup_paths['fasta_file']))
        return sequences
    
    def test_expand_annotation_ar1(self, grn_processor, opsin_sequences):
        """Test expand_annotation on AR1_A sequence."""
        # Get AR1_A sequence
        ar1_seq = opsin_sequences['AR1_A']
        
        # Get a reference sequence from mo_ref (e.g., 7BMH)
        ref_row = grn_processor.data.loc['7BMH']
        
        # Create alignment between AR1 and reference
        aligner = init_aligner(open_gap_score=-22)
        
        # Get reference sequence from the GRN table
        ref_seq_parts = []
        for grn, aa in ref_row.items():
            if aa != '-':
                ref_seq_parts.append(aa)
        ref_seq = ''.join(ref_seq_parts)
        
        # Align sequences
        raw_alignment = align_blosum62(ar1_seq, ref_seq, aligner)
        alignment = format_alignment(raw_alignment)
        
        # Create new_row with initial GRN assignments from alignment
        # This would typically come from init_row_from_alignment
        # For testing, we'll create a minimal mapping
        new_row = pd.Series({
            '1.50': 'L21',  # L at position 21 maps to GRN 1.50
            '2.50': 'A65',  # A at position 65 maps to GRN 2.50
            '3.50': 'R107', # R at position 107 maps to GRN 3.50
            '4.50': 'G149', # G at position 149 maps to GRN 4.50
            '5.50': 'F189', # F at position 189 maps to GRN 5.50
            '6.50': 'W230', # W at position 230 maps to GRN 6.50
            '7.50': 'K268', # K at position 268 maps to GRN 7.50
        })
        
        # Run expand_annotation
        grns, rns, missing = expand_annotation(
            new_row=new_row,
            query_seq=ar1_seq,
            alignment=alignment,
            protein_family='microbial_opsins',
            max_alignment_gap=1,
            verbose=0
        )
        
        # Verify results
        assert len(grns) > 0, "No GRNs were assigned"
        assert len(rns) > 0, "No residue numbers were assigned"
        assert len(grns) == len(rns), "GRNs and residue numbers don't match"
        
        # Check that we have assignments for key positions
        grn_rn_dict = dict(zip(grns, rns))
        assert '1.50' in grn_rn_dict, "Missing assignment for TM1 position 1.50"
        assert '7.50' in grn_rn_dict, "Missing assignment for TM7 position 7.50"
        
        # Save the results
        output_file = grn_processor.data_root / "grn" / "tables" / "ar1_expanded.csv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Create DataFrame with results
        result_df = pd.DataFrame({grn: [rn] for grn, rn in zip(grns, rns)}, index=['AR1_A'])
        result_df.to_csv(output_file)
        
        return grns, rns, missing
    
    def test_expand_annotation_multiple_sequences(self, grn_processor, opsin_sequences, setup_paths):
        """Test expand_annotation on multiple opsin sequences."""
        # Test sequences
        test_sequences = ['AR1_A', 'AR2_A', 'AR3_A']
        
        all_results = {}
        
        for seq_name in test_sequences:
            if seq_name not in opsin_sequences:
                continue
                
            query_seq = opsin_sequences[seq_name]
            
            # Get reference sequence (7BMH)
            ref_row = grn_processor.data.loc['7BMH']
            ref_seq_parts = []
            for grn, aa in ref_row.items():
                if aa != '-':
                    ref_seq_parts.append(aa)
            ref_seq = ''.join(ref_seq_parts)
            
            # Create alignment
            aligner = init_aligner(open_gap_score=-22)
            raw_alignment = align_blosum62(query_seq, ref_seq, aligner)
            alignment = format_alignment(raw_alignment)
            
            # Create minimal new_row for testing
            # In real usage, this would come from proper alignment analysis
            new_row = pd.Series({
                '1.50': f'L{20}',  # Approximate positions
                '2.50': f'A{60}',
                '3.50': f'R{100}',
                '4.50': f'G{140}',
                '5.50': f'F{180}',
                '6.50': f'W{220}',
                '7.50': f'K{260}',
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
        output_file = setup_paths['output_dir'] / "multiple_sequences_expanded.csv"
        
        # Create combined DataFrame
        combined_data = {}
        for seq_name, result in all_results.items():
            grn_rn_dict = dict(zip(result['grns'], result['rns']))
            combined_data[seq_name] = grn_rn_dict
        
        # Convert to DataFrame with sequences as rows, GRNs as columns
        df = pd.DataFrame.from_dict(combined_data, orient='index')
        df.fillna('-', inplace=True)
        df.to_csv(output_file)
        
        return all_results
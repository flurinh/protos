"""
Tests for the GRN assignment workflow.

This module contains integration old_tests for the complete GRN assignment process,
from FASTA input to annotated GRN table output.
"""

import os
import pytest
import pandas as pd
from tempfile import NamedTemporaryFile

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string
)

# Test fixtures - protein sequences
@pytest.fixture
def test_sequences():
    """Provide test sequences for various protein families."""
    return {
        'gpcr_a': {
            'adrb2_human': "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL",
            'oprd_human': "MEPAPSAGAELQPPLFANASDAYPSACPSAGANASGPPGARSASSLALAIAITALYSAVCAVGLLGNVLVMFGIVRYTKMKTATNIYIFNLALADALATSTLPFQSAKYLMETWPFGELLCKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVKALDFRTPAKAKLINICIWVLASGVGVPIMVMAVTRPRDGAVVCMLQFPSPAWYKVLGSRAGSSLALAIAITALYSAVCAVGLLGNVLVMFGIVRYTKMKTATNIYIFNLALADALATSTLPFQSAKYLMETWPFGELLCKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVKALDFRTPAKAKLINICIWVLASGVGVPIMVMAVTRPRDGAVVCMLQFPSP"
        },
        'microbial_opsins': {
            'bacteriorhodopsin': "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD",
            'halorhodopsin': "MVGVKVVADGPAYAQSADISHGSVPLFIGLALVIAGIALAGVMSGQAASGTASALLSKLTNSLFTWIGALPGAGLSALLGKIRAVVDWLLVTLGALLVVLGLSGLAFEKRPASGALATLSGWAVSGPAATGAVGLVIAYFLGAGAMVWGPMMDAALAAGAGTVGGALAVFGYTSWGLVWLVLGEGAGAVVSGVVATVLWAAIGHSLILASVGAGGLEPGHHHHHH"
        }
    }

@pytest.fixture
def create_test_fasta(tmp_path):
    """Create a test FASTA file from sequences."""
    def _create_fasta(sequences, family):
        fasta_path = tmp_path / f"{family}_test.fasta"
        with open(fasta_path, "w") as f:
            for name, seq in sequences[family].items():
                f.write(f">{name}\n{seq}\n")
        return str(fasta_path)
    return _create_fasta

# Mock for the assign_grns_to_fasta function - we'll implement this later
def mock_assign_grns_to_fasta(
    fasta_path, 
    reference_dataset, 
    protein_family, 
    output_path, 
    num_cores=1
):
    """Mock implementation for testing."""
    # This is a placeholder - we'll implement the real function later
    # For now, just create a minimal output file for testing
    sequences = {}
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        current_id = None
        current_seq = []
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                    current_seq = []
                current_id = line[1:]
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    # Create a simple GRN table
    grn_data = {}
    for protein_id in sequences.keys():
        grn_data[protein_id] = {
            '1x50': 'L',
            '2x50': 'A',
            '3x50': 'R',
            '4x50': 'Y',
            '5x50': 'F',
            '6x50': 'W',
            '7x50': 'N',
            'n.5': 'M',
            'c.3': 'L',
            '12.003': 'G'
        }
    
    # Save to CSV
    df = pd.DataFrame.from_dict(grn_data, orient='index')
    df.to_csv(output_path, index=True)
    
    return {
        "success": True,
        "assigned_sequences": len(sequences),
        "reference_dataset": reference_dataset,
        "protein_family": protein_family
    }

class TestGRNAssignmentWorkflow:
    """Tests for the complete GRN assignment workflow."""
    
    def test_assignment_gpcr(self, test_sequences, create_test_fasta, tmp_path):
        """Test GRN assignment for GPCR sequences."""
        # Setup
        fasta_path = create_test_fasta(test_sequences, 'gpcr_a')
        output_path = str(tmp_path / "gpcr_test_grn.csv")
        
        # Execute
        result = mock_assign_grns_to_fasta(
            fasta_path=fasta_path,
            reference_dataset="ref",
            protein_family="gpcr_a",
            output_path=output_path,
            num_cores=1  # For testing
        )
        
        # Verify
        assert result["success"] is True
        assert result["assigned_sequences"] == len(test_sequences['gpcr_a'])
        
        # Load and validate output
        assert os.path.exists(output_path), f"Output file {output_path} not created"
        
        # Check if it can be loaded as a GRN table
        df = pd.read_csv(output_path, index_col=0)
        assert not df.empty, "GRN table is empty"
        
        # Validate GRN formats in column names
        for column in df.columns:
            is_valid, message = validate_grn_string(column)
            assert is_valid, f"Invalid GRN format: {column}, {message}"
    
    def test_grn_validation(self):
        """Test GRN format validation functions."""
        # Standard formats
        assert validate_grn_string('1x50')[0] is True
        assert validate_grn_string('n.10')[0] is True
        assert validate_grn_string('c.5')[0] is True
        assert validate_grn_string('12.003')[0] is True
        
        # Invalid formats
        assert validate_grn_string('9x50')[0] is False  # Invalid helix number
        assert validate_grn_string('1x100')[0] is False  # Position out of range
        assert validate_grn_string('12.1000')[0] is False  # Distance too long
        
        # Legacy formats that should normalize
        is_valid, message = validate_grn_string('12x05')
        assert is_valid is True
        assert "normalization" in message.lower()
        
    def test_grn_normalization(self):
        """Test GRN format normalization functions."""
        # Legacy loop format with x
        assert normalize_grn_format('12x05') == '12.005'
        
        # Legacy loop format without zero padding
        assert normalize_grn_format('12.5') == '12.005'
        
        # Standard GRN with dot instead of x
        assert normalize_grn_format('1.50') == '1x50'
        
        # Already normalized formats should remain unchanged
        assert normalize_grn_format('1x50') == '1x50'
        assert normalize_grn_format('n.10') == 'n.10'
        assert normalize_grn_format('c.5') == 'c.5'
        assert normalize_grn_format('12.003') == '12.003'
        
    def test_grn_round_trip(self):
        """Test round-trip conversion between string and float GRN formats."""
        test_cases = [
            '1x50',         # Standard TM
            'n.10',         # N-terminal
            'c.5',          # C-terminal
            '12.003',       # Loop between 1-2, closer to 1
            '65.011'        # Loop between 5-6, closer to 6
        ]
        
        for grn in test_cases:
            # Convert to float and back to string
            float_val = parse_grn_str2float(grn)
            round_trip = parse_grn_float2str(float_val)
            
            # Normalize both for comparison (accounts for format differences)
            assert normalize_grn_format(grn) == normalize_grn_format(round_trip), \
                f"Round-trip conversion failed for {grn}: {float_val} -> {round_trip}"
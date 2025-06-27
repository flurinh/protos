"""
Pytest configuration for sequence processing tests
"""

import pytest
import pandas as pd


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>seq1
MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG
>seq2
TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITIL
>seq3
MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKK
"""


@pytest.fixture
def mock_grn_table():
    """Create a mock GRN reference table."""
    # Simplified GRN table for testing
    data = {
        '1.50': ['M', 'L', 'T'],
        '2.50': ['V', 'I', 'A'],
        '3.50': ['T', 'S', 'T'],
        '7.50': ['K', 'K', 'R']
    }
    
    index = ['BR', 'HR', 'bPR']
    
    # Add sequence columns
    sequences = {
        'BR': list("MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY"),
        'HR': list("MKKPLGLLGAPGENAWVDMAAVTGVSAALGVAGGVLGVLATVGAAAVAADPLLARTTRPGEWICLAFALVLVLVGVLY"),
        'bPR': list("MALPDTFFDLVAADAERQWWLVVGILVAVIGTAFTALSVGNGFNFGKPTDDIFNVFKTVFEIVLGSALVEVIIGTLSY")
    }
    
    # Create full table
    df = pd.DataFrame(index=index)
    
    # Add some sequence positions
    for i in range(10):
        col_name = f'seq_{i+1}'
        df[col_name] = [sequences[idx][i] if i < len(sequences[idx]) else '-' for idx in index]
    
    # Add GRN positions
    for grn, residues in data.items():
        df[grn] = residues
    
    return df


# Skip markers for tests requiring external tools
def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line(
        "markers", "requires_mmseqs2: mark test as requiring MMseqs2"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )
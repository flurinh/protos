"""
Pytest configuration for sequence processing tests
"""

import pytest
import os
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent.parent.parent / "test-data"


@pytest.fixture(scope="session")
def project_root():
    """Path to project root."""
    return Path(__file__).parent.parent.parent.parent


@pytest.fixture(autouse=True)
def setup_environment():
    """Set up test environment."""
    # Try to import and use the mmseqs_utils module
    try:
        from protos.processing.sequence.mmseqs_utils import detect_mmseqs2
        mmseqs_path, use_wsl = detect_mmseqs2()
        # The detect function already sets environment variables
    except ImportError:
        # Fallback to manual detection if module not available
        import subprocess
        import platform
        
        is_windows = platform.system() == 'Windows'
        
        if is_windows:
            # On Windows, check if we can use WSL
            try:
                # Try WSL first
                result = subprocess.run(['wsl', 'which', 'mmseqs'], capture_output=True, text=True, shell=True)
                if result.returncode == 0:
                    mmseqs_path = result.stdout.strip()
                    os.environ["MMSEQS_PATH"] = mmseqs_path
                    os.environ["USE_WSL_MMSEQS"] = "true"
                else:
                    # Try common WSL paths
                    wsl_paths = [
                        '~/MMseqs2/build/bin/mmseqs',
                        '~/mmseqs2/bin/mmseqs',
                        '/usr/local/bin/mmseqs',
                        '/usr/bin/mmseqs'
                    ]
                    for path in wsl_paths:
                        result = subprocess.run(['wsl', 'test', '-f', path], capture_output=True, shell=True)
                        if result.returncode == 0:
                            os.environ["MMSEQS_PATH"] = path
                            os.environ["USE_WSL_MMSEQS"] = "true"
                            break
            except:
                pass
        else:
            # On Linux/Unix
            try:
                result = subprocess.run(['which', 'mmseqs'], capture_output=True, text=True)
                if result.returncode == 0:
                    mmseqs_path = result.stdout.strip()
                    os.environ["MMSEQS_PATH"] = mmseqs_path
            except:
                pass
    
    # Set paths for tests
    project_root = Path(__file__).parent.parent.parent.parent
    data_dir = project_root / "src" / "protos" / "reference_data"
    
    if data_dir.exists():
        os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
        os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())
    
    yield
    
    # Cleanup
    if "USE_WSL_MMSEQS" in os.environ:
        del os.environ["USE_WSL_MMSEQS"]
    if "MMSEQS_PATH" in os.environ:
        del os.environ["MMSEQS_PATH"]


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_content = """>seq1
MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG
>seq2
TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITIL
>seq3
MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKK
"""
    
    fasta_file = tmp_path / "test_sequences.fasta"
    fasta_file.write_text(fasta_content)
    
    return str(fasta_file)


@pytest.fixture
def mock_grn_table():
    """Create a mock GRN reference table."""
    import pandas as pd
    
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
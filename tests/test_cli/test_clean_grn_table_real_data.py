"""
Tests for clean_grn_table using real microbial opsin GRN data.
"""

import pytest
import pandas as pd
from pathlib import Path

from protos.cli.grn.clean_grn_table import clean_grn_table, process_table
from protos.io.paths.path_config import ProtosPaths


class TestCleanGRNTableRealData:
    """Test clean_grn_table with real mo_ref data."""
    
    @pytest.fixture
    def setup_paths(self, tmp_path):
        """Set up paths for testing."""
        # Set global data root
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Get test data paths
        test_data_dir = Path(__file__).parent.parent / "test-data"
        return {
            'grn_ref_file': test_data_dir / "grn" / "ref" / "mo_ref.csv",
            'output_dir': tmp_path / "grn" / "tables"
        }
    
    def test_clean_mo_ref_table(self, setup_paths):
        """Test cleaning the mo_ref table."""
        # Create output directory
        setup_paths['output_dir'].mkdir(parents=True, exist_ok=True)
        
        # Define output path
        output_file = setup_paths['output_dir'] / "mo_ref_cleaned.csv"
        
        # Process the table without using processor (direct file I/O)
        erroneous_report = process_table(
            str(setup_paths['grn_ref_file']),
            str(output_file),
            use_processor=False
        )
        
        # Check that the output file was created
        assert output_file.exists(), "Cleaned file was not created"
        
        # Load and verify the cleaned data
        cleaned_df = pd.read_csv(output_file, index_col=0)
        original_df = pd.read_csv(setup_paths['grn_ref_file'], index_col=0)
        
        # Basic validation
        assert len(cleaned_df) > 0, "Cleaned table is empty"
        assert len(cleaned_df.columns) == len(original_df.columns), "Column count changed"
        
        # The erroneous_report should contain any problematic sequences
        print(f"Erroneous sequences found: {list(erroneous_report.keys())}")
    
    def test_clean_with_processor(self, setup_paths, tmp_path):
        """Test cleaning using dataset ID approach."""
        # Create the necessary directory structure
        grn_dir = tmp_path / "grn"
        grn_dir.mkdir(exist_ok=True)
        
        # Copy the test file to the processor's expected location
        import shutil
        tables_dir = grn_dir / "tables"
        tables_dir.mkdir(exist_ok=True)
        
        # Copy test file to tables directory with a simple name
        test_file = tables_dir / "test_mo_ref.csv"
        shutil.copy(setup_paths['grn_ref_file'], test_file)
        
        # Output file
        output_file = tables_dir / "test_mo_ref_cleaned.csv"
        
        # Clean the table using dataset IDs
        erroneous_report = clean_grn_table(
            "test_mo_ref",  # Input dataset ID
            "test_mo_ref_cleaned"  # Output dataset ID
        )
        
        # Verify the output exists
        assert output_file.exists(), "Cleaned file was not created"
        
        # Load and check the result
        cleaned_df = pd.read_csv(output_file, index_col=0)
        assert len(cleaned_df) > 0, "Cleaned table is empty"
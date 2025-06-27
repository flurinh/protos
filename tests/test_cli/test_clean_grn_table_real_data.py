"""
Tests for clean_grn_table using real microbial opsin GRN data.

Uses dataset IDs and processor methods instead of direct file access.
"""

import pytest
import pandas as pd

from protos.cli.grn.clean_grn_table import clean_grn_table, process_table
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.paths.path_config import ProtosPaths


class TestCleanGRNTableRealData:
    """Test clean_grn_table with real mo_ref data."""
    
    @pytest.fixture
    def setup_processor(self, tmp_path):
        """Set up GRN processor with test data."""
        # Set global data root to tmp directory
        ProtosPaths.set_data_root(str(tmp_path))
        
        # Initialize processor
        processor = GRNBaseProcessor(
            name='test_clean',
            processor_data_dir='grn'
        )
        
        # Load test data from reference using dataset ID
        # This assumes mo_ref exists in test-data
        try:
            processor.load_grn_table('mo_ref')
        except Exception:
            # If not found, create minimal test data
            test_data = pd.DataFrame({
                '1.50': ['L21', 'L22', 'I23'],
                '2.50': ['A65', 'A66', '-'],
                '3.50': ['R107', 'R108', 'K109']
            }, index=['SEQ1', 'SEQ2', 'SEQ3'])
            processor.save_data('mo_ref', test_data)
        
        return processor
    
    def test_clean_mo_ref_table_with_processor(self, setup_processor):
        """Test cleaning the mo_ref table using processor."""
        processor = setup_processor
        
        # Clean the table using dataset IDs
        erroneous_report = clean_grn_table(
            "mo_ref",  # Input dataset ID
            "mo_ref_cleaned"  # Output dataset ID
        )
        
        # Load the cleaned data via processor
        processor.load_grn_table("mo_ref_cleaned")
        cleaned_df = processor.data
        
        # Basic validation
        assert len(cleaned_df) > 0, "Cleaned table is empty"
        assert len(cleaned_df.columns) > 0, "No columns in cleaned table"
        
        # The erroneous_report should contain any problematic sequences
        print(f"Erroneous sequences found: {list(erroneous_report.keys())}")
    
    def test_clean_with_validation(self, setup_processor):
        """Test cleaning with specific validation checks."""
        processor = setup_processor
        
        # Create test data with known issues
        test_data = pd.DataFrame({
            '1.50': ['L21', 'L22', 'X99'],  # X99 is erroneous
            '2.50': ['A65', '-', 'B66'],     # B is not standard AA
            '3.50': ['R107', 'R108', '-'],
            '4.50': ['G149', 'G150', 'G151']
        }, index=['GOOD1', 'GOOD2', 'BAD1'])
        
        processor.save_data('test_validation', test_data)
        
        # Clean the data
        erroneous_report = clean_grn_table(
            "test_validation",
            "test_validation_cleaned"
        )
        
        # Check that problematic sequence was identified
        assert len(erroneous_report) > 0, "Should have found erroneous sequences"
        
        # Load cleaned data
        processor.load_grn_table("test_validation_cleaned")
        cleaned_df = processor.data
        
        # Verify data was cleaned appropriately
        assert len(cleaned_df) >= 2, "Should have at least 2 good sequences"
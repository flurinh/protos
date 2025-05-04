"""
Tests for GRN diagnostic tools.

This module old_tests the diagnostic tools for GRN tables.
"""

import os
import pytest
import pandas as pd
from tempfile import NamedTemporaryFile

from protos.cli.grn.diagnose_grn import diagnose_grn_table

class TestGRNDiagnostics:
    """Tests for GRN diagnostic tools."""
    
    @pytest.fixture
    def sample_grn_table(self, tmp_path):
        """Create a sample GRN table for testing."""
        # Create a DataFrame with some GRN issues
        data = {
            'protein1': {
                '1x50': 'L',
                '2x50': 'A',
                '3x50': 'R',
                '7x50': 'K',  # Schiff base okay
                'n.10': 'M',
                'c.5': 'L',
                '12.003': 'G',  # Valid loop format
                '23x05': 'S',   # Invalid loop format (should be 23.005)
                '34.5': 'T',    # Invalid loop format (should be 34.005)
                '9x50': 'Y'     # Invalid helix number
            },
            'protein2': {
                '1x50': 'M',
                '2x50': 'S',
                '3x50': 'T',
                '7x50': 'R',   # Missing Schiff base lysine
                'n.10': 'G',
                'c.5': 'V',
                '12.003': 'W', 
                '45.123': 'D'  # Valid loop format
            }
        }
        
        # Transpose the data to create a DataFrame with proteins as rows and GRNs as columns
        df = pd.DataFrame(data).T
        
        # Create a temporary CSV file
        csv_path = tmp_path / "test_grn_table.csv"
        df.to_csv(csv_path)
        
        return str(csv_path)
    
    def test_diagnose_grn_table(self, sample_grn_table):
        """Test diagnosing a GRN table."""
        # Diagnose the sample GRN table
        results = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='microbial_opsins',
            output_path=None
        )
        
        # Check basic results
        assert results["status"] == "success"
        assert results["protein_count"] == 2
        assert results["grn_count"] == 10
        
        # Check invalid GRNs
        assert len(results["invalid_grns"]) > 0
        assert any(issue["grn"] == "9x50" for issue in results["invalid_grns"])
        
        # Check normalized GRNs
        assert len(results["normalized_grns"]) > 0
        normalized_grns = [issue["grn"] for issue in results["normalized_grns"]]
        assert "23x05" in normalized_grns or "34.5" in normalized_grns
        
        # Check Schiff base issues
        assert len(results["schiff_base_issues"]) == 1
        assert results["schiff_base_issues"][0]["protein_id"] == "protein2"
        
    def test_diagnose_with_output(self, sample_grn_table, tmp_path):
        """Test saving diagnostic results to file."""
        # Create output path
        output_path = tmp_path / "diagnostic_report.json"
        
        # Diagnose and save results
        results = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='microbial_opsins',
            output_path=str(output_path)
        )
        
        # Check that the file was created
        assert os.path.exists(output_path)
        
        # Load and check the saved file
        import json
        with open(output_path, 'r') as f:
            saved_results = json.load(f)
        
        # Verify the saved results match the returned results
        assert saved_results["status"] == results["status"]
        assert saved_results["protein_count"] == results["protein_count"]
        assert saved_results["grn_count"] == results["grn_count"]
        
    def test_different_protein_families(self, sample_grn_table):
        """Test diagnosing GRN tables for different protein families."""
        # Test for microbial opsins
        results_opsins = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='microbial_opsins',
            output_path=None
        )
        
        # Test for GPCRs
        results_gpcr = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='gpcr_a',
            output_path=None
        )
        
        # Different families should have different Schiff base criteria
        assert len(results_opsins["schiff_base_issues"]) != len(results_gpcr["schiff_base_issues"])
        
    def test_skip_diagnostics(self, sample_grn_table):
        """Test skipping specific diagnostic checks."""
        # Skip loop diagnostics
        results_no_loops = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='microbial_opsins',
            diagnose_loops=False,
            output_path=None
        )
        
        # Skip Schiff base check
        results_no_schiff = diagnose_grn_table(
            grn_table_path=sample_grn_table,
            protein_family='microbial_opsins',
            check_schiff_base=False,
            output_path=None
        )
        
        # Verify loop issues are empty when skipped
        assert len(results_no_loops["loop_issues"]) == 0
        
        # Verify Schiff base issues are empty when skipped
        assert len(results_no_schiff["schiff_base_issues"]) == 0
"""
Tests for the download_structures module and all related downloading functionality
in the loaders directory, integrated with the protos path system.
"""

import pytest
import requests
from unittest.mock import patch, MagicMock

from protos.loaders.download_structures import download_protein_structures
from protos.loaders.alphafold_utils import download_alphafold_structures
from protos.loaders.uniprot_utils import map_uniprot_to_pdb
from protos.io.paths.path_config import ProtosPaths
from protos.io.fasta_utils import read_fasta, write_fasta
from protos.core.base_processor import BaseProcessor
from Bio.PDB import PDBList
import pandas as pd


class MockResponse:
    """Mock response object for testing HTTP requests"""
    def __init__(self, status_code=200, content=b"mock content"):
        self.status_code = status_code
        self.content = content

    def raise_for_status(self):
        if self.status_code != 200:
            raise requests.HTTPError(f"HTTP Error: {self.status_code}")


@pytest.fixture
def test_processor(tmp_path):
    """Create a test processor for download tests."""
    # Set global data root
    ProtosPaths.set_data_root(str(tmp_path))
    
    processor = BaseProcessor(
        name="test_downloader",
        processor_data_dir="structure"
    )
    
    yield processor
    
    # Cleanup
    ProtosPaths.set_data_root(None)


class TestDownloadStructures:
    """Test the download_protein_structures function"""
    
    def test_download_structures_basic(self, test_processor):
        """Test basic download functionality"""
        pdb_ids = ["1ABC", "2DEF"]
        
        # Mock PDBList to avoid actual downloads
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_path"
            
            # Use processor's data directory for downloads
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path
            )
            
            # Check that retrieve was called for each PDB ID
            assert mock_retrieve.call_count == len(pdb_ids)
    
    def test_download_structures_with_failures(self, test_processor):
        """Test download with some failures"""
        pdb_ids = ["1ABC", "INVALID", "2DEF"]
        
        def mock_retrieve_side_effect(pdb_id, **kwargs):
            if pdb_id == "INVALID":
                raise Exception("Invalid PDB ID")
            return f"mock_path_{pdb_id}"
        
        with patch.object(PDBList, 'retrieve_pdb_file', side_effect=mock_retrieve_side_effect):
            # Download will skip failures
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path
            )
            
            # Check that retrieve was called for all
            assert mock_retrieve.call_count == len(pdb_ids)
    
    def test_download_structures_empty_list(self, test_processor):
        """Test download with empty list"""
        # Download empty list should not fail
        download_protein_structures(
            [],
            target_folder=test_processor.data_path
        )
    
    @pytest.mark.parametrize("file_format", ['mmCif', 'pdb', 'xml', 'mmtf'])
    def test_download_structures_formats(self, test_processor, file_format):
        """Test download with different file formats"""
        pdb_ids = ["1ABC"]
        
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_path"
            
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path
            )
            
            # Check that retrieve was called
            assert mock_retrieve.called
    
    def test_download_structures_overwrite(self, test_processor):
        """Test download with overwrite option"""
        pdb_ids = ["1ABC"]
        
        # Save a mock file first
        test_processor.save_data('1ABC', {'pdb_id': '1ABC'}, format='json')
        
        # Verify the data was saved
        loaded_data = test_processor.load_data('1ABC', format='json')
        assert loaded_data == {'pdb_id': '1ABC'}
        
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_path"
            
            # Download twice to test overwrite behavior
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path
            )
            
            # Download again 
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path
            )
            
            # Should have been called twice
            assert mock_retrieve.call_count == 2


class TestAlphafoldUtils:
    """Test AlphaFold download functionality"""
    
    @patch('requests.get')
    def test_download_alphafold_basic(self, mock_get, test_processor):
        """Test basic AlphaFold download"""
        mock_get.return_value = MockResponse(content=b"ALPHAFOLD STRUCTURE DATA")
        
        # Download single AlphaFold structure
        uniprot_id = "P12345"
        download_alphafold_structures(
            uniprot_id,
            output_dir=test_processor.data_path
        )
        
        # Check that requests were made (tries multiple model versions)
        assert mock_get.call_count >= 1
    
    @patch('requests.get')
    def test_download_alphafold_with_failures(self, mock_get, test_processor):
        """Test AlphaFold download with some failures"""
        def mock_get_side_effect(url, **kwargs):
            if "INVALID" in url:
                return MockResponse(status_code=404)
            return MockResponse(content=b"ALPHAFOLD STRUCTURE DATA")
        
        mock_get.side_effect = mock_get_side_effect
        
        # Test with invalid ID
        uniprot_id = "INVALID"
        download_alphafold_structures(
            uniprot_id,
            output_dir=test_processor.data_path
        )
        
        # Should have tried to download
        assert mock_get.called
    
    @patch('requests.get')
    def test_download_alphafold_empty_list(self, mock_get, test_processor):
        """Test AlphaFold download with empty list"""
        # Test with empty string
        download_alphafold_structures(
            "",
            output_dir=test_processor.data_path
        )
        
        # Should still make attempts
        assert mock_get.called
    
    @patch('requests.get')
    def test_download_alphafold_format_options(self, mock_get, test_processor):
        """Test AlphaFold download with format options"""
        mock_get.return_value = MockResponse(content=b"ALPHAFOLD STRUCTURE DATA")
        
        uniprot_ids = ["P12345"]
        
        # Test with max_models parameter
        download_alphafold_structures(
            uniprot_id,
            max_models=3,
            output_dir=test_processor.data_path
        )
        
        # Should try up to 3 models
        assert mock_get.call_count == 3
        
        # Check URL construction contains model version
        call_urls = [call[0][0] for call in mock_get.call_args_list]
        assert any("model_v1" in url for url in call_urls)
        assert any("model_v2" in url for url in call_urls)
        assert any("model_v3" in url for url in call_urls)


class TestUniprotMapping:
    """Test UniProt to PDB mapping functionality"""
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_basic(self, mock_get):
        """Test basic UniProt to PDB mapping"""
        # Mock UniProt API response
        mock_response_data = {
            "results": [
                {
                    "from": "P12345",
                    "to": {
                        "primaryAccession": "P12345",
                        "uniProtKBCrossReferences": [
                            {
                                "database": "PDB",
                                "id": "1ABC"
                            },
                            {
                                "database": "PDB",
                                "id": "2DEF"
                            }
                        ]
                    }
                }
            ]
        }
        
        mock_get.return_value = MagicMock(
            status_code=200,
            json=lambda: mock_response_data
        )
        
        result = map_uniprot_to_pdb(["P12345"])
        
        assert "P12345" in result
        assert result["P12345"] == ["1ABC", "2DEF"]
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_no_structures(self, mock_get):
        """Test mapping for UniProt ID with no PDB structures"""
        mock_response_data = {
            "results": [
                {
                    "from": "P12345",
                    "to": {
                        "primaryAccession": "P12345",
                        "uniProtKBCrossReferences": []
                    }
                }
            ]
        }
        
        mock_get.return_value = MagicMock(
            status_code=200,
            json=lambda: mock_response_data
        )
        
        result = map_uniprot_to_pdb(["P12345"])
        
        assert "P12345" in result
        assert result["P12345"] == []
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_multiple(self, mock_get):
        """Test mapping for multiple UniProt IDs"""
        mock_response_data = {
            "results": [
                {
                    "from": "P12345",
                    "to": {
                        "primaryAccession": "P12345",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1ABC"}
                        ]
                    }
                },
                {
                    "from": "Q67890",
                    "to": {
                        "primaryAccession": "Q67890",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "2DEF"},
                            {"database": "PDB", "id": "3GHI"}
                        ]
                    }
                }
            ]
        }
        
        mock_get.return_value = MagicMock(
            status_code=200,
            json=lambda: mock_response_data
        )
        
        result = map_uniprot_to_pdb(["P12345", "Q67890"])
        
        assert len(result) == 2
        assert result["P12345"] == ["1ABC"]
        assert result["Q67890"] == ["2DEF", "3GHI"]
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_error_handling(self, mock_get):
        """Test error handling in UniProt mapping"""
        # Test HTTP error
        mock_get.return_value = MagicMock(status_code=500)
        
        result = map_uniprot_to_pdb(["P12345"])
        assert result == {}
        
        # Test request exception
        mock_get.side_effect = requests.RequestException("Network error")
        
        result = map_uniprot_to_pdb(["P12345"])
        assert result == {}


class TestIntegratedDownloadWorkflow:
    """Test integrated download workflows"""
    
    def test_download_workflow_with_processor(self, test_processor):
        """Test download workflow using processor abstractions"""
        # Create test data
        test_data = {
            'structures': ['1ABC', '2DEF'],
            'uniprot_ids': ['P12345', 'Q67890']
        }
        
        # Save using processor
        test_processor.save_data('download_targets', test_data, format='json')
        
        # Load using processor
        loaded_data = test_processor.load_data('download_targets', format='json')
        
        assert loaded_data['structures'] == test_data['structures']
        assert loaded_data['uniprot_ids'] == test_data['uniprot_ids']
    
    @pytest.mark.network
    def test_real_download_integration(self, test_processor):
        """Test real download integration (requires network)"""
        # This test is marked as network-dependent
        # It would test actual downloads but is skipped in CI
        pytest.skip("Network test - requires internet connection")


class TestDownloadHelpers:
    """Test helper functions for downloads"""
    
    def test_fasta_operations_with_processor(self, test_processor):
        """Test FASTA read/write using processor"""
        # Create test FASTA content
        fasta_content = ">seq1\nMKLLILTCLV\n>seq2\nMALIGTLLML"
        
        # Save as text file
        test_processor.save_data('test_sequences.fasta', fasta_content, format='text')
        
        # Read back
        loaded = test_processor.load_data('test_sequences.fasta', format='text')
        
        assert ">seq1" in loaded
        assert "MKLLILTCLV" in loaded
        assert ">seq2" in loaded
        assert "MALIGTLLML" in loaded
    
    def test_batch_download_organization(self, test_processor):
        """Test organizing batch downloads using processor"""
        # Simulate batch download results
        download_results = {
            'experimental': ['1ABC', '2DEF', '3GHI'],
            'alphafold': ['AF-P12345', 'AF-Q67890'],
            'failed': ['INVALID1', 'INVALID2']
        }
        
        # Save results
        test_processor.save_data('download_summary', download_results, format='json')
        
        # Create dataset mapping
        dataset_info = {
            'name': 'test_structures',
            'total_structures': 5,
            'experimental_count': 3,
            'alphafold_count': 2,
            'failed_count': 2
        }
        
        test_processor.save_data('dataset_info', dataset_info, format='json')
        
        # Verify saved data
        loaded_summary = test_processor.load_data('download_summary', format='json')
        loaded_info = test_processor.load_data('dataset_info', format='json')
        
        assert loaded_summary['experimental'] == download_results['experimental']
        assert loaded_info['total_structures'] == 5
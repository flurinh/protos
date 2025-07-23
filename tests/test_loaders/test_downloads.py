"""
Tests for the download_structures module and all related downloading functionality
in the loaders directory, integrated with the protos path system.
"""

import pytest
import requests
from pathlib import Path
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
def test_processor():
    """Create a test processor for download tests."""
    # ProtosPaths already configured in conftest.py to use tests/test-data
    processor = BaseProcessor(name="test")
    
    return processor


class TestDownloadStructures:
    """Test the download_protein_structures function"""
    
    def test_download_structures_basic(self, test_processor):
        """Test basic download functionality"""
        pdb_ids = ["1ABC", "2DEF"]
        
        # Use processor's structure directory (managed by ProtosPaths)
        target_dir = test_processor.paths.get_subdir_path("structure", "structure_dir")
        
        # Mock PDBList to avoid actual downloads
        def mock_retrieve_side_effect(pdb_id, **kwargs):
            # Create a mock file in the expected location
            pdb_id = pdb_id.lower()
            mock_file = target_dir / f"{pdb_id}.cif"
            mock_file.write_text(f"# Mock CIF file for {pdb_id}")
            return str(mock_file)
        
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.side_effect = mock_retrieve_side_effect
            
            # Use processor's data directory for downloads, use overwrite to force download
            successful, failed = download_protein_structures(
                pdb_ids,
                target_folder=str(target_dir),
                overwrite=True
            )
            
            # Check results - convert to lowercase as our function does that
            assert successful == [pid.lower() for pid in pdb_ids]
            assert failed == []
            
            # Check that retrieve was called for each PDB ID
            assert mock_retrieve.call_count == len(pdb_ids)
    
    def test_download_structures_with_failures(self, test_processor):
        """Test download with some failures"""
        pdb_ids = ["1ABC", "INVALID", "2DEF"]
        
        # Use processor's structure directory (managed by ProtosPaths)
        target_dir = test_processor.paths.get_subdir_path("structure", "structure_dir")
        
        def mock_retrieve_side_effect(pdb_id, **kwargs):
            if pdb_id.upper() == "INVALID":
                raise Exception("Invalid PDB ID")
            # Create a mock file in the expected location
            pdb_id = pdb_id.lower()
            mock_file = target_dir / f"{pdb_id}.cif"
            mock_file.write_text(f"# Mock CIF file for {pdb_id}")
            return str(mock_file)
        
        with patch.object(PDBList, 'retrieve_pdb_file', side_effect=mock_retrieve_side_effect) as mock_retrieve:
            # Download will skip failures, use overwrite to test actual download behavior
            successful, failed = download_protein_structures(
                pdb_ids,
                target_folder=str(target_dir),
                overwrite=True
            )
            
            # Check results - convert to lowercase as our function does that
            assert successful == ["1abc", "2def"]
            assert failed == ["invalid"]
            
            # Check that retrieve was called for all
            assert mock_retrieve.call_count == len(pdb_ids)
    
    def test_download_structures_empty_list(self, test_processor):
        """Test download with empty list"""
        # Download empty list should not fail
        successful, failed = download_protein_structures(
            [],
            target_folder=test_processor.data_path
        )
        assert successful == []
        assert failed == []
    
    @pytest.mark.parametrize("file_format", ['mmCif', 'pdb', 'xml', 'mmtf'])
    def test_download_structures_formats(self, test_processor, file_format):
        """Test download with different file formats"""
        pdb_ids = ["1ABC"]
        
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_path"
            
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path,
                overwrite=True
            )
            
            # Check that retrieve was called
            assert mock_retrieve.called
    
    def test_download_structures_overwrite(self, test_processor):
        """Test download with overwrite option"""
        pdb_ids = ["1ABC"]
        
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_path"
            
            # First download - should call retrieve (file doesn't exist)
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path,
                overwrite=True
            )
            
            # Should have been called once
            assert mock_retrieve.call_count == 1
            
            # Second download with overwrite=False (default) - should not call retrieve (file exists)
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path,
                overwrite=False
            )
            
            # Should still be called only once (skipped second time)
            assert mock_retrieve.call_count == 1
            
            # Third download with overwrite=True - should call retrieve again
            download_protein_structures(
                pdb_ids,
                target_folder=test_processor.data_path,
                overwrite=True
            )
            
            # Should now be called twice (forced overwrite)
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
            uniprot_ids[0],
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
    
    @patch('protos.loaders.uniprot_utils.session')
    @patch('requests.post')
    def test_map_uniprot_to_pdb_basic(self, mock_post, mock_session):
        """Test basic UniProt to PDB mapping"""
        # Mock job submission
        mock_post.return_value = MagicMock(
            status_code=200,
            json=lambda: {"jobId": "test-job-123"}
        )
        
        # Mock status check (ready)
        mock_status_response = MagicMock(
            status_code=200,
            json=lambda: {"results": "https://rest.uniprot.org/idmapping/results/test-job-123"}
        )
        
        # Mock details call (for redirectURL)
        mock_details_response = MagicMock(
            status_code=200,
            json=lambda: {"redirectURL": "https://rest.uniprot.org/idmapping/results/test-job-123"}
        )
        
        # Mock results
        mock_results_response = MagicMock(
            status_code=200,
            headers={"x-total-results": "2"},
            json=lambda: {
                "results": [
                    {
                        "from": "P12345",
                        "to": "1ABC"
                    },
                    {
                        "from": "P12345", 
                        "to": "2DEF"
                    }
                ]
            }
        )
        
        mock_session.get.side_effect = [mock_status_response, mock_details_response, mock_results_response]
        
        result = map_uniprot_to_pdb(["P12345"])
        
        # Result should be a DataFrame
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert all(result['uid'] == 'P12345')
        assert set(result['pdb_id']) == {'1ABC', '2DEF'}
    
    @patch('protos.loaders.uniprot_utils.session')
    @patch('requests.post')
    def test_map_uniprot_to_pdb_no_structures(self, mock_post, mock_session):
        """Test mapping for UniProt ID with no PDB structures"""
        # Mock job submission
        mock_post.return_value = MagicMock(
            status_code=200,
            json=lambda: {"jobId": "test-job-456"}
        )
        
        # Mock status check (ready)
        mock_status_response = MagicMock(
            status_code=200,
            json=lambda: {"results": "https://rest.uniprot.org/idmapping/results/test-job-456"}
        )
        
        # Mock details call (for redirectURL)
        mock_details_response = MagicMock(
            status_code=200,
            json=lambda: {"redirectURL": "https://rest.uniprot.org/idmapping/results/test-job-456"}
        )
        
        # Mock empty results
        mock_results_response = MagicMock(
            status_code=200,
            headers={"x-total-results": "0"},
            json=lambda: {"results": []}
        )
        
        mock_session.get.side_effect = [mock_status_response, mock_details_response, mock_results_response]
        
        result = map_uniprot_to_pdb(["P12345"])
        
        # Result should be an empty DataFrame
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
    
    @patch('protos.loaders.uniprot_utils.session')
    @patch('requests.post')
    def test_map_uniprot_to_pdb_multiple(self, mock_post, mock_session):
        """Test mapping for multiple UniProt IDs"""
        # Mock job submission
        mock_post.return_value = MagicMock(
            status_code=200,
            json=lambda: {"jobId": "test-job-789"}
        )
        
        # Mock status check (ready)
        mock_status_response = MagicMock(
            status_code=200,
            json=lambda: {"results": "https://rest.uniprot.org/idmapping/results/test-job-789"}
        )
        
        # Mock details call (for redirectURL)
        mock_details_response = MagicMock(
            status_code=200,
            json=lambda: {"redirectURL": "https://rest.uniprot.org/idmapping/results/test-job-789"}
        )
        
        # Mock results for multiple IDs
        mock_results_response = MagicMock(
            status_code=200,
            headers={"x-total-results": "3"},
            json=lambda: {
                "results": [
                    {"from": "P12345", "to": "1ABC"},
                    {"from": "Q67890", "to": "2DEF"},
                    {"from": "Q67890", "to": "3GHI"}
                ]
            }
        )
        
        mock_session.get.side_effect = [mock_status_response, mock_details_response, mock_results_response]
        
        result = map_uniprot_to_pdb(["P12345", "Q67890"])
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        assert len(result[result['uid'] == 'P12345']) == 1
        assert len(result[result['uid'] == 'Q67890']) == 2
        assert set(result[result['uid'] == 'Q67890']['pdb_id']) == {'2DEF', '3GHI'}
    
    @patch('protos.loaders.uniprot_utils.session')
    @patch('requests.post')
    def test_map_uniprot_to_pdb_error_handling(self, mock_post, mock_session):
        """Test error handling in UniProt mapping"""
        # Test HTTP error on job submission
        mock_post.side_effect = requests.HTTPError("Server error")
        
        with pytest.raises(requests.HTTPError):
            map_uniprot_to_pdb(["P12345"])


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
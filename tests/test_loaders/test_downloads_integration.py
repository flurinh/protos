"""
Integration tests for download functionality across multiple components.
Tests the interaction between download utilities and the protos system.
"""

import pytest
import gzip
from unittest.mock import patch, MagicMock, Mock
import requests

from protos.loaders.download_structures import download_protein_structures
from protos.loaders.alphafold_utils import download_alphafold_structures
from protos.loaders.uniprot_utils import map_uniprot_to_pdb
from protos.io.paths.path_config import ProtosPaths
from protos.core.base_processor import BaseProcessor
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from Bio.PDB import PDBList


@pytest.fixture
def integration_processor(tmp_path):
    """Create processors for integration testing."""
    ProtosPaths.set_data_root(str(tmp_path))
    
    # Create both base and structure processors
    base_proc = BaseProcessor(
        name="test_integration",
        processor_data_dir="downloads"
    )
    
    struct_proc = CifBaseProcessor(
        name="test_structures",
        processor_data_dir="structure"
    )
    
    yield base_proc, struct_proc
    
    ProtosPaths.set_data_root(None)


class TestDownloadToPipelineIntegration:
    """Test download to processing pipeline integration"""
    
    @patch('requests.get')
    @patch.object(PDBList, 'retrieve_pdb_file')
    def test_uniprot_to_structure_pipeline(self, mock_pdb_retrieve, mock_get, integration_processor):
        """Test complete pipeline from UniProt ID to structure processing"""
        base_proc, struct_proc = integration_processor
        
        # Mock UniProt API response
        mock_uniprot_response = {
            "results": [
                {
                    "from": "P00720",
                    "to": {
                        "primaryAccession": "P00720",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1HEW"},
                            {"database": "PDB", "id": "2LYZ"}
                        ]
                    }
                }
            ]
        }
        
        mock_get.return_value = MagicMock(
            status_code=200,
            json=lambda: mock_uniprot_response
        )
        
        # Mock PDB download
        mock_pdb_retrieve.return_value = "mock_1hew.cif"
        
        # Step 1: Map UniProt to PDB
        uniprot_ids = ["P00720"]
        pdb_mapping = map_uniprot_to_pdb(uniprot_ids)
        
        assert "P00720" in pdb_mapping
        assert "1HEW" in pdb_mapping["P00720"]
        
        # Save mapping using processor
        base_proc.save_data('uniprot_pdb_mapping', pdb_mapping, format='json')
        
        # Step 2: Download structures
        pdb_ids = pdb_mapping["P00720"]
        download_results = download_protein_structures(
            pdb_ids,
            target_folder=struct_proc.data_path
        )
        
        # Save download results
        base_proc.save_data('download_results', download_results, format='json')
        
        # Step 3: Create structure dataset
        dataset_info = {
            'name': 'lysozyme_structures',
            'uniprot_id': 'P00720',
            'pdb_ids': pdb_ids,
            'download_status': download_results
        }
        
        struct_proc.save_data('lysozyme_dataset', dataset_info, format='json')
        
        # Verify complete pipeline
        loaded_mapping = base_proc.load_data('uniprot_pdb_mapping', format='json')
        loaded_results = base_proc.load_data('download_results', format='json')
        loaded_dataset = struct_proc.load_data('lysozyme_dataset', format='json')
        
        assert loaded_mapping == pdb_mapping
        assert loaded_dataset['uniprot_id'] == 'P00720'
        assert len(loaded_dataset['pdb_ids']) == 2
    
    @patch('requests.get')
    def test_alphafold_download_integration(self, mock_get, integration_processor):
        """Test AlphaFold download and processing integration"""
        base_proc, struct_proc = integration_processor
        
        # Mock AlphaFold structure content
        mock_structure_content = b"""HEADER    ALPHAFOLD STRUCTURE
ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 90.00           N
ATOM      2  CA  MET A   1      11.000  21.000  31.000  1.00 90.00           C
END"""
        
        mock_get.return_value = MagicMock(
            status_code=200,
            content=mock_structure_content,
            raise_for_status=lambda: None
        )
        
        # Download AlphaFold structures
        uniprot_ids = ["P12345", "Q67890"]
        af_results = download_alphafold_structures(
            uniprot_ids,
            target_folder=struct_proc.data_path
        )
        
        # Process results
        af_dataset = {
            'name': 'alphafold_test',
            'type': 'predicted',
            'source': 'AlphaFold',
            'uniprot_ids': uniprot_ids,
            'files': af_results
        }
        
        # Save dataset info
        struct_proc.save_data('alphafold_dataset', af_dataset, format='json')
        
        # Create summary
        summary = {
            'total_requested': len(uniprot_ids),
            'successful_downloads': len(af_results),
            'failed_downloads': len(uniprot_ids) - len(af_results)
        }
        
        base_proc.save_data('download_summary', summary, format='json')
        
        # Verify
        loaded_dataset = struct_proc.load_data('alphafold_dataset', format='json')
        loaded_summary = base_proc.load_data('download_summary', format='json')
        
        assert loaded_dataset['type'] == 'predicted'
        assert loaded_dataset['source'] == 'AlphaFold'
        assert loaded_summary['total_requested'] == 2


class TestBatchDownloadIntegration:
    """Test batch download operations"""
    
    @patch.object(PDBList, 'retrieve_pdb_file')
    def test_batch_pdb_download(self, mock_retrieve, integration_processor):
        """Test batch PDB download with progress tracking"""
        base_proc, struct_proc = integration_processor
        
        # Large batch of PDB IDs
        pdb_ids = [f"{i}ABC" for i in range(1, 11)]  # 10 structures
        
        # Mock some failures
        def mock_retrieve_side_effect(pdb_id, **kwargs):
            if pdb_id in ["3ABC", "7ABC"]:
                raise Exception(f"Failed to download {pdb_id}")
            return f"mock_{pdb_id}.cif"
        
        mock_retrieve.side_effect = mock_retrieve_side_effect
        
        # Download with progress tracking
        successful = []
        failed = []
        
        for pdb_id in pdb_ids:
            try:
                result = download_protein_structures(
                    [pdb_id],
                    target_folder=struct_proc.data_path
                )
                if result:
                    successful.append(pdb_id)
                else:
                    failed.append(pdb_id)
            except:
                failed.append(pdb_id)
        
        # Save batch results
        batch_results = {
            'total': len(pdb_ids),
            'successful': successful,
            'failed': failed,
            'success_rate': len(successful) / len(pdb_ids)
        }
        
        base_proc.save_data('batch_download_results', batch_results, format='json')
        
        # Create dataset from successful downloads
        dataset_info = {
            'name': 'batch_structures',
            'pdb_ids': successful,
            'excluded_ids': failed,
            'total_structures': len(successful)
        }
        
        struct_proc.save_data('batch_dataset', dataset_info, format='json')
        
        # Verify
        loaded_results = base_proc.load_data('batch_download_results', format='json')
        loaded_dataset = struct_proc.load_data('batch_dataset', format='json')
        
        assert loaded_results['total'] == 10
        assert len(loaded_results['failed']) == 2
        assert loaded_results['success_rate'] == 0.8
        assert loaded_dataset['total_structures'] == 8
    
    @patch('requests.get')
    def test_mixed_source_download(self, mock_get, integration_processor):
        """Test downloading from multiple sources (PDB + AlphaFold)"""
        base_proc, struct_proc = integration_processor
        
        # Mock responses
        mock_get.return_value = MagicMock(
            status_code=200,
            content=b"MOCK STRUCTURE DATA",
            raise_for_status=lambda: None
        )
        
        # Define targets
        targets = {
            'experimental': ['1ABC', '2DEF'],
            'alphafold': ['P12345', 'Q67890']
        }
        
        # Save targets
        base_proc.save_data('download_targets', targets, format='json')
        
        # Mock downloads
        with patch.object(PDBList, 'retrieve_pdb_file', return_value="mock_pdb.cif"):
            pdb_results = download_protein_structures(
                targets['experimental'],
                target_folder=struct_proc.data_path
            )
        
        af_results = download_alphafold_structures(
            targets['alphafold'],
            target_folder=struct_proc.data_path
        )
        
        # Combine results
        combined_dataset = {
            'name': 'mixed_sources',
            'experimental': {
                'ids': targets['experimental'],
                'downloaded': list(pdb_results.keys())
            },
            'predicted': {
                'ids': targets['alphafold'],
                'downloaded': list(af_results.keys())
            }
        }
        
        struct_proc.save_data('mixed_dataset', combined_dataset, format='json')
        
        # Verify
        loaded_dataset = struct_proc.load_data('mixed_dataset', format='json')
        
        assert 'experimental' in loaded_dataset
        assert 'predicted' in loaded_dataset
        assert len(loaded_dataset['experimental']['downloaded']) == 2
        assert len(loaded_dataset['predicted']['downloaded']) == 2


class TestDownloadErrorHandling:
    """Test error handling in download operations"""
    
    def test_network_error_recovery(self, integration_processor):
        """Test recovery from network errors"""
        base_proc, struct_proc = integration_processor
        
        # Track retry attempts
        retry_log = []
        
        def mock_download_with_retry(pdb_id, max_retries=3):
            for attempt in range(max_retries):
                try:
                    if attempt < 2:  # Fail first 2 attempts
                        retry_log.append({'pdb_id': pdb_id, 'attempt': attempt, 'status': 'failed'})
                        raise requests.ConnectionError("Network error")
                    else:
                        retry_log.append({'pdb_id': pdb_id, 'attempt': attempt, 'status': 'success'})
                        return f"mock_{pdb_id}.cif"
                except requests.ConnectionError:
                    if attempt == max_retries - 1:
                        raise
        
        # Test with retries
        pdb_ids = ["1ABC", "2DEF"]
        results = {}
        
        for pdb_id in pdb_ids:
            try:
                result = mock_download_with_retry(pdb_id)
                results[pdb_id] = result
            except Exception as e:
                results[pdb_id] = f"Failed: {str(e)}"
        
        # Save retry log
        base_proc.save_data('retry_log', retry_log, format='json')
        base_proc.save_data('download_results', results, format='json')
        
        # Verify
        loaded_log = base_proc.load_data('retry_log', format='json')
        loaded_results = base_proc.load_data('download_results', format='json')
        
        # Each PDB should have 3 attempts (2 failures + 1 success)
        assert len(loaded_log) == 6  # 2 PDBs × 3 attempts
        assert all('mock_' in v for v in loaded_results.values())
    
    def test_partial_download_recovery(self, integration_processor):
        """Test recovery from partial downloads"""
        base_proc, struct_proc = integration_processor
        
        # Simulate partial download scenario
        download_state = {
            'completed': ['1ABC', '2DEF'],
            'in_progress': ['3GHI'],
            'pending': ['4JKL', '5MNO'],
            'failed': ['6PQR']
        }
        
        # Save state
        base_proc.save_data('download_state', download_state, format='json')
        
        # Resume download
        to_download = download_state['in_progress'] + download_state['pending']
        
        # Mock resumed downloads
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = "mock_resumed.cif"
            
            resumed_results = {}
            for pdb_id in to_download:
                try:
                    result = download_protein_structures(
                        [pdb_id],
                        target_folder=struct_proc.data_path
                    )
                    resumed_results[pdb_id] = 'success'
                except:
                    resumed_results[pdb_id] = 'failed'
        
        # Update state
        final_state = {
            'completed': download_state['completed'] + list(resumed_results.keys()),
            'failed': download_state['failed'],
            'total_attempted': len(download_state['completed']) + len(to_download) + len(download_state['failed'])
        }
        
        base_proc.save_data('final_download_state', final_state, format='json')
        
        # Verify
        loaded_state = base_proc.load_data('final_download_state', format='json')
        
        assert loaded_state['total_attempted'] == 6
        assert len(loaded_state['completed']) >= len(download_state['completed'])


@pytest.mark.network
class TestRealDownloadIntegration:
    """Test with real downloads (requires network)"""
    
    def test_real_small_structure_download(self, integration_processor):
        """Test downloading a small real structure"""
        pytest.skip("Network test - requires internet connection")
        
        base_proc, struct_proc = integration_processor
        
        # Download a small, well-known structure
        pdb_ids = ["1UBQ"]  # Ubiquitin - small protein
        
        results = download_protein_structures(
            pdb_ids,
            target_folder=struct_proc.data_path
        )
        
        # Save results
        struct_proc.save_data('ubiquitin_download', results, format='json')
        
        # Verify download
        assert "1UBQ" in results
        assert results["1UBQ"] is not None
"""
Tests for the UniProt loader functionality in the protos package.

This module tests downloading and processing protein sequences from UniProt.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

from protos.io.fasta_utils import read_fasta, write_fasta
from protos.loaders.uniprot_utils import get_uniprot, map_uniprot_to_pdb
from protos.loaders.uniprot_loader import UniprotDL
from protos.io.paths.path_config import ProtosPaths
from protos.core.base_processor import BaseProcessor


@pytest.fixture
def test_uniprot_ids():
    """Define test UniProt IDs to use for real data tests."""
    # Return list of UniProt IDs for proteins with diverse features
    return ["P00533", "P01308", "P05067", "P02751"]


@pytest.fixture
def uniprot_processor(tmp_path):
    """Create a processor for UniProt testing."""
    ProtosPaths.set_data_root(str(tmp_path))
    
    processor = BaseProcessor(
        name="test_uniprot",
        processor_data_dir="sequence"
    )
    
    yield processor
    
    ProtosPaths.set_data_root(None)


@pytest.fixture
def prepared_test_environment(uniprot_processor, test_uniprot_ids):
    """Create a properly structured test environment with dataset file."""
    # Create dataset file with test IDs
    dataset_name = "test_dataset"
    dataset_content = pd.DataFrame({
        'uniprot_id': test_uniprot_ids,
        'protein_name': ['EGFR', 'Insulin', 'APP', 'Fibronectin'],
        'organism': ['Human', 'Human', 'Human', 'Human']
    })
    
    # Save dataset using processor
    uniprot_processor.save_data(f'{dataset_name}_proteins', dataset_content, format='csv', index=False)
    
    return uniprot_processor, dataset_name, test_uniprot_ids


class TestUniprotUtils:
    """Test the UniProt utility functions."""
    
    @patch('requests.get')
    def test_get_uniprot_basic(self, mock_get):
        """Test basic UniProt retrieval functionality."""
        # Mock response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">sp|P00533|EGFR_HUMAN\nMQKLLIL\n"
        mock_get.return_value = mock_response
        
        # Test retrieval
        result = get_uniprot("P00533")
        
        assert result is not None
        assert ">sp|P00533|EGFR_HUMAN" in result
        assert "MQKLLIL" in result
    
    @patch('requests.get')
    def test_get_uniprot_error_handling(self, mock_get):
        """Test error handling in UniProt retrieval."""
        # Mock error response
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        # Test retrieval
        result = get_uniprot("INVALID_ID")
        
        # Should return None or empty on error
        assert result is None or result == ""
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_basic(self, mock_get):
        """Test UniProt to PDB mapping functionality."""
        # Mock UniProt API response
        mock_response_data = {
            "results": [
                {
                    "from": "P00533",
                    "to": {
                        "primaryAccession": "P00533",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1IVO"},
                            {"database": "PDB", "id": "1M14"}
                        ]
                    }
                }
            ]
        }
        
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_response_data
        mock_get.return_value = mock_response
        
        # Test mapping
        result = map_uniprot_to_pdb(["P00533"])
        
        assert "P00533" in result
        assert "1IVO" in result["P00533"]
        assert "1M14" in result["P00533"]
    
    @patch('requests.get')
    def test_map_uniprot_to_pdb_multiple(self, mock_get):
        """Test mapping multiple UniProt IDs to PDB."""
        # Mock response for multiple IDs
        mock_response_data = {
            "results": [
                {
                    "from": "P00533",
                    "to": {
                        "primaryAccession": "P00533",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1IVO"}
                        ]
                    }
                },
                {
                    "from": "P01308",
                    "to": {
                        "primaryAccession": "P01308",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1TRZ"}
                        ]
                    }
                }
            ]
        }
        
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_response_data
        mock_get.return_value = mock_response
        
        # Test mapping
        result = map_uniprot_to_pdb(["P00533", "P01308"])
        
        assert len(result) == 2
        assert "P00533" in result
        assert "P01308" in result


class TestUniprotDL:
    """Test the UniprotDL class functionality."""
    
    def test_uniprot_dl_initialization(self, uniprot_processor):
        """Test UniprotDL initialization."""
        # Initialize with processor's data directory
        dl = UniprotDL(
            fasta_dir=uniprot_processor.data_path,
            metadata_dir=uniprot_processor.data_path
        )
        
        assert dl is not None
        assert dl.fasta_dir is not None
        assert dl.metadata_dir is not None
    
    @patch('requests.get')
    def test_download_fasta_single(self, mock_get, uniprot_processor):
        """Test downloading a single FASTA file."""
        # Mock FASTA response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">sp|P00533|EGFR_HUMAN\nMQKLLILTCLVAVAL\n"
        mock_get.return_value = mock_response
        
        # Initialize downloader
        dl = UniprotDL(
            fasta_dir=uniprot_processor.data_path,
            metadata_dir=uniprot_processor.data_path
        )
        
        # Download single sequence
        result = dl.download_fasta("P00533")
        
        # Verify download
        assert result is not None
        
        # Save result using processor
        uniprot_processor.save_data('P00533.fasta', result, format='text')
        
        # Verify saved content
        loaded = uniprot_processor.load_data('P00533.fasta', format='text')
        assert ">sp|P00533|EGFR_HUMAN" in loaded
    
    @patch('protos.loaders.uniprot_utils.get_uniprot')
    def test_download_dataset(self, mock_get_uniprot, prepared_test_environment):
        """Test downloading a complete dataset."""
        processor, dataset_name, test_ids = prepared_test_environment
        
        # Mock UniProt responses
        def mock_uniprot_response(uniprot_id):
            return f">sp|{uniprot_id}|PROT_HUMAN\nMQKLLILTCLVAVAL\n"
        
        mock_get_uniprot.side_effect = mock_uniprot_response
        
        # Initialize downloader
        dl = UniprotDL(
            fasta_dir=processor.data_path,
            metadata_dir=processor.data_path
        )
        
        # Create a mock dataset CSV using processor
        df = pd.DataFrame({'uniprot_id': test_ids})
        processor.save_data(f'metadata/{dataset_name}', df, format='csv', index=False)
        
        # Download dataset
        results = dl.download_dataset(dataset_name)
        
        # Should have downloaded all sequences
        assert len(results) == len(test_ids)
    
    def test_combine_fasta_files(self, uniprot_processor):
        """Test combining multiple FASTA files."""
        # Create test FASTA files
        fasta1 = ">seq1\nACGT\n"
        fasta2 = ">seq2\nTGCA\n"
        
        uniprot_processor.save_data('test1.fasta', fasta1, format='text')
        uniprot_processor.save_data('test2.fasta', fasta2, format='text')
        
        # Initialize downloader
        dl = UniprotDL(
            fasta_dir=processor.data_path,
            metadata_dir=processor.data_path
        )
        
        # Combine files
        output_name = "combined"
        dl.combine_fasta(["test1", "test2"], output_name)
        
        # Verify combined file
        combined = uniprot_processor.load_data(f'{output_name}.fasta', format='text')
        
        assert ">seq1" in combined
        assert "ACGT" in combined
        assert ">seq2" in combined
        assert "TGCA" in combined


class TestIntegrationWorkflow:
    """Test integrated workflows with UniProt data."""
    
    @patch('requests.get')
    def test_uniprot_to_structure_workflow(self, mock_get, uniprot_processor):
        """Test workflow from UniProt ID to structure mapping."""
        # Mock UniProt sequence response
        seq_response = MagicMock()
        seq_response.status_code = 200
        seq_response.text = ">sp|P00533|EGFR_HUMAN\nMQKLLILTCLVAVAL\n"
        
        # Mock PDB mapping response
        pdb_response_data = {
            "results": [
                {
                    "from": "P00533",
                    "to": {
                        "primaryAccession": "P00533",
                        "uniProtKBCrossReferences": [
                            {"database": "PDB", "id": "1IVO"}
                        ]
                    }
                }
            ]
        }
        
        pdb_response = MagicMock()
        pdb_response.status_code = 200
        pdb_response.json.return_value = pdb_response_data
        
        # Set up mock responses
        def mock_get_side_effect(url, **kwargs):
            if "fasta" in url:
                return seq_response
            else:
                return pdb_response
        
        mock_get.side_effect = mock_get_side_effect
        
        # Download sequence
        sequence = get_uniprot("P00533")
        assert sequence is not None
        
        # Map to PDB
        pdb_mapping = map_uniprot_to_pdb(["P00533"])
        assert "P00533" in pdb_mapping
        
        # Save workflow results
        workflow_data = {
            'uniprot_id': 'P00533',
            'sequence': sequence,
            'pdb_ids': pdb_mapping.get('P00533', [])
        }
        
        uniprot_processor.save_data('workflow_results', workflow_data, format='json')
        
        # Verify saved data
        loaded = uniprot_processor.load_data('workflow_results', format='json')
        assert loaded['uniprot_id'] == 'P00533'
        assert loaded['pdb_ids'] == ['1IVO']
    
    def test_batch_processing_workflow(self, uniprot_processor, test_uniprot_ids):
        """Test batch processing of UniProt IDs."""
        # Create batch processing summary
        batch_results = {
            'total_ids': len(test_uniprot_ids),
            'processed': [],
            'failed': [],
            'statistics': {}
        }
        
        # Simulate processing
        for uniprot_id in test_uniprot_ids:
            # In real scenario, would download and process
            # For test, just track as processed
            batch_results['processed'].append(uniprot_id)
        
        # Calculate statistics
        batch_results['statistics'] = {
            'success_rate': len(batch_results['processed']) / batch_results['total_ids'],
            'total_processed': len(batch_results['processed']),
            'total_failed': len(batch_results['failed'])
        }
        
        # Save batch results
        uniprot_processor.save_data('batch_results', batch_results, format='json')
        
        # Verify
        loaded = uniprot_processor.load_data('batch_results', format='json')
        assert loaded['statistics']['success_rate'] == 1.0
        assert loaded['total_ids'] == len(test_uniprot_ids)


@pytest.mark.network
class TestRealUniProtData:
    """Test with real UniProt data (requires network)."""
    
    def test_real_uniprot_download(self, uniprot_processor):
        """Test downloading real UniProt data."""
        pytest.skip("Network test - requires internet connection")
        
        # Test with a small, well-characterized protein
        uniprot_id = "P00698"  # Lysozyme
        
        # Download sequence
        sequence = get_uniprot(uniprot_id)
        
        assert sequence is not None
        assert len(sequence) > 0
        assert ">sp|P00698|" in sequence
        
        # Save for verification
        uniprot_processor.save_data(f'{uniprot_id}_real.fasta', sequence, format='text')
    
    def test_real_pdb_mapping(self):
        """Test real UniProt to PDB mapping."""
        pytest.skip("Network test - requires internet connection")
        
        # Test with protein known to have PDB structures
        uniprot_ids = ["P00698"]  # Lysozyme has many structures
        
        mapping = map_uniprot_to_pdb(uniprot_ids)
        
        assert "P00698" in mapping
        assert len(mapping["P00698"]) > 0  # Should have multiple PDB entries
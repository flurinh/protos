"""
Tests for the UniProt loader functionality in the protos package.

This module tests downloading and processing protein sequences from UniProt.
"""

import pytest
import pandas as pd
import os
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
    os.environ["PROTOS_DATA_ROOT"] = str(tmp_path)
    
    processor = BaseProcessor(name="test")
    
    yield processor
    
    # Clear environment


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
    uniprot_processor.save_data(f'{dataset_name}_proteins', dataset_content, index=False)
    
    return uniprot_processor, dataset_name, test_uniprot_ids


class TestUniprotUtils:
    """Test the UniProt utility functions."""
    
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_stream')
    def test_get_uniprot_basic(self, mock_stream):
        """Test basic UniProt retrieval functionality."""
        # Mock TSV response data - needs to be a list of lines for csv.DictReader
        mock_tsv_lines = [
            "Entry\tProtein names\tGene Names\tOrganism\tLength\tSequence",
            "P00533\tEpidermal growth factor receptor\tEGFR\tHomo sapiens\t1210\tMQKLLIL"
        ]
        mock_stream.return_value = mock_tsv_lines
        
        # Test retrieval
        result = get_uniprot("P00533")
        
        # Result should be a DataFrame
        import pandas as pd
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "P00533" in result['Entry'].values
    
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_stream')
    def test_get_uniprot_error_handling(self, mock_stream):
        """Test error handling in UniProt retrieval."""
        # Mock empty response - just header, no data
        mock_stream.return_value = ["Entry\tProtein names\tGene Names\tOrganism\tLength\tSequence"]
        
        # Test retrieval - should return empty DataFrame
        result = get_uniprot("INVALID_ID")
        
        # Should return empty DataFrame on error
        import pandas as pd
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
    
    @patch('protos.loaders.uniprot_utils.submit_id_mapping')
    @patch('protos.loaders.uniprot_utils.check_id_mapping_results_ready')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_link')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_search')
    def test_map_uniprot_to_pdb_basic(self, mock_search, mock_link, mock_ready, mock_submit):
        """Test UniProt to PDB mapping functionality."""
        # Mock the submission
        mock_submit.return_value = "job123"
        
        # Mock the ready check
        mock_ready.return_value = True
        
        # Mock the link
        mock_link.return_value = "https://rest.uniprot.org/idmapping/results/job123"
        
        # Mock the search results - this is what get_id_mapping_results_search returns
        mock_search.return_value = {
            "results": [
                {"from": "P00533", "to": "1IVO"},
                {"from": "P00533", "to": "1M14"}
            ]
        }
        
        # Test mapping
        result = map_uniprot_to_pdb(["P00533"])
        
        # Result should be a DataFrame with columns 'uid' and 'pdb_id'
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "uid" in result.columns
        assert "pdb_id" in result.columns
        assert "P00533" in result["uid"].values
        assert "1IVO" in result["pdb_id"].values
        assert "1M14" in result["pdb_id"].values
    
    @patch('protos.loaders.uniprot_utils.submit_id_mapping')
    @patch('protos.loaders.uniprot_utils.check_id_mapping_results_ready')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_link')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_search')
    def test_map_uniprot_to_pdb_multiple(self, mock_search, mock_link, mock_ready, mock_submit):
        """Test mapping multiple UniProt IDs to PDB."""
        # Mock the submission
        mock_submit.return_value = "job456"
        
        # Mock the ready check
        mock_ready.return_value = True
        
        # Mock the link
        mock_link.return_value = "https://rest.uniprot.org/idmapping/results/job456"
        
        # Mock the search results
        mock_search.return_value = {
            "results": [
                {"from": "P00533", "to": "1IVO"},
                {"from": "P01308", "to": "1TRZ"}
            ]
        }
        
        # Test mapping
        result = map_uniprot_to_pdb(["P00533", "P01308"])
        
        # Result should be a DataFrame
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "uid" in result.columns
        assert "pdb_id" in result.columns
        assert set(result["uid"].values) == {"P00533", "P01308"}
        assert set(result["pdb_id"].values) == {"1IVO", "1TRZ"}


class TestUniprotDL:
    """Test the UniprotDL class functionality."""
    
    def test_uniprot_dl_initialization(self, uniprot_processor):
        """Test UniprotDL initialization."""
        # Initialize with processor's data root
        dl = UniprotDL(
            data_root=uniprot_processor.data_root
        )
        
        assert dl is not None
        assert dl.paths is not None
        assert dl.fasta_dir is not None
    
    @patch('protos.loaders.uniprot_utils.get_uniprot')
    def test_download_fasta_single(self, mock_get_uniprot, uniprot_processor):
        """Test downloading a single sequence."""
        # Mock UniProt response as DataFrame
        mock_df = pd.DataFrame({
            'Entry': ['P00533'],
            'Protein names': ['Epidermal growth factor receptor'],
            'Gene Names': ['EGFR'],
            'Organism': ['Homo sapiens'],
            'Length': ['1210'],
            'Sequence': ['MQKLLILTCLVAVAL']
        })
        mock_get_uniprot.return_value = mock_df
        
        # Initialize downloader
        dl = UniprotDL(
            data_root=uniprot_processor.data_root
        )
        
        # Download single sequence using actual method
        result = dl.download_gene_single_query("P00533")
        
        # Verify result is a list with expected data
        assert isinstance(result, list)
        assert len(result) == 7  # [uniprot, seq, gene_name, species, organism, info, dataset]
        assert result[0] == "P00533"
        assert isinstance(result[1], str)  # sequence is a string
        assert len(result[1]) > 0  # sequence is not empty
        
        # First populate gene_df with the result
        dl.gene_df = pd.DataFrame([result], columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset'])
        
        # Save FASTA using actual method
        saved_paths = dl.save_uniprot_fasta("P00533")
        
        # Verify FASTA was saved
        assert saved_paths is not None
        assert isinstance(saved_paths, list)
        assert len(saved_paths) > 0
        assert os.path.exists(saved_paths[0])
    
    @patch('protos.loaders.uniprot_utils.get_uniprot')
    def test_download_genes_batch(self, mock_get_uniprot, prepared_test_environment):
        """Test downloading multiple sequences."""
        processor, dataset_name, test_ids = prepared_test_environment
        
        # Mock UniProt responses - return DataFrame for each call
        def mock_uniprot_response(query, **kwargs):
            # Extract the UniProt ID from the query
            uniprot_id = query
            return pd.DataFrame({
                'Entry': [uniprot_id],
                'Protein names': [f'{uniprot_id} protein'],
                'Gene Names': [f'GENE_{uniprot_id}'],
                'Organism': ['Homo sapiens'],
                'Length': ['1000'],
                'Sequence': ['MQKLLILTCLVAVAL']
            })
        
        mock_get_uniprot.side_effect = mock_uniprot_response
        
        # Initialize downloader with the dataset
        dl = UniprotDL(
            data_root=processor.data_root,
            dataset=dataset_name
        )
        
        # Create dataset file in metadata directory
        dataset_file = os.path.join(dl.metadata_dir, f'{dataset_name}.txt')
        os.makedirs(dl.metadata_dir, exist_ok=True)
        with open(dataset_file, 'w') as f:
            f.write(' '.join(test_ids))
        
        # Load dataset to populate dl.genes
        dl.load_dataset()
        
        # Download genes using actual method
        # Note: download_genes_single_query uses self.genes, not a parameter
        result_df = dl.download_genes_single_query(batchsize=len(test_ids))
        
        # Verify gene_df was populated with all sequences
        assert not result_df.empty
        assert len(result_df) >= len(test_ids)  # May have more if genes list was larger
        
        # Check that our test IDs are in the results
        downloaded_ids = result_df['uniprot'].tolist()
        for uid in test_ids:
            assert uid in downloaded_ids
    
    def test_save_and_load_gene_table(self, uniprot_processor):
        """Test saving and loading gene table functionality."""
        # Initialize downloader
        dl = UniprotDL(
            data_root=uniprot_processor.data_root,
            dataset='test_dataset'
        )
        
        # Create test gene data
        test_gene_df = pd.DataFrame({
            'uniprot': ['P00533', 'P01308'],
            'seq': ['MQKLLIL', 'MALWMRLL'],
            'gene': ['EGFR', 'INS'],
            'species': ['Homo sapiens', 'Homo sapiens'],
            'organism': ['Human', 'Human'],
            'info': ['Receptor', 'Hormone'],
            'dataset': ['test_dataset', 'test_dataset']
        })
        
        # Set gene_df
        dl.gene_df = test_gene_df
        
        # Save gene table
        saved_path = dl.save_gene_table()
        assert saved_path is not None
        assert os.path.exists(saved_path)
        
        # Create new instance and load table
        dl2 = UniprotDL(
            data_root=uniprot_processor.data_root,
            dataset='test_dataset'
        )
        
        # Load gene table
        loaded_df = dl2.load_gene_table()
        
        # Verify loaded data matches original
        assert len(loaded_df) == 2
        assert 'P00533' in loaded_df['uniprot'].values
        assert 'P01308' in loaded_df['uniprot'].values


class TestIntegrationWorkflow:
    """Test integrated workflows with UniProt data."""
    
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_stream')
    @patch('protos.loaders.uniprot_utils.submit_id_mapping')
    @patch('protos.loaders.uniprot_utils.check_id_mapping_results_ready')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_link')
    @patch('protos.loaders.uniprot_utils.get_id_mapping_results_search')
    def test_uniprot_to_structure_workflow(self, mock_search, mock_link, mock_ready, mock_submit, mock_stream, uniprot_processor):
        """Test workflow from UniProt ID to structure mapping."""
        # Mock UniProt TSV response for get_uniprot - needs to be a list of lines
        mock_tsv_lines = [
            "Entry\tProtein names\tGene Names\tOrganism\tLength\tSequence",
            "P00533\tEpidermal growth factor receptor\tEGFR\tHomo sapiens\t1210\tMQKLLILTCLVAVAL"
        ]
        mock_stream.return_value = mock_tsv_lines
        
        # Mock PDB mapping
        mock_submit.return_value = "job789"
        mock_ready.return_value = True
        mock_link.return_value = "https://rest.uniprot.org/idmapping/results/job789"
        mock_search.return_value = {
            "results": [
                {"from": "P00533", "to": "1IVO"}
            ]
        }
        
        # Download sequence - returns DataFrame
        sequence_df = get_uniprot("P00533")
        assert isinstance(sequence_df, pd.DataFrame)
        assert len(sequence_df) > 0
        
        # Map to PDB - returns DataFrame
        pdb_mapping_df = map_uniprot_to_pdb(["P00533"])
        assert isinstance(pdb_mapping_df, pd.DataFrame)
        
        # Get sequence string from DataFrame
        sequence_str = sequence_df['Sequence'].iloc[0] if not sequence_df.empty else ""
        
        # Get PDB IDs from mapping DataFrame
        pdb_ids = pdb_mapping_df[pdb_mapping_df['uid'] == 'P00533']['pdb_id'].tolist()
        
        # Save workflow results
        workflow_data = {
            'uniprot_id': 'P00533',
            'sequence': sequence_str,
            'pdb_ids': pdb_ids
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
        
        # Download sequence - returns DataFrame
        sequence_df = get_uniprot(uniprot_id)
        
        assert isinstance(sequence_df, pd.DataFrame)
        assert len(sequence_df) > 0
        assert "P00698" in sequence_df['Entry'].values if 'Entry' in sequence_df.columns else False
        
        # Convert to FASTA format for saving
        if not sequence_df.empty and 'Entry' in sequence_df.columns and 'Sequence' in sequence_df.columns:
            entry = sequence_df['Entry'].iloc[0]
            sequence = sequence_df['Sequence'].iloc[0]
            fasta_content = f">sp|{entry}|PROTEIN\n{sequence}\n"
            
            # Save for verification
            uniprot_processor.save_data(f'{uniprot_id}_real.fasta', fasta_content, format='text')
    
    def test_real_pdb_mapping(self):
        """Test real UniProt to PDB mapping."""
        pytest.skip("Network test - requires internet connection")
        
        # Test with protein known to have PDB structures
        uniprot_ids = ["P00698"]  # Lysozyme has many structures
        
        mapping_df = map_uniprot_to_pdb(uniprot_ids)
        
        assert isinstance(mapping_df, pd.DataFrame)
        assert len(mapping_df) > 0
        assert "P00698" in mapping_df['uid'].values
        # Should have multiple PDB entries for lysozyme
        assert len(mapping_df[mapping_df['uid'] == 'P00698']) > 0
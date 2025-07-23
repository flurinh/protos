"""
Tests for StructureProcessor entity handling functionality.
"""

import os
import pytest
import pandas as pd
from pathlib import Path

from protos.processing.structure import StructureProcessor
from protos.io.paths.path_config import ProtosPaths
from protos.io.data_access import generate_entity_id


class TestCifBaseProcessorEntityHandling:
    """Test entity management in StructureProcessor."""
    
    @pytest.fixture
    def processor(self, configure_test_paths):
        """Create a test StructureProcessor with entity support."""
        processor = StructureProcessor(
            name="test_cif_entity"
        )
        yield processor
    
    def test_register_structure_entity(self, processor):
        """Test that loading a structure registers it as an entity."""
        # Use a real small protein structure
        pdb_id = '1ubq'
        
        # Download the structure if not already present
        from protos.loaders.download_structures import download_structures_with_processor
        successful, failed = download_structures_with_processor([pdb_id], processor)
        
        if pdb_id not in successful:
            pytest.skip(f"Could not download {pdb_id} for testing")
        
        # Load the structure
        structure_df = processor.load_structure(pdb_id)
        
        # Register the structure
        entity_id = processor._register_structure_entity(pdb_id, structure_df)
        
        # Verify entity ID is the PDB ID (human-readable name)
        assert entity_id == pdb_id
        
        # Check that entity was registered
        assert processor.entity_registry.entity_exists(pdb_id)
    
    def test_get_entity_id_for_pdb(self, processor):
        """Test getting entity ID for a PDB ID."""
        # Get entity ID for a PDB
        entity_id1 = processor.get_entity_id_for_pdb('2gb1')
        assert len(entity_id1) == 10
        assert entity_id1.isalnum()
        
        # Same PDB should give same entity ID
        entity_id2 = processor.get_entity_id_for_pdb('2gb1')
        assert entity_id1 == entity_id2
        
        # Different PDB should give different entity ID
        entity_id3 = processor.get_entity_id_for_pdb('1crn')
        assert entity_id1 != entity_id3
    
    def test_list_structure_entities(self, processor):
        """Test listing structure entities."""
        # Use real small protein structures
        pdb_ids = ['1ubq', '2gb1', '1crn']
        
        # Download and register real structures
        from protos.loaders.download_structures import download_structures_with_processor
        successful, failed = download_structures_with_processor(pdb_ids, processor)
        
        # Ensure we downloaded at least some structures
        assert len(successful) > 0, f"Failed to download any structures. Failed: {failed}"
        
        # List all entities - should return PDB IDs
        entities = processor.list_entities()
        
        # Should have at least our downloaded ones
        assert len(entities) >= len(successful)
        for pdb_id in successful:
            assert pdb_id in entities
    
    def test_create_dataset_with_entity_ids(self, processor):
        """Test creating a dataset that uses entity IDs."""
        # Set up some PDB IDs
        pdb_ids = ['1ubq', '2gb1', '1crn']
        processor.pdb_ids = pdb_ids
        
        # Create dataset with entity IDs
        dataset_id = processor.create_dataset(
            dataset_id="test_entity_dataset",
            pdb_ids=pdb_ids,
            metadata={
                "name": "Test Entity Dataset",
                "description": "Dataset using entity IDs",
                "use_entity_ids": True
            }
        )
        
        # Check the dataset was created
        assert dataset_id == "test_entity_dataset"
        
        # Load the dataset info to check its content
        dataset_info = processor.get_dataset_info("test_entity_dataset")
        if dataset_info:
            # Check metadata
            assert dataset_info.get("metadata", {}).get("use_entity_ids") is True
            
            # The dataset should store PDB IDs, not entity hashes
            # (entity conversion happens at runtime when needed)
            assert "pdb_ids" in dataset_info or "entities" in dataset_info
    
    def test_load_dataset_with_entity_ids(self, processor):
        """Test loading a dataset that contains entity IDs."""
        # Create a dataset with entity IDs (no prefix for universal IDs)
        entity_ids = [
            generate_entity_id('1ubq'),
            generate_entity_id('2gb1'),
            generate_entity_id('1crn')
        ]
        
        # Create test dataset
        if processor.dataset_manager is not None:
            dataset_id = processor.dataset_manager.create_dataset(
                name="entity_test",  # This becomes the dataset_id
                entities=entity_ids,
                metadata={
                    "description": "Test dataset with entity IDs",
                    "entity_based": True
                }
            )
            
            # Register entities with the processor's entity registry
            for i, entity_id in enumerate(entity_ids):
                pdb_id = ['1ubq', '2gb1', '1crn'][i]
                processor.entity_registry.register_entity(
                    pdb_id,  # Human-readable name (positional)
                    "structure",  # Format type (positional)
                    f"structure/mmcif/{pdb_id}.cif",  # File path (positional)
                    {"pdb_id": pdb_id, "datasets": ["entity_test"]}  # Metadata (positional)
                )
            
            # Mock load_structures to avoid actual file loading
            processor.load_structures = lambda pdb_ids, **kwargs: None
            
            # Load the dataset - should convert entity IDs to PDB IDs
            loaded_pdb_ids = processor.load_dataset(dataset_id, apply_dtypes=False)
            
            # Should have extracted the PDB IDs
            assert len(loaded_pdb_ids) == 3
            assert '1ubq' in loaded_pdb_ids
            assert '2gb1' in loaded_pdb_ids
            assert '1crn' in loaded_pdb_ids
    
    def test_save_structure_as_entity(self, processor):
        """Test saving a structure as an entity."""
        # Use a real small protein structure
        pdb_id = '1l2y'
        
        # Download the structure if not already present
        from protos.loaders.download_structures import download_structures_with_processor
        successful, failed = download_structures_with_processor([pdb_id], processor)
        
        if pdb_id not in successful:
            pytest.skip(f"Could not download {pdb_id} for testing")
        
        # Load the real structure
        test_structure = processor.load_structure(pdb_id)
        
        # Save as entity (this will overwrite the existing file)
        entity_id = processor.save_structure_as_entity(
            structure_df=test_structure,
            pdb_id=pdb_id,
            datasets=['test_dataset'],
            metadata={'resolution': 2.0, 'method': 'X-RAY'}
        )
        
        # Verify entity ID
        assert len(entity_id) == 10
        assert entity_id.isalnum()
        
        # Check file was created
        cif_path = Path(processor.paths.get_subdir_path("structure", "structure_dir")) / f'{pdb_id}.cif'
        assert cif_path.exists()
        
        # Check entity mapping
        assert processor.get_entity_id_for_pdb(pdb_id) == entity_id
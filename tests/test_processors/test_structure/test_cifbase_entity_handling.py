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
        # Create a minimal test structure
        test_structure = pd.DataFrame({
            'group': ['ATOM'] * 3,
            'pdb_id': ['1ABC'] * 3,
            'auth_chain_id': ['A', 'A', 'A'],
            'auth_seq_id': [1, 1, 2],
            'auth_comp_id': ['ALA', 'ALA', 'GLY'],
            'atom_name': ['N', 'CA', 'N'],
            'x': [1.0, 2.0, 3.0],
            'y': [4.0, 5.0, 6.0],
            'z': [7.0, 8.0, 9.0]
        })
        
        # Register the structure
        entity_id = processor._register_structure_entity('1ABC', test_structure)
        
        # Verify entity ID is the PDB ID (human-readable name)
        assert entity_id == '1ABC'
        
        # Check that entity was registered
        assert processor.entity_registry.entity_exists('1ABC')
    
    def test_get_entity_id_for_pdb(self, processor):
        """Test getting entity ID for a PDB ID."""
        # Get entity ID for a PDB
        entity_id1 = processor.get_entity_id_for_pdb('2DEF')
        assert len(entity_id1) == 10
        assert entity_id1.isalnum()
        
        # Same PDB should give same entity ID
        entity_id2 = processor.get_entity_id_for_pdb('2DEF')
        assert entity_id1 == entity_id2
        
        # Different PDB should give different entity ID
        entity_id3 = processor.get_entity_id_for_pdb('3GHI')
        assert entity_id1 != entity_id3
    
    def test_list_structure_entities(self, processor):
        """Test listing structure entities."""
        # Register some entities properly
        from protos.io.data_access import GlobalRegistry, generate_entity_id
        global_registry = GlobalRegistry()
        
        pdb_ids = ['1ABC', '2DEF', '3GHI']
        for pdb_id in pdb_ids:
            entity_id = generate_entity_id(pdb_id)
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="structure",
                original_id=pdb_id,
                file_path=None,
                metadata={"pdb_id": pdb_id},
                datasets=[]
            )
        
        # List all entities - should return PDB IDs, not hashes
        entities = processor.list_entities()
        
        # Should get the PDB IDs, not entity hashes
        assert len(entities) >= 3  # At least our 3 registered ones
        for pdb_id in pdb_ids:
            assert pdb_id in entities
    
    def test_create_dataset_with_entity_ids(self, processor):
        """Test creating a dataset that uses entity IDs."""
        # Set up some PDB IDs
        pdb_ids = ['1ABC', '2DEF', '3GHI']
        processor.pdb_ids = pdb_ids
        
        # Create dataset with entity IDs
        dataset = processor.create_dataset(
            dataset_id="test_entity_dataset",
            name="Test Entity Dataset",
            description="Dataset using entity IDs",
            use_entity_ids=True
        )
        
        if dataset is not None:
            # Dataset content should be entity IDs
            assert len(dataset.content) == 3
            for entity_id in dataset.content:
                assert len(entity_id) == 10
                assert entity_id.isalnum()
            
            # Metadata should indicate entity-based
            assert dataset.metadata.get("entity_based") is True
    
    def test_load_dataset_with_entity_ids(self, processor):
        """Test loading a dataset that contains entity IDs."""
        # Create a dataset with entity IDs (no prefix for universal IDs)
        entity_ids = [
            generate_entity_id('1ABC'),
            generate_entity_id('2DEF'),
            generate_entity_id('3GHI')
        ]
        
        # Create test dataset
        if processor.dataset_manager is not None:
            dataset = processor.dataset_manager.create_dataset(
                dataset_id="entity_test",
                name="Entity Test Dataset",
                description="Test dataset with entity IDs",
                content=entity_ids,
                metadata={"entity_based": True}
            )
            
            # Register entities with original IDs using global registry
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            for i, entity_id in enumerate(entity_ids):
                pdb_id = ['1ABC', '2DEF', '3GHI'][i]
                global_registry.entity_registry.register_entity(
                    entity_id=entity_id,
                    entity_type="structure",
                    original_id=pdb_id,
                    file_path=f"/data/structures/{pdb_id}.cif",
                    datasets=["entity_test"]
                )
            
            # Mock load_structures to avoid actual file loading
            processor.load_structures = lambda pdb_ids, **kwargs: None
            
            # Load the dataset - should convert entity IDs to PDB IDs
            loaded_pdb_ids = processor.load_dataset("entity_test", apply_dtypes=False)
            
            # Should have extracted the PDB IDs
            assert len(loaded_pdb_ids) == 3
            assert '1ABC' in loaded_pdb_ids
            assert '2DEF' in loaded_pdb_ids
            assert '3GHI' in loaded_pdb_ids
    
    def test_save_structure_as_entity(self, processor):
        """Test saving a structure as an entity."""
        # Create test structure
        test_structure = pd.DataFrame({
            'group': ['ATOM'] * 3,
            'model_number': [1] * 3,
            'pdb_id': ['4JKL'] * 3,
            'auth_chain_id': ['A', 'A', 'A'],
            'entity_id': [1, 1, 1],
            'auth_seq_id': [1, 1, 2],
            'auth_comp_id': ['ALA', 'ALA', 'GLY'],
            'label_comp_id': ['ALA', 'ALA', 'GLY'],
            'atom_name': ['N', 'CA', 'N'],
            'label_atom_id': ['N', 'CA', 'N'],
            'element': ['N', 'C', 'N'],
            'x': [1.0, 2.0, 3.0],
            'y': [4.0, 5.0, 6.0],
            'z': [7.0, 8.0, 9.0],
            'b_factor': [20.0, 20.0, 20.0],
            'occupancy': [1.0, 1.0, 1.0],
            'atom_id': [1, 2, 3]
        })
        
        # Save as entity
        entity_id = processor.save_structure_as_entity(
            structure_df=test_structure,
            pdb_id='4JKL',
            datasets=['test_dataset'],
            metadata={'resolution': 2.0, 'method': 'X-RAY'}
        )
        
        # Verify entity ID
        assert len(entity_id) == 10
        assert entity_id.isalnum()
        
        # Check file was created
        cif_path = Path(processor.path_structure_dir) / '4JKL.cif'
        assert cif_path.exists()
        
        # Check entity mapping
        assert processor.get_entity_id_for_pdb('4JKL') == entity_id
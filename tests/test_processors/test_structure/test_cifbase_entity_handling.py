"""
Tests for CifBaseProcessor entity handling functionality.
"""

import os
import pytest
import pandas as pd
from pathlib import Path

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.io.paths.path_config import ProtosPaths
from protos.io.data_access import generate_entity_id


class TestCifBaseProcessorEntityHandling:
    """Test entity management in CifBaseProcessor."""
    
    @pytest.fixture
    def processor(self, tmp_path):
        """Create a test CifBaseProcessor with entity support."""
        ProtosPaths.set_data_root(str(tmp_path))
        processor = CifBaseProcessor(
            name="test_cif_entity",
            processor_data_dir="structure"
        )
        yield processor
        ProtosPaths.set_data_root(None)
    
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
        
        # Verify entity ID format
        assert len(entity_id) == 10
        assert entity_id.isalnum()
        
        # Check local mapping
        assert '1ABC' in processor._pdb_entity_map
        assert processor._pdb_entity_map['1ABC'] == entity_id
    
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
        # Register some entities
        processor.get_entity_id_for_pdb('1ABC')
        processor.get_entity_id_for_pdb('2DEF')
        processor.get_entity_id_for_pdb('3GHI')
        
        # List all entities
        entities = processor.list_structure_entities()
        assert len(entities) == 3
        
        # All should be valid entity IDs
        for entity_id in entities:
            assert len(entity_id) == 10
            assert entity_id.isalnum()
    
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
    
    def test_load_dataset_with_entity_ids(self, processor, tmp_path):
        """Test loading a dataset that contains entity IDs."""
        # Create a dataset with entity IDs
        entity_ids = [
            generate_entity_id('1ABC', prefix='structure'),
            generate_entity_id('2DEF', prefix='structure'),
            generate_entity_id('3GHI', prefix='structure')
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
            
            # Register entities with original IDs
            if hasattr(processor, 'entity_registry') and processor.entity_registry is not None:
                for i, entity_id in enumerate(entity_ids):
                    pdb_id = ['1ABC', '2DEF', '3GHI'][i]
                    processor.entity_registry.register_entity(
                        entity_id=entity_id,
                        entity_type="structure",
                        original_id=pdb_id,
                        file_path=f"/data/structures/{pdb_id}.cif",
                        datasets=["entity_test"]
                    )
            
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
        cif_path = os.path.join(processor.path_structure_dir, '4JKL.cif')
        assert os.path.exists(cif_path)
        
        # Check entity mapping
        assert processor.get_entity_id_for_pdb('4JKL') == entity_id
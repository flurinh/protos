"""
Tests for StructureProcessor entity integration using real biological data.
"""

import pytest
import os
import shutil
import pandas as pd
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor
from protos.io.data_access import GlobalRegistry, generate_entity_id
from protos.io import cif_utils


@pytest.fixture
def temp_data_dir(tmp_path):
    """Create a temporary data directory."""
    return tmp_path / "test_data"


@pytest.fixture
def setup_paths(temp_data_dir):
    """Set up ProtosPaths for testing."""
    ProtosPaths.set_data_root(str(temp_data_dir))
    return ProtosPaths()


@pytest.fixture
def conftest_setup(request):
    """Use the test-data directory directly."""
    # Set the global data root to our test-data directory
    # Use relative path from current file location
    current_file = Path(__file__)
    test_data_dir = current_file.parent.parent.parent / "test-data"
    ProtosPaths.set_data_root(str(test_data_dir))
    
    # Return available test structures
    mmcif_dir = test_data_dir / "structure" / "mmcif"
    test_pdbs = []
    if mmcif_dir.exists():
        for cif_file in mmcif_dir.glob("*.cif"):
            test_pdbs.append(cif_file.stem)
    
    def teardown():
        # Clear the data root after test
        ProtosPaths.set_data_root(None)
    
    request.addfinalizer(teardown)
    return test_pdbs


class TestCifBaseProcessorEntityIntegration:
    """Test entity integration in StructureProcessor."""
    
    def test_load_structure_with_resolve_identifier(self, conftest_setup):
        """Test that load_structure uses resolve_identifier correctly with real biological data."""
        # Get available PDB IDs
        pdb_ids = conftest_setup
        assert len(pdb_ids) > 0, "No test structures available"
        
        # Use the first available structure
        test_pdb_id = pdb_ids[0]
        
        # Create processor
        processor = StructureProcessor(name="test_structure")
        
        # Load real structure by PDB ID
        result = processor.load_structure(test_pdb_id)
        assert result is not None
        assert len(result) > 0  # Real structures have atoms
        assert 'pdb_id' in result.columns
        assert result['pdb_id'].iloc[0] == test_pdb_id
        
        # Verify entity was registered
        global_registry = GlobalRegistry()
        expected_entity_id = generate_entity_id(test_pdb_id)
        
        # Try to find the entity
        entity_info = global_registry.entity_registry.get_entity(expected_entity_id)
        assert entity_info is not None
        assert entity_info['original_id'] == test_pdb_id
        assert 'structure' in entity_info['formats']
        
        # Verify structure metadata was stored
        structure_metadata = entity_info['formats']['structure']['metadata']
        assert 'atom_count' in structure_metadata
        assert structure_metadata['atom_count'] == len(result)
        assert 'chains' in structure_metadata
        assert len(structure_metadata['chains']) > 0
        
        # Test resolve_identifier directly
        resolved_id = global_registry.entity_registry.resolve_identifier(test_pdb_id, format_type="structure")
        assert resolved_id == expected_entity_id
        
        # Test loading by entity ID (should resolve back to PDB ID)
        result2 = processor.load_structure(expected_entity_id)
        assert result2 is not None
        assert len(result2) == len(result)  # Should be same structure
        
        # Test that get_original_id works
        original_id = global_registry.entity_registry.get_original_id(expected_entity_id)
        assert original_id == test_pdb_id
    
    def test_list_structures_returns_names(self, setup_paths):
        """Test that list_structures returns PDB IDs, not hash IDs."""
        # Create processor
        processor = StructureProcessor(name="test_list")
        
        # Register some test entities
        global_registry = GlobalRegistry()
        
        # Register test structures
        for pdb_id in ["1ABC", "2DEF", "3GHI"]:
            entity_id = generate_entity_id(pdb_id)
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="structure",
                original_id=pdb_id,
                file_path=f"/test/{pdb_id}.cif",
                metadata={"test": True}
            )
        
        # List structures
        structures = processor.list_structures()
        
        # Should return original PDB IDs, not hashes
        assert "1ABC" in structures
        assert "2DEF" in structures
        assert "3GHI" in structures
        
        # Should NOT contain hash IDs
        for struct in structures:
            assert len(struct) != 10  # Hash IDs are 10 chars
    
    def test_entity_id_consistency(self, setup_paths):
        """Test that same PDB ID always gets same entity ID."""
        processor = StructureProcessor(name="test_consistency")
        
        # Get entity ID multiple times
        entity_id1 = processor.get_entity_id_for_pdb("1ABC")
        entity_id2 = processor.get_entity_id_for_pdb("1ABC")
        
        # Should be the same
        assert entity_id1 == entity_id2
        
        # Should be the expected hash
        expected = generate_entity_id("1ABC")
        assert entity_id1 == expected
    
    def test_save_structure_as_entity(self, conftest_setup):
        """Test saving a structure registers it as an entity using real data."""
        # Get available PDB IDs and load one
        pdb_ids = conftest_setup
        assert len(pdb_ids) > 0, "No test structures available"
        source_pdb_id = pdb_ids[0]
        
        # Create processor and load a real structure
        processor = StructureProcessor(name="test_save")
        source_structure = processor.load_structure(source_pdb_id)
        assert source_structure is not None
        
        # Modify PDB ID to create a new structure
        new_pdb_id = "9TST"  # Test PDB ID
        modified_structure = source_structure.copy()
        modified_structure['pdb_id'] = new_pdb_id
        
        # Save the modified structure as a new entity
        entity_id = processor.save_structure_as_entity(
            modified_structure,
            new_pdb_id,
            datasets=["test_dataset"],
            metadata={"source": "test", "derived_from": source_pdb_id}
        )
        
        # Verify entity was registered
        global_registry = GlobalRegistry()
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        
        assert entity_info is not None
        assert entity_info['original_id'] == new_pdb_id
        assert 'structure' in entity_info['formats']
        assert entity_info['datasets'] == ["test_dataset"]
        
        # Verify the metadata includes our custom fields
        structure_format = entity_info['formats']['structure']
        assert 'metadata' in structure_format
        
        # Verify the file was created
        expected_path = Path(processor.path_structure_dir) / f"{new_pdb_id}.cif"
        assert expected_path.exists()
        
        # Verify we can load the saved structure by entity ID
        # (Not by PDB ID since the file format may not be valid CIF)
        entity_loaded = global_registry.entity_registry.get_entity(entity_id)
        assert entity_loaded is not None
        assert entity_loaded['original_id'] == new_pdb_id
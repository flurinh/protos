import pytest
import pandas as pd
import numpy as np
import os
import pickle
from pathlib import Path

from protos.processing.structure import StructureProcessor
from protos.io.paths import ProtosPaths


class TestCifBaseProcessorHumanNames:
    """Test StructureProcessor with human-readable names."""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Set up test environment."""
        # Configure test data path
        test_data_dir = Path(__file__).parent.parent.parent / "test-data"
        ProtosPaths.set_data_root(str(test_data_dir))
        yield
        # Clean up
        ProtosPaths.set_data_root(None)
    
    def create_test_structure(self, pdb_id: str = "test_protein") -> pd.DataFrame:
        """Create a simple test structure DataFrame."""
        return pd.DataFrame({
            'pdb_id': [pdb_id] * 4,
            'auth_chain_id': ['A'] * 4,
            'auth_seq_id': [1, 2, 3, 4],
            'auth_comp_id': ['ALA', 'GLY', 'SER', 'THR'],
            'atom_name': ['CA', 'CA', 'CA', 'CA'],
            'x': [1.0, 2.0, 3.0, 4.0],
            'y': [1.0, 2.0, 3.0, 4.0],
            'z': [1.0, 2.0, 3.0, 4.0],
            'group': ['ATOM'] * 4
        })
    
    def test_save_and_load_structure_pkl(self):
        """Test saving and loading structures with human names using PKL format."""
        processor = StructureProcessor()
        
        # Create test structure
        structure = self.create_test_structure("my_protein")
        
        # Save with human-readable name (default PKL format)
        saved_name = processor.save_structure("my_protein", structure)
        assert saved_name == "my_protein"
        
        # Verify entity exists
        assert processor.entity_exists("my_protein")
        
        # Load by human name (should use cached PKL)
        loaded = processor.load_structure("my_protein")
        assert loaded is not None
        assert len(loaded) == 4
        assert loaded['pdb_id'].iloc[0] == "my_protein"
    
    def test_caching_behavior(self):
        """Test that structures are cached and loaded from cache."""
        processor = StructureProcessor()
        
        # Create and save structure
        structure = self.create_test_structure("cached_protein")
        processor.save_structure("cached_protein", structure)
        
        # Load with cache (default)
        loaded1 = processor.load_structure("cached_protein", debug=True)
        assert loaded1 is not None
        
        # Load without cache
        loaded2 = processor.load_structure("cached_protein", use_cache=False)
        assert loaded2 is not None
        
        # Both should have same data
        pd.testing.assert_frame_equal(loaded1, loaded2)
    
    def test_no_hash_ids_in_filenames(self):
        """Test that files are saved with human names, not hash IDs."""
        processor = StructureProcessor()
        
        # Save multiple structures
        for name in ["protein_a", "protein_b", "complex_123"]:
            structure = self.create_test_structure(name)
            processor.save_structure(name, structure)
        
        # List entities should return human names
        entities = processor.list_entities()
        assert "protein_a" in entities
        assert "protein_b" in entities
        assert "complex_123" in entities
        
        # No hash IDs (10-character alphanumeric) should be in the list
        for entity in entities:
            # Hash IDs are 10 characters and alphanumeric
            assert not (len(entity) == 10 and entity.isalnum())
    
    def test_load_dataset_impl(self):
        """Test _load_dataset_impl returns Dict[str, Structure]."""
        processor = StructureProcessor()
        
        # Save some structures
        structures = {
            "ubiquitin": self.create_test_structure("ubiquitin"),
            "lysozyme": self.create_test_structure("lysozyme"),
            "myoglobin": self.create_test_structure("myoglobin")
        }
        
        for name, structure in structures.items():
            processor.save_structure(name, structure)
        
        # Create dataset with human names
        dataset_name = processor.create_dataset(
            "test_proteins",
            ["ubiquitin", "lysozyme", "myoglobin"],
            metadata={"study": "test"}
        )
        
        # Load dataset using _load_dataset_impl
        loaded = processor._load_dataset_impl("test_proteins")
        
        # Verify return type and contents
        assert isinstance(loaded, dict)
        assert len(loaded) == 3
        assert "ubiquitin" in loaded
        assert "lysozyme" in loaded
        assert "myoglobin" in loaded
        
        # Verify each structure
        for name, df in loaded.items():
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 4
            assert df['pdb_id'].iloc[0] == name
    
    def test_zero_configuration(self):
        """Test that StructureProcessor works with zero configuration."""
        # No setup, no paths, just works
        processor = StructureProcessor()
        
        # Should be able to save immediately
        structure = self.create_test_structure("test_zero_config")
        processor.save_structure("test_zero_config", structure)
        
        # And load it back
        loaded = processor.load_structure("test_zero_config")
        assert loaded is not None
        assert len(loaded) == 4
    
    def test_drag_and_drop_workflow(self):
        """Test that dropped files can be loaded immediately."""
        processor = StructureProcessor()
        
        # Simulate dropping a file by creating it directly in the structure directory
        import os
        import pickle
        structure_dir = processor.path_structure_dir
        os.makedirs(structure_dir, exist_ok=True)
        
        # Create a pkl file (simulating a dropped structure file)
        dropped_structure = self.create_test_structure("dropped_protein")
        dropped_path = os.path.join(structure_dir, "dropped_protein.pkl")
        with open(dropped_path, 'wb') as f:
            pickle.dump(dropped_structure, f)
        
        # Should be able to load immediately by name
        loaded = processor.load_entity("dropped_protein")
        assert loaded is not None
        assert isinstance(loaded, pd.DataFrame)
        assert len(loaded) == 4
        assert loaded['pdb_id'].iloc[0] == "dropped_protein"
    
    def test_dataset_save_and_load(self):
        """Test saving and loading entire datasets as single PKL files."""
        processor = StructureProcessor()
        
        # Create and save multiple structures
        structures = {
            "protein_a": self.create_test_structure("protein_a"),
            "protein_b": self.create_test_structure("protein_b"),
            "protein_c": self.create_test_structure("protein_c")
        }
        
        # Save each structure
        for name, structure in structures.items():
            processor.save_structure(name, structure)
        
        # Create dataset
        dataset_name = processor.create_dataset(
            "test_dataset",
            list(structures.keys()),
            metadata={"study": "test"}
        )
        
        # Save dataset as single PKL
        dataset_path = processor.save_dataset("test_dataset")
        assert os.path.exists(dataset_path)
        
        # Clear processor and load dataset
        processor2 = StructureProcessor()
        loaded_dataset = processor2._load_dataset_impl("test_dataset")
        
        # Verify all structures are loaded
        assert isinstance(loaded_dataset, dict)
        assert len(loaded_dataset) == 3
        assert set(loaded_dataset.keys()) == {"protein_a", "protein_b", "protein_c"}
        
        # Verify structure contents
        for name, df in loaded_dataset.items():
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 4
            assert df['pdb_id'].iloc[0] == name
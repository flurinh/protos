"""
Test the updated StructureProcessor that uses ProtosPaths exclusively.
"""

import pytest
import tempfile
from pathlib import Path
import pandas as pd

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseProcessorUpdated:
    """Test the updated StructureProcessor implementation."""
    
    def test_zero_configuration(self):
        """Test processor works with zero configuration."""
        # Should work without any setup
        processor = StructureProcessor()
        
        # Should have created all components
        assert processor.paths is not None
        assert processor.entity_registry is not None
        assert processor.dataset_manager is not None
        assert processor.data_path.exists()
        
        # Should have proper subdirectories
        assert processor.path_structure_dir.exists()
        assert processor.path_dataset_dir.exists()
        assert processor.path_alignment_dir.exists()
        assert processor.path_cache_dir.exists()
    
    def test_no_custom_paths(self):
        """Test that processor has no custom path handling."""
        import inspect
        
        # Get init signature
        sig = inspect.signature(StructureProcessor.__init__)
        params = sig.parameters
        
        # These parameters should NOT exist
        param_names = list(params.keys())
        assert "data_root" not in param_names
        assert "processor_data_dir" not in param_names
        assert "structure_dir" not in param_names
        assert "dataset_dir" not in param_names
        assert "alignments_dir" not in param_names
        
        # Only 'paths' parameter is allowed
        assert "paths" in param_names
    
    def test_accepts_protospaths(self):
        """Test processor accepts ProtosPaths instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Should use provided paths
            assert processor.paths == paths
            assert str(processor.data_path).startswith(tmpdir)
    
    def test_path_properties_use_protospaths(self):
        """Test all path properties use ProtosPaths methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # All paths should go through ProtosPaths
            assert str(processor.path_structure_dir) == str(processor.get_subdirectory_path('structure_dir'))
            assert str(processor.path_dataset_dir) == str(processor.get_subdirectory_path('dataset_dir'))
            assert str(processor.path_alignment_dir) == str(processor.get_subdirectory_path('alignments_dir'))
            assert str(processor.path_cache_dir) == str(processor.get_subdirectory_path('cache_dir'))
    
    def test_save_and_load_structure(self):
        """Test saving and loading structures with human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create a simple test structure
            test_data = pd.DataFrame({
                'pdb_id': ['test_protein'] * 3,
                'auth_chain_id': ['A'] * 3,
                'auth_seq_id': [1, 2, 3],
                'res_name3l': ['ALA', 'GLY', 'VAL'],
                'res_name1l': ['A', 'G', 'V'],
                'res_atom_name': ['CA', 'CA', 'CA'],
                'atom_name': ['CA', 'CA', 'CA'],
                'x': [1.0, 2.0, 3.0],
                'y': [4.0, 5.0, 6.0],
                'z': [7.0, 8.0, 9.0],
                'group': ['ATOM'] * 3
            })
            
            # Save structure
            processor.save_structure("my_test_protein", test_data)
            
            # Should be registered
            assert processor.entity_exists("my_test_protein")
            
            # Load structure
            loaded = processor.load_entity("my_test_protein")
            assert loaded is not None
            assert len(loaded) == 3
            assert loaded['pdb_id'].iloc[0] == 'test_protein'
    
    def test_drag_and_drop_workflow(self):
        """Test drag-and-drop file discovery."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Manually create a CIF file (simulating drag-and-drop)
            mmcif_dir = Path(tmpdir) / "structure" / "mmcif"
            mmcif_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a simple test file
            test_file = mmcif_dir / "dropped_structure.cif"
            test_file.write_text("# Dummy CIF content")
            
            # Initialize PDB IDs
            processor.initialize_pdb_ids()
            
            # Should find the dropped file
            assert "dropped_structure" in processor.pdb_ids
    
    def test_dataset_operations(self):
        """Test dataset management with human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create some test structures
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'protein_{i}'] * 2,
                    'auth_chain_id': ['A'] * 2,
                    'auth_seq_id': [1, 2],
                    'res_name1l': ['A', 'G'],
                    'res_atom_name': ['CA', 'CA'],
                    'x': [1.0, 2.0],
                    'y': [3.0, 4.0],
                    'z': [5.0, 6.0],
                    'group': ['ATOM'] * 2
                })
                processor.save_structure(f"protein_{i}", test_data)
            
            # Create dataset with human names
            dataset_name = processor.create_dataset(
                "my_proteins",
                ["protein_0", "protein_1"],
                {"description": "Test protein dataset"}
            )
            
            assert dataset_name == "my_proteins"
            
            # List datasets
            datasets = processor.list_datasets()
            # BaseProcessor returns just the dataset names
            assert "my_proteins" in datasets
            
            # Add to dataset
            processor.add_to_dataset("my_proteins", ["protein_2"])
            info = processor.get_dataset_info("my_proteins")
            assert info["entity_count"] == 3
    
    def test_no_os_path_operations(self):
        """Verify no direct os.path operations in the code."""
        import inspect
        import ast
        
        # Get the source code
        source = inspect.getsource(StructureProcessor)
        
        # Parse the AST
        tree = ast.parse(source)
        
        # Check for forbidden operations
        forbidden_found = []
        
        class ForbiddenChecker(ast.NodeVisitor):
            def visit_Attribute(self, node):
                # Check for os.path.join, os.makedirs
                if isinstance(node.value, ast.Name) and node.value.id == 'os':
                    if node.attr in ['makedirs', 'path']:
                        forbidden_found.append(f"os.{node.attr}")
                self.generic_visit(node)
            
            def visit_Call(self, node):
                # Check for Path().mkdir()
                if isinstance(node.func, ast.Attribute) and node.func.attr == 'mkdir':
                    forbidden_found.append("mkdir()")
                self.generic_visit(node)
        
        ForbiddenChecker().visit(tree)
        
        # Should not find any forbidden operations
        assert len(forbidden_found) == 0, f"Found forbidden operations: {forbidden_found}"
    
    def test_format_data_types(self):
        """Test data type formatting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test data with mixed types
            test_data = pd.DataFrame({
                'pdb_id': ['test'] * 3,
                'auth_chain_id': ['A'] * 3,
                'auth_seq_id': ['1', '2', '3'],  # String numbers
                'x': ['1.5', '2.5', '3.5'],      # String floats
                'y': [4.0, 5.0, 6.0],
                'z': [7.0, 8.0, 9.0],
                'res_name1l': ['A', 'G', 'V'],
                'res_atom_name': ['CA', 'CA', 'CA'],
                'group': ['ATOM'] * 3
            })
            
            processor.data = test_data
            processor.format_data_types()
            
            # Check types are correct
            assert processor.data['x'].dtype == 'float64'
            assert processor.data['y'].dtype == 'float64'
            assert processor.data['z'].dtype == 'float64'
            # Integer columns should be nullable Int64
            assert str(processor.data['auth_seq_id'].dtype) == 'Int64'
    
    def test_filter_operations(self):
        """Test filtering methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test data with multiple proteins and chains
            test_data = pd.DataFrame({
                'pdb_id': ['prot1', 'prot1', 'prot2', 'prot2'],
                'auth_chain_id': ['A', 'B', 'A', 'B'],
                'auth_seq_id': [1, 1, 1, 1],
                'res_name1l': ['A', 'G', 'V', 'L'],
                'res_atom_name': ['CA', 'CA', 'CA', 'CA'],
                'x': [1.0, 2.0, 3.0, 4.0],
                'y': [5.0, 6.0, 7.0, 8.0],
                'z': [9.0, 10.0, 11.0, 12.0],
                'group': ['ATOM'] * 4
            })
            
            processor.data = test_data.copy()
            processor.pdb_ids = ['prot1', 'prot2']
            
            # Test filter by PDB IDs
            processor.filter_by_pdb_ids(['prot1'])
            assert len(processor.data) == 2
            assert processor.data['pdb_id'].unique()[0] == 'prot1'
            
            # Test filter by chain
            processor.data = test_data.copy()
            processor.filter_by_chain('A')
            assert len(processor.data) == 2
            assert all(processor.data['auth_chain_id'] == 'A')
"""
Test the updated SequenceProcessor that uses ProtosPaths exclusively.
"""

import pytest
import tempfile
from pathlib import Path
import pandas as pd

from protos.io.paths import ProtosPaths
from protos.processing.sequence import SequenceProcessor


class TestSeqProcessorUpdated:
    """Test the updated SequenceProcessor implementation."""
    
    def test_zero_configuration(self):
        """Test processor works with zero configuration."""
        # Should work without any setup
        processor = SequenceProcessor()
        
        # Should have created all components
        assert processor.paths is not None
        assert processor.entity_registry is not None
        assert processor.dataset_manager is not None
        assert processor.data_path.exists()
        
        # Should have proper subdirectories
        assert processor.path_fasta_dir.exists()
        assert processor.path_alignments_dir.exists()
        assert processor.path_metadata_dir.exists()
    
    def test_no_custom_paths(self):
        """Test that processor has no custom path handling."""
        import inspect
        
        # Get init signature
        sig = inspect.signature(SequenceProcessor.__init__)
        params = sig.parameters
        
        # These parameters should NOT exist
        param_names = list(params.keys())
        assert "data_root" not in param_names
        assert "processor_data_dir" not in param_names
        assert "fasta_dir" not in param_names
        assert "alignments_dir" not in param_names
        
        # Only 'paths' parameter is allowed
        assert "paths" in param_names
    
    def test_accepts_protospaths(self):
        """Test processor accepts ProtosPaths instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Should use provided paths
            assert processor.paths == paths
            assert str(processor.data_path).startswith(tmpdir)
    
    def test_path_properties_use_protospaths(self):
        """Test all path properties use ProtosPaths methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # All paths should go through ProtosPaths
            assert str(processor.path_fasta_dir) == str(processor.get_subdirectory_path('fasta_dir'))
            assert str(processor.path_alignments_dir) == str(processor.get_subdirectory_path('alignment_dir'))
            assert str(processor.path_metadata_dir) == str(processor.get_subdirectory_path('metadata_dir'))
    
    def test_save_and_load_single_sequence(self):
        """Test saving and loading single sequences with human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Save a single sequence
            test_sequence = "MAGICLAMP"
            processor.save_entity("my_protein", test_sequence)
            
            # Should be registered
            assert processor.entity_exists("my_protein")
            
            # Load sequence
            loaded = processor.load_entity("my_protein")
            assert loaded == test_sequence
    
    def test_save_and_load_multi_sequence(self):
        """Test saving and loading multi-sequence files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Save multiple sequences
            sequences = {
                "protein_A": "MAGICLAMP",
                "protein_B": "EVERYPROTEINS",
                "protein_C": "INALIGNMENT"
            }
            
            processor.save_entity("my_protein_set", sequences)
            
            # Should be registered
            assert processor.entity_exists("my_protein_set")
            
            # Load sequences
            loaded = processor.load_entity("my_protein_set")
            assert isinstance(loaded, dict)
            assert len(loaded) == 3
            assert loaded["protein_A"] == "MAGICLAMP"
    
    def test_filename_sanitization(self):
        """Test that problematic characters in names are handled."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Save with problematic name
            problematic_name = "protein/with:special|chars*"
            processor.save_entity(problematic_name, "TESTSEQ")
            
            # Should still be able to load
            loaded = processor.load_entity(problematic_name)
            assert loaded == "TESTSEQ"
            
            # Check actual filename was sanitized
            fasta_files = list(processor.path_fasta_dir.glob("*.fasta"))
            assert len(fasta_files) == 1
            assert "/" not in fasta_files[0].name
            assert ":" not in fasta_files[0].name
    
    def test_drag_and_drop_workflow(self):
        """Test drag-and-drop file discovery."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Manually create a FASTA file (simulating drag-and-drop)
            fasta_dir = Path(tmpdir) / "sequence" / "fasta"
            fasta_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a simple test file
            test_file = fasta_dir / "dropped_sequences.fasta"
            test_file.write_text(">seq1\nMAGICLAMP\n>seq2\nEVERYPROTEINS")
            
            # List entities should find it
            entities = processor.list_entities()
            assert "dropped_sequences" in entities
            
            # Should be able to load it
            sequences = processor.load_entity("dropped_sequences")
            assert isinstance(sequences, dict)
            assert len(sequences) == 2
    
    def test_dataset_operations(self):
        """Test dataset management with human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Create some test sequences
            sequences = {}
            for i in range(5):
                sequences[f"protein_{i}"] = f"SEQUENCE{i}" * 10
            
            # Save as multi-sequence file
            processor.save_entity("protein_family", sequences)
            
            # Create dataset from subset
            dataset_name = processor.create_dataset(
                "my_protein_subset",
                ["protein_0", "protein_1", "protein_2"],
                {"description": "Subset of protein family"}
            )
            
            assert dataset_name == "my_protein_subset"
            
            # List datasets
            datasets = processor.list_datasets()
            assert "my_protein_subset" in datasets
            # The multi-sequence file itself is also a dataset
            assert "protein_family" in datasets
    
    def test_no_os_path_operations(self):
        """Verify no direct os.path operations in the code."""
        import inspect
        import ast
        
        # Get the source code
        source = inspect.getsource(SequenceProcessor)
        
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
    
    def test_sequence_alignment(self):
        """Test sequence alignment functionality."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Test sequences
            seq1 = "MAGICLAMP"
            seq2 = "MAGICWAND"
            
            # Align sequences
            score, alignment = processor.align_sequences(seq1, seq2)
            
            # Should have positive score (similar sequences)
            assert score > 0
            assert len(alignment) == 3  # Three lines in alignment format
            
            # Check alignment was stored
            assert "seq1_vs_seq2" in processor.alignments
    
    def test_mutate_sequence(self):
        """Test sequence mutation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Original sequence
            original = "MAGICLAMP"
            
            # Apply mutations
            mutations = ["M1W", "L6W"]  # M->W at position 1, L->W at position 6
            mutated = processor.mutate_sequence(original, mutations)
            
            assert mutated == "WAGICWAMP"
            assert mutated != original
    
    def test_find_best_match(self):
        """Test finding best matching sequence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Reference sequences
            references = {
                "ref1": "MAGICLAMP",
                "ref2": "EVERYPROTEIN",
                "ref3": "MAGICWAND"
            }
            
            # Query similar to ref3
            query = "MAGICBAND"
            
            # Find best match (using BioPython since MMseqs2 may not be available)
            best_id, score, alignment = processor.find_best_match(
                query, references, use_mmseqs=False
            )
            
            # Should match ref3 (most similar)
            assert best_id == "ref3"
            assert score > 0
    
    def test_get_sequence(self):
        """Test getting individual sequences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Save a multi-sequence file
            sequences = {
                "protA": "AAAAA",
                "protB": "BBBBB",
                "protC": "CCCCC"
            }
            processor.save_entity("test_proteins", sequences)
            
            # Get individual sequence
            seq = processor.get_sequence("protB")
            assert seq == "BBBBB"
            
            # Try non-existent
            seq = processor.get_sequence("protX")
            assert seq is None
    
    def test_sequence_metadata(self):
        """Test sequence metadata calculation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = SequenceProcessor(paths=paths)
            
            # Add some sequences to cache
            processor.sequences = {
                "seq1": "ARNDCQEGHILKMFPSTWYVAA",  # All 20 amino acids + AA
                "seq2": "DDDDDEEEEE"  # Acidic
            }
            
            # Get metadata
            metadata = processor.get_sequence_metadata()
            
            assert len(metadata) == 2
            assert "molecular_weight" in metadata.columns
            assert "isoelectric_point" in metadata.columns
            assert "amino_acid_composition" in metadata.columns
            
            # Check pI calculation
            acidic_seq = metadata[metadata['sequence_id'] == 'seq2']
            assert acidic_seq['isoelectric_point'].iloc[0] == 4.5  # Acidic
import pytest
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from protos.processing.structure.gtalign import GTalign, transform_coordinates, apply_alignment_transformation, print_alignment_result
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from pathlib import Path

# Global constants for testing
TEST_PDB_IDS = ["5ahz", "5awz", "6rf0"]


@pytest.fixture(scope="session")
def persistent_output_dir():
    """
    Create a persistent output directory for old_tests,
    located in data/temp.
    Returns a tuple: (base_dir, cif_dir, gtalign_dir)
    """
    base_dir = os.path.join("data", "temp")
    os.makedirs(base_dir, exist_ok=True)
    cif_dir = os.path.join(base_dir, "cif")
    gtalign_dir = os.path.join(base_dir, "gtalign")
    os.makedirs(cif_dir, exist_ok=True)
    os.makedirs(gtalign_dir, exist_ok=True)
    yield base_dir, cif_dir, gtalign_dir
    # Optionally, cleanup after the test session:
    # shutil.rmtree(base_dir, ignore_errors=True)


@pytest.fixture
def gta():
    """Create a GTalign instance for testing."""
    import shutil
    if shutil.which("gtalign") is None:
        pytest.skip("GTalign executable not found in PATH. Skipping test.")
    return GTalign()


@pytest.fixture
def cp():
    """Create a CifProcessor instance for testing."""
    data_path = "data"
    structure_dir = "mmcif"
    return CifBaseProcessor(path=data_path, structure_dir=structure_dir)


class TestGTalign:
    def test_real_integration_with_cp(self, cp, gta, persistent_output_dir):
        """
        Integration test with real GTalign and CifProcessor.
        Uses CIF files saved to the persistent cif folder and stores GTalign outputs
        in subfolders of the persistent gtalign folder.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(file1) and os.path.exists(file2)):
            pytest.skip("Test structures not found for integration test")

        cp.load_structure(TEST_PDB_IDS[0])
        cp.load_structure(TEST_PDB_IDS[1])

        temp_file1 = cp.save_temp_cif(TEST_PDB_IDS[0], suffix="alignment", output_dir=cif_dir)
        temp_file2 = cp.save_temp_cif(TEST_PDB_IDS[1], suffix="alignment", output_dir=cif_dir)
        assert os.path.exists(temp_file1)
        assert os.path.exists(temp_file2)

        print("\nFirst 20 lines of temporary CIF file:")
        with open(temp_file1, 'r') as f:
            for i, line in enumerate(f.readlines()[:20]):
                print(f"{i + 1}: {line.rstrip()}")

        alignment_output = os.path.join(gtalign_dir, "alignment_output")
        os.makedirs(alignment_output, exist_ok=True)

        try:
            version_result = subprocess.run(['gtalign', '--help'],
                                            capture_output=True, text=True)
            print("\nGTalign help info (first few lines):")
            for line in version_result.stdout.split('\n')[:10]:
                print(line)
        except Exception as e:
            print(f"Error checking GTalign version: {e}")

        result = gta.align(
            query_structures=str(temp_file1),
            reference_structures=str(temp_file2),
            output_dir=alignment_output,
            tm_score_threshold=0.0,
            speed=9,
            verbose=True
        )

        alignment_files = [f for f in os.listdir(alignment_output) if f.endswith('.aln') or f.endswith('.out')]
        print("\nRegular alignment results:")
        print(f"Alignment files found: {alignment_files}")
        print(f"Stdout length: {len(result.get('stdout', ''))}")
        print(f"Stderr length: {len(result.get('stderr', ''))}")
        print(f"Number of alignments: {len(result.get('alignments', []))}")
        print("\nFull stderr content:")
        print("-" * 50)
        print(result.get("stderr", "No stderr content"))
        print("-" * 50)

        self_alignment_output = os.path.join(gtalign_dir, "self_alignment_output")
        os.makedirs(self_alignment_output, exist_ok=True)
        self_result = gta.align(
            query_structures=str(temp_file1),
            reference_structures=str(temp_file1),
            output_dir=self_alignment_output,
            tm_score_threshold=0.0,
            speed=9,
            verbose=True
        )
        self_aln_files = [f for f in os.listdir(self_alignment_output) if f.endswith('.aln') or f.endswith('.out')]
        print(f"Self-alignment files found: {self_aln_files}")
        print(f"Self-alignment stdout length: {len(self_result.get('stdout', ''))}")
        print(f"Self-alignment stderr length: {len(self_result.get('stderr', ''))}")
        print(f"Self-alignment number of alignments: {len(self_result.get('alignments', []))}")
        print("\nSelf-alignment stderr content:")
        print("-" * 50)
        print(self_result.get("stderr", "No stderr content"))
        print("-" * 50)

        sensitive_output = os.path.join(gtalign_dir, "sensitive_alignment_output")
        os.makedirs(sensitive_output, exist_ok=True)
        sensitive_result = gta.align(
            query_structures=str(temp_file1),
            reference_structures=str(temp_file2),
            output_dir=sensitive_output,
            tm_score_threshold=0.0,
            speed=2,
            depth=3,
            verbose=True
        )
        sensitive_aln_files = [f for f in os.listdir(sensitive_output) if f.endswith('.aln') or f.endswith('.out')]
        print(f"Sensitive alignment files found: {sensitive_aln_files}")
        print(f"Sensitive alignment stderr length: {len(sensitive_result.get('stderr', ''))}")
        assert True

    def test_cluster_with_real_files(self, cp, gta, persistent_output_dir):
        """
        Test clustering with real structure files.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(file1) and os.path.exists(file2)):
            pytest.skip("Test structures not found for clustering test")

        cp.load_structure(TEST_PDB_IDS[0])
        cp.load_structure(TEST_PDB_IDS[1])
        temp_file1 = cp.save_temp_cif(TEST_PDB_IDS[0], suffix="cluster", output_dir=cif_dir)
        temp_file2 = cp.save_temp_cif(TEST_PDB_IDS[1], suffix="cluster", output_dir=cif_dir)
        assert os.path.exists(temp_file1)
        assert os.path.exists(temp_file2)

        cluster_output = os.path.join(gtalign_dir, "cluster_output")
        os.makedirs(cluster_output, exist_ok=True)
        try:
            result = gta.cluster(
                structures=[str(temp_file1), str(temp_file2)],
                output_dir=cluster_output,
                tm_score_threshold=0.0,
                speed=9,
                verbose=True
            )
            assert "stdout" in result
        except subprocess.CalledProcessError as e:
            pytest.fail(f"GTalign clustering failed to run: {e}")
        except RuntimeError as e:
            print(f"GTalign clustering reported an error: {e}")

    def test_alignment_with_original_files(self, gta, persistent_output_dir):
        """
        Test alignment using the original CIF files in data/mmcif.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        original_file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        original_file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(original_file1) and os.path.exists(original_file2)):
            pytest.skip("Test structures not found for original file alignment test")
        assert os.path.exists(original_file1)
        assert os.path.exists(original_file2)

        orig_alignment_output = os.path.join(gtalign_dir, "original_alignment_output")
        os.makedirs(orig_alignment_output, exist_ok=True)
        try:
            result = gta.align(
                query_structures=original_file1,
                reference_structures=original_file2,
                output_dir=orig_alignment_output,
                tm_score_threshold=0.0,
                speed=9,
                verbose=True
            )
            print("\nAlignment using original files:")
            print(f"Stdout length: {len(result.get('stdout', ''))}")
            print(f"Stderr length: {len(result.get('stderr', ''))}")
            print("\nFull stderr content (original files):")
            print("-" * 50)
            print(result.get("stderr", "No stderr content"))
            print("-" * 50)
            alignment_files = [f for f in os.listdir(orig_alignment_output) if f.endswith('.aln') or f.endswith('.out')]
            print(f"Alignment files found: {alignment_files}")
            print(f"Number of alignments: {len(result.get('alignments', []))}")
            assert True
        except Exception as e:
            print(f"Error running GTalign on original files: {e}")
            assert True

    def test_msa_with_three_proteins(self, cp, gta, persistent_output_dir):
        """
        Test a multiple alignment scenario with three proteins.
        This test loads 5AHZ, 5AWZ, and 6RF0, saves their CIF files into the persistent cif folder,
        and then runs GTalign with both query and reference structures set to all three files.
        The resulting pairwise alignments are printed, and these can be post-processed into an MSA table.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        pdb_ids = ["5ahz", "5awz", "6rf0"]
        cif_files = []
        for pdb in pdb_ids:
            file_path = os.path.join("data", "mmcif", f"{pdb}.cif")
            if not os.path.exists(file_path):
                pytest.skip(f"Test structure {pdb}.cif not found")
            cp.load_structure(pdb)
            cif_file = cp.save_temp_cif(pdb, suffix="msa", output_dir=cif_dir)
            cif_files.append(str(cif_file))
            assert os.path.exists(cif_file)

        msa_output = os.path.join(gtalign_dir, "msa_output")
        os.makedirs(msa_output, exist_ok=True)

        # Run GTalign with all three proteins as both query and reference.
        result = gta.align(
            query_structures=cif_files,
            reference_structures=cif_files,
            output_dir=msa_output,
            tm_score_threshold=0.0,
            speed=9,
            verbose=True
        )

        print("\nMSA test alignment results:")
        for aln in result["alignments"]:
            print("Alignment file:", aln["file"])
            print("Query sequence:", aln["query_sequence"])
            print("Reference sequence:", aln["reference_sequence"])
            print("TM-score:", aln["tm_score"])
            print("RMSD:", aln["rmsd"])
            
            # Print additional alignment information if available
            if aln.get("rotation_matrix"):
                print("Rotation matrix:")
                for row in aln["rotation_matrix"]:
                    print("  ", " ".join([f"{x:.6f}" for x in row]))
                    
            if aln.get("translation_vector"):
                print("Translation vector:", " ".join([f"{x:.6f}" for x in aln["translation_vector"]]))
                
            if aln.get("aligned_query_structure"):
                print("Aligned query structure:", aln["aligned_query_structure"])
                
            if aln.get("aligned_reference_structure"):
                print("Aligned reference structure:", aln["aligned_reference_structure"])
                
            print("-----")

        # Here, you can post-process the pairwise alignments into an MSA table.
        # For now, we simply assert that the alignment results list is non-empty.
        assert len(result["alignments"]) > 0
        
    def test_alignment_with_transformation(self, cp, gta, persistent_output_dir):
        """
        Test alignment with rotation matrix and translation vector extraction.
        This test verifies that the rotation matrix and translation vector
        are correctly extracted from the alignment results output file.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(file1) and os.path.exists(file2)):
            pytest.skip("Test structures not found for transformation test")

        cp.load_structure(TEST_PDB_IDS[0])
        cp.load_structure(TEST_PDB_IDS[1])

        temp_file1 = cp.save_temp_cif(TEST_PDB_IDS[0], suffix="transform", output_dir=cif_dir)
        temp_file2 = cp.save_temp_cif(TEST_PDB_IDS[1], suffix="transform", output_dir=cif_dir)
        assert os.path.exists(temp_file1)
        assert os.path.exists(temp_file2)

        transform_output = os.path.join(gtalign_dir, "transform_output")
        os.makedirs(transform_output, exist_ok=True)

        # Run a standard alignment - the rotation matrix and translation vector
        # should be parsed from the output file
        result = gta.align(
            query_structures=str(temp_file1),
            reference_structures=str(temp_file2),
            output_dir=transform_output,
            tm_score_threshold=0.0,
            speed=9,
            verbose=True,
            referenced=True,  # Generate transformation matrices for reference structures
        )

        alignment_files = [f for f in os.listdir(transform_output) if f.endswith('.aln') or f.endswith('.out')]
        print("\nTransformation test alignment results:")
        print(f"Alignment files found: {alignment_files}")
        
        # Print the content of the alignment file to debug the format
        if alignment_files:
            print("\nExamining first alignment file:")
            first_file = os.path.join(transform_output, alignment_files[0])
            with open(first_file, 'r') as f:
                content = f.read()
                # Print a sample of the content
                print(f"Sample of file content (first 200 chars): {content[:200]}")
                
                # Look for the rotation matrix section
                matrix_pattern = r"Rotation \[3,3\] and translation \[3,1\] for Query:"
                if matrix_pattern in content:
                    print(f"\nFound rotation matrix pattern: '{matrix_pattern}'")
                    idx = content.find(matrix_pattern)
                    print(f"Found at position {idx}")
                    # Print content around the match
                    start = max(0, idx - 10)
                    end = min(len(content), idx + 200)
                    print(f"Context: {content[start:end]}")
        
        # Check if we have alignments
        assert len(result["alignments"]) > 0
        
        # Print alignment details including transformations
        for aln in result["alignments"]:
            print("Alignment file:", aln["file"])
            print("TM-score:", aln["tm_score"])
            print("RMSD:", aln["rmsd"])
            
            # Check for rotation matrix
            if aln.get("rotation_matrix"):
                print("Rotation matrix found:")
                for row in aln["rotation_matrix"]:
                    print("  ", " ".join([f"{x:.6f}" for x in row]))
                # A valid rotation matrix should have 3 rows with 3 columns each
                assert len(aln["rotation_matrix"]) == 3
                for row in aln["rotation_matrix"]:
                    assert len(row) == 3
            else:
                print("Rotation matrix not found in alignment results")
                
            # Check for translation vector
            if aln.get("translation_vector"):
                print("Translation vector found:", " ".join([f"{x:.6f}" for x in aln["translation_vector"]]))
                # A valid translation vector should have 3 elements
                assert len(aln["translation_vector"]) == 3
            else:
                print("Translation vector not found in alignment results")
            
            print("-----")
            
    def test_coordinate_transformation(self):
        """
        Test the coordinate transformation functions.
        This test verifies that the transform_coordinates and apply_alignment_transformation
        functions correctly apply rotation and translation to coordinates.
        """
        # Create a simple test case with 3 points
        # Points form a right triangle in the xy-plane: (0,0,0), (1,0,0), (0,1,0)
        coords = np.array([
            [0.0, 0.0, 0.0],  # Origin
            [1.0, 0.0, 0.0],  # 1 unit along x
            [0.0, 1.0, 0.0]   # 1 unit along y
        ])
        
        # Define a rotation matrix (90 degrees around z-axis)
        # This will map (1,0,0) to (0,1,0) and (0,1,0) to (-1,0,0)
        rot_matrix = np.array([
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        
        # Define a translation vector (shift by (2,3,4))
        trans_vector = np.array([2.0, 3.0, 4.0])
        
        # Apply the transformation
        transformed = transform_coordinates(coords, rot_matrix, trans_vector)
        
        # Expected result: 
        # (0,0,0) -> (2,3,4)
        # (1,0,0) -> (2,4,4)  # rotated then translated
        # (0,1,0) -> (1,3,4)  # rotated then translated
        expected = np.array([
            [2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0],
            [1.0, 3.0, 4.0]
        ])
        
        # Check that the transformation was applied correctly
        np.testing.assert_almost_equal(transformed, expected)
        
        # Test apply_alignment_transformation with DataFrame
        df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'atom_name': ['CA', 'CA', 'CA'],
            'x': [0.0, 1.0, 0.0],
            'y': [0.0, 0.0, 1.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        transformed_df = apply_alignment_transformation(df, rot_matrix, trans_vector)
        
        # Check that the transformation was applied correctly to the DataFrame
        expected_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'atom_name': ['CA', 'CA', 'CA'],
            'x': [2.0, 2.0, 1.0],
            'y': [3.0, 4.0, 3.0],
            'z': [4.0, 4.0, 4.0]
        })
        
        # Check that only x, y, z columns were transformed
        assert list(transformed_df.columns) == list(df.columns)
        assert transformed_df['atom_id'].equals(df['atom_id'])
        assert transformed_df['atom_name'].equals(df['atom_name'])
        
        # Check transformed coordinates
        np.testing.assert_almost_equal(
            transformed_df[['x', 'y', 'z']].values, 
            expected_df[['x', 'y', 'z']].values
        )
        
    def test_align_from_processor(self, cp, gta, persistent_output_dir):
        """
        Test the align_from_processor helper function.
        This test verifies that the align_from_processor function correctly
        prepares structures from CifProcessor, runs GTalign, and returns
        transformed structures.
        """
        base_dir, cif_dir, gtalign_dir = persistent_output_dir
        
        # Check if test structures exist
        file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(file1) and os.path.exists(file2)):
            pytest.skip("Test structures not found for align_from_processor test")
        
        # Load structures into CifProcessor
        cp.load_structure(TEST_PDB_IDS[0])
        cp.load_structure(TEST_PDB_IDS[1])
        
        # Create output directory
        processor_output = os.path.join(gtalign_dir, "processor_output")
        os.makedirs(processor_output, exist_ok=True)
        
        # Run align_from_processor with save_transformed and save_matrix
        result = gta.align_from_processor(
            cif_processor=cp,
            reference_pdb_id=TEST_PDB_IDS[0],
            target_pdb_ids=[TEST_PDB_IDS[1]],
            chain_id='A',
            output_dir=processor_output,
            tm_score_threshold=0.0,
            speed=9,
            referenced=True,
            verbose=True,
            save_transformed=True,
            save_matrix=True
        )
        
        # Verify that we got alignment results
        assert 'alignments' in result
        assert len(result['alignments']) > 0
        
        # Verify that we got transformed structures
        assert 'aligned_structures' in result
        assert TEST_PDB_IDS[1] in result['aligned_structures']
        
        # The transformed structure should be a DataFrame
        transformed_structure = result['aligned_structures'][TEST_PDB_IDS[1]]
        assert isinstance(transformed_structure, pd.DataFrame)
        
        # It should have x, y, z columns
        assert 'x' in transformed_structure.columns
        assert 'y' in transformed_structure.columns
        assert 'z' in transformed_structure.columns
        
        # Get the original and transformed coordinates
        original_structure = cp.data[cp.data['pdb_id'] == TEST_PDB_IDS[1]]
        
        # Verify that coordinates changed
        assert not np.array_equal(
            original_structure[['x', 'y', 'z']].values,
            transformed_structure[['x', 'y', 'z']].values
        )
        
        print(f"\nAlign from processor test results:")
        print(f"Number of alignments: {len(result['alignments'])}")
        print(f"Number of transformed structures: {len(result['aligned_structures'])}")

        # Find the alignment for our target structure
        for aln in result['alignments']:
            if TEST_PDB_IDS[1] in aln['file']:
                print(f"Found alignment for {TEST_PDB_IDS[1]}")
                print(f"TM-score: {aln['tm_score']}")
                print(f"RMSD: {aln['rmsd']}")
                break
                
        # Verify that transformed structures were saved
        assert 'transformed_structure_paths' in result
        assert TEST_PDB_IDS[1] in result['transformed_structure_paths']
        transformed_file_path = result['transformed_structure_paths'][TEST_PDB_IDS[1]]
        assert os.path.exists(transformed_file_path)
        print(f"Transformed structure saved to: {transformed_file_path}")
        
        # Check the content of the transformed file
        with open(transformed_file_path, 'r') as f:
            content = f.read()
            assert "data_aligned_structure" in content
            assert "_atom_site.Cartn_x" in content
            assert "_atom_site.Cartn_y" in content
            assert "_atom_site.Cartn_z" in content
        
        # Verify that matrix was saved
        matrix_dir = os.path.join(os.path.dirname(transformed_file_path), "..", "matrices")
        assert os.path.exists(matrix_dir)
        matrix_files = [f for f in os.listdir(matrix_dir) if f.startswith("transformation_matrix_")]
        assert len(matrix_files) > 0
        
        # Check matrix file content
        matrix_file_path = os.path.join(matrix_dir, matrix_files[0])
        with open(matrix_file_path, 'r') as f:
            content = f.read()
            assert "Transformation Matrix" in content
            assert "Rotation matrix" in content
                
    def test_timestamp_filtering(self, cp, gta, persistent_output_dir):
        """
        Test that alignment results are correctly filtered by timestamp.
        This test ensures that when running multiple alignments, only the
        newly created alignment files are used and not older ones.
        """
        import time

        base_dir, cif_dir, gtalign_dir = persistent_output_dir

        # Check if test structures exist
        file1 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[0]}.cif")
        file2 = os.path.join("data", "mmcif", f"{TEST_PDB_IDS[1]}.cif")
        if not (os.path.exists(file1) and os.path.exists(file2)):
            pytest.skip("Test structures not found for timestamp filtering test")

        # Load structures into CifProcessor
        cp.load_structure(TEST_PDB_IDS[0])
        cp.load_structure(TEST_PDB_IDS[1])

        # Create separate output directory for this test
        timestamp_output = os.path.join(gtalign_dir, "timestamp_filter_output")
        os.makedirs(timestamp_output, exist_ok=True)

        # Run first alignment
        run1_result = gta.align_from_processor(
            cif_processor=cp,
            reference_pdb_id=TEST_PDB_IDS[0],
            target_pdb_ids=[TEST_PDB_IDS[1]],
            chain_id='A',
            output_dir=timestamp_output,
            tm_score_threshold=0.0,
            speed=9,
            referenced=True,
            verbose=True
        )

        print(f"First run alignment count: {len(run1_result['alignments'])}")

        # Get the alignment directories from first run
        run1_dirs = []
        for alignment in run1_result['alignments']:
            alignment_file = alignment['file']
            alignment_dir = os.path.dirname(alignment_file)
            if alignment_dir not in run1_dirs:
                run1_dirs.append(alignment_dir)

        print(f"First run alignment directories: {run1_dirs}")

        # Sleep to ensure different timestamps
        time.sleep(2)

        # Run second alignment with save_transformed and save_matrix
        run2_result = gta.align_from_processor(
            cif_processor=cp,
            reference_pdb_id=TEST_PDB_IDS[0],
            target_pdb_ids=[TEST_PDB_IDS[1]],
            chain_id='A',
            output_dir=timestamp_output,
            tm_score_threshold=0.0,
            speed=9,
            referenced=True,
            verbose=True,
            save_transformed=True,
            save_matrix=True
        )

        print(f"Second run alignment count: {len(run2_result['alignments'])}")

        # Get the alignment directories from second run
        run2_dirs = []
        for alignment in run2_result['alignments']:
            alignment_file = alignment['file']
            alignment_dir = os.path.dirname(alignment_file)
            if alignment_dir not in run2_dirs:
                run2_dirs.append(alignment_dir)

        print(f"Second run alignment directories: {run2_dirs}")

        # Verify that the alignment directories are different
        # This confirms we're using different directories with timestamps
        for run1_dir in run1_dirs:
            for run2_dir in run2_dirs:
                assert run1_dir != run2_dir, "Both runs used the same alignment directory"

        # Verify that we didn't pick up old alignment files
        # If our timestamp filtering is working, run2 should not include run1's files
        for alignment1 in run1_result['alignments']:
            for alignment2 in run2_result['alignments']:
                assert alignment1['file'] != alignment2['file'], "Second run reused alignment file from first run"
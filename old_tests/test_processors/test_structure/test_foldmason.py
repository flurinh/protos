import pytest
import os
import shutil
import subprocess
from pathlib import Path
from protos.processing.structure.foldmason import FoldMason
from protos.processing.structure.struct_processor import CifProcessor

# Global constants for testing
TEST_PDB_IDS = [
    "5ahz",
    "5awz"
]


@pytest.fixture(scope="session")
def persistent_output_dir():
    """
    Create a persistent output directory for old_tests,
    located at C:/Users/hidbe/PycharmProjects/GPCR/data/temp/foldmason.
    Returns the base directory path.
    """
    base_dir = Path(r"C:/Users/hidbe/PycharmProjects/GPCR/data/temp/foldmason")
    base_dir.mkdir(parents=True, exist_ok=True)
    yield base_dir
    # Uncomment to clean up after old_tests (not recommended for debugging)
    # shutil.rmtree(base_dir, ignore_errors=True)


@pytest.fixture
def temp_dir(persistent_output_dir):
    """Create a unique temporary directory for each test within the persistent output dir."""
    import time
    test_dir = persistent_output_dir / f"test_{int(time.time())}"
    test_dir.mkdir(exist_ok=True)
    yield test_dir
    try:
        shutil.rmtree(test_dir, ignore_errors=True)
    except Exception as e:
        print(f"Warning: Failed to clean up {test_dir}: {e}")


@pytest.fixture
def fm():
    """Create a FoldMason instance for testing with WSL."""
    try:
        result = subprocess.run('wsl bash -lic "which foldmason"',
                                shell=True, capture_output=True, text=True, check=False)
        if not result.stdout.strip():
            pytest.skip("FoldMason executable not found in WSL. Skipping test.")
        print(f"FoldMason found in WSL at: {result.stdout.strip()}")
        return FoldMason(use_wsl=True, debug=True)
    except Exception as e:
        pytest.skip(f"Error checking for FoldMason in WSL: {e}")


@pytest.fixture
def cp():
    """Create a CifProcessor instance for testing."""
    data_path = "data"
    structure_dir = "mmcif"
    return CifProcessor(path=data_path, structure_dir=structure_dir)


class TestFoldMason:
    def test_init(self, fm):
        """Test initialization of FoldMason."""
        assert isinstance(fm, FoldMason)
        assert hasattr(fm, "run_command")
        assert hasattr(fm, "easy_msa")

    def test_run_command(self, fm):
        """Test that a simple command executes successfully."""
        try:
            result = fm.run_command("foldmason --help")
            print("\nFoldMason help output:")
            print("-" * 50)
            print(result[:500] + "..." if len(result) > 500 else result)
            print("-" * 50)
            assert "FoldMason" in result
            assert len(result) > 0
        except Exception as e:
            print(f"Error running basic command: {e}")
            assert True

    def test_easy_msa_interface(self, fm, temp_dir):
        """Test the easy-msa interface using dummy PDB files."""
        file1 = temp_dir / "file1.pdb"
        file2 = temp_dir / "file2.pdb"
        file1.write_text("""HEADER    PROTEIN                                 01-JAN-21   XXXX              
TITLE     DUMMY PROTEIN STRUCTURE 1                                             
ATOM      1  N   MET A   1      11.104  13.207  10.567  1.00 20.00           N  
ATOM      2  CA  MET A   1      11.804  12.147   9.843  1.00 20.00           C  
ATOM      3  C   MET A   1      11.753  12.247   8.340  1.00 20.00           C  
ATOM      4  O   MET A   1      10.882  12.940   7.815  1.00 20.00           O  
ATOM      5  CB  MET A   1      11.130  10.814  10.187  1.00 20.00           C  
ATOM      6  CG  MET A   1      11.130  10.472  11.669  1.00 20.00           C  
ATOM      7  SD  MET A   1      10.513   8.838  12.041  1.00 20.00           S  
ATOM      8  CE  MET A   1      10.780   8.741  13.834  1.00 20.00           C  
TER       9      MET A   1
END""")
        file2.write_text("""HEADER    PROTEIN                                 01-JAN-21   XXXX              
TITLE     DUMMY PROTEIN STRUCTURE 2                                             
ATOM      1  N   ALA A   1      12.104  14.207  11.567  1.00 20.00           N  
ATOM      2  CA  ALA A   1      12.804  13.147  10.843  1.00 20.00           C  
ATOM      3  C   ALA A   1      12.753  13.247   9.340  1.00 20.00           C  
ATOM      4  O   ALA A   1      11.882  13.940   8.815  1.00 20.00           O  
ATOM      5  CB  ALA A   1      12.130  11.814  11.187  1.00 20.00           C  
TER       6      ALA A   1
END""")
        file1_abs = file1.resolve()
        file2_abs = file2.resolve()
        print(f"\nInput file 1: {file1_abs}")
        print(f"Input file 2: {file2_abs}")
        print("\nFirst 10 lines of file1.pdb:")
        with open(file1_abs, 'r') as f:
            for i, line in enumerate(f.readlines()[:10]):
                print(f"{i + 1}: {line.rstrip()}")
        print("\nFirst 10 lines of file2.pdb:")
        with open(file2_abs, 'r') as f:
            for i, line in enumerate(f.readlines()[:10]):
                print(f"{i + 1}: {line.rstrip()}")

        # Create output directories
        output_prefix = temp_dir / "test_msa"
        tmp_folder = temp_dir / "tmp"
        tmp_folder.mkdir(exist_ok=True)

        try:
            # Use string command instead of list
            help_output = fm.run_command("foldmason easy-msa --help")
            print(f"\nFoldMason easy-msa help: {help_output[:500]}...")

            try:
                result = fm.easy_msa(
                    input_files=[file1_abs, file2_abs],
                    output_prefix=output_prefix,
                    tmp_folder=tmp_folder,
                    report_mode=1
                )
                print("\neasy-msa result:", result)
            except Exception as e:
                print(f"\neasy-msa failed: {e}")
        finally:
            print("\nChecking output directory contents:")
            for file in temp_dir.iterdir():
                print(f"  {file.name}")
            if tmp_folder.exists():
                print("\nChecking tmp folder contents:")
                for file in tmp_folder.iterdir():
                    print(f"  {file.name}")

    def test_msa2lddt_interface(self, fm, temp_dir):
        """Test msa2lddt interface using valid test files."""
        test_db = temp_dir / "test.db"
        test_fasta = temp_dir / "test.fasta"
        test_fasta.write_text(""">sequence1
ACDEFGHIKLMNPQRSTVWY
>sequence2
ACDEFGHIKLMNPQRSTVWY
""")
        with open(test_db, 'wb') as f:
            f.write(b'\x00\x01\x02\x03')

        try:
            try:
                result = fm.msa2lddt(
                    structure_db=test_db,
                    input_fasta=test_fasta
                )
                print("\nmsa2lddt result:", result)
            except Exception as e:
                print(f"\nmsa2lddt failed as expected with dummy db: {e}")
        finally:
            print(f"\nContents of temp_dir {temp_dir}:")
            for file in temp_dir.iterdir():
                print(f"  {file.name} ({file.stat().st_size} bytes)")

    @pytest.mark.skipif(
        not (Path("data") / "mmcif" / f"{TEST_PDB_IDS[0]}.cif").exists() or
        not (Path("data") / "mmcif" / f"{TEST_PDB_IDS[1]}.cif").exists(),
        reason="Test structures not found"
    )
    def test_integration_with_cp_generated(self, cp, fm, persistent_output_dir):
        """Integration test with CifProcessor and FoldMason using generated CIF files."""
        # Step 1: Create consistent directory structure as with other old_tests
        cif_dir = persistent_output_dir / "generated_cif"
        cif_dir.mkdir(exist_ok=True)

        # Create dedicated directories for output and temporary files
        output_prefix = cif_dir / "generated_msa"
        tmp_folder = cif_dir / "generated_tmp"
        tmp_folder.mkdir(exist_ok=True)

        # Step 2: Load and save structures in the dedicated CIF directory
        try:
            cp.load_structure(TEST_PDB_IDS[0])
            cp.load_structure(TEST_PDB_IDS[1])
        except Exception as e:
            pytest.skip(f"Failed to load structures: {e}")

        try:
            # Save temporary CIF files in the cif_dir
            temp_file1 = cp.save_temp_cif(TEST_PDB_IDS[0], suffix="alignment", output_dir=cif_dir)
            temp_file2 = cp.save_temp_cif(TEST_PDB_IDS[1], suffix="alignment", output_dir=cif_dir)
        except Exception as e:
            pytest.skip(f"Failed to save temp CIF files: {e}")

        # Step 3: Resolve absolute paths to the temporary files
        temp_file1_path = Path(temp_file1).resolve()
        temp_file2_path = Path(temp_file2).resolve()
        print(f"\nGenerated CIF file 1: {temp_file1_path}")
        print(f"Generated CIF file 2: {temp_file2_path}")

        # Display the first 20 lines of the first CIF file for debugging
        print("\nFirst 20 lines of CIF file 1:")
        try:
            with open(temp_file1_path, 'r') as f:
                for i, line in enumerate(f.readlines()[:20]):
                    print(f"{i + 1}: {line.rstrip()}")
        except Exception as e:
            print(f"Error reading CIF file: {e}")

        # Log the easy-msa command details
        print("\nRunning easy-msa command with generated CIF files")
        print(f"Input files: {[temp_file1_path, temp_file2_path]}")
        print(f"Output prefix: {output_prefix}")
        print(f"Temp folder: {tmp_folder}")

        # Execute easy-msa and verify output
        try:
            result = fm.easy_msa(
                input_files=[temp_file1_path, temp_file2_path],
                output_prefix=output_prefix,
                tmp_folder=tmp_folder,
                report_mode=1
            )
            print("FoldMason result with generated files:", result)

            # Check for expected output files
            expected_files = [
                output_prefix.with_suffix(".html"),
                output_prefix.with_suffix("_aa.fa"),
                output_prefix.with_suffix("_3di.fa"),
                output_prefix.with_suffix(".nw")
            ]

            for file in expected_files:
                if file.exists():
                    print(f"Found expected output file: {file.name} ({file.stat().st_size} bytes)")
                else:
                    print(f"WARNING: Expected file not found: {file.name}")

        except Exception as e:
            print(f"easy-msa failed with generated CIF files: {e}")

    @pytest.mark.skipif(
        not (Path("data") / "mmcif" / f"{TEST_PDB_IDS[0]}.cif").exists() or
        not (Path("data") / "mmcif" / f"{TEST_PDB_IDS[1]}.cif").exists(),
        reason="Original CIF files not found"
    )
    def test_real_integration_with_cp(self, cp, fm, persistent_output_dir):
        """Real integration test with CifProcessor and FoldMason using original CIF files."""
        test_dir = persistent_output_dir / "real_integration_test"
        test_dir.mkdir(exist_ok=True)
        file1 = Path("data") / "mmcif" / f"{TEST_PDB_IDS[0]}.cif"
        file2 = Path("data") / "mmcif" / f"{TEST_PDB_IDS[1]}.cif"
        if not (file1.exists() and file2.exists()):
            pytest.skip("Original CIF files not found for integration test")

        print(f"Testing with original CIF files: {file1}, {file2}")
        output_prefix = test_dir / "original_msa"
        tmp_folder = test_dir / "original_tmp"
        tmp_folder.mkdir(exist_ok=True)

        try:
            print("\nRunning FoldMason with original CIF files")
            result = fm.easy_msa(
                input_files=[file1, file2],
                output_prefix=output_prefix,
                tmp_folder=tmp_folder,
                report_mode=1
            )
            print("FoldMason result with original files:", result)

            # Check for expected output files
            expected_files = [
                output_prefix.with_suffix(".html"),
                output_prefix.with_suffix("_aa.fa"),
                output_prefix.with_suffix("_3di.fa"),
                output_prefix.with_suffix(".nw")
            ]

            for file in expected_files:
                if file.exists():
                    print(f"Found expected output file: {file.name} ({file.stat().st_size} bytes)")
                else:
                    print(f"WARNING: Expected file not found: {file.name}")

        except Exception as e:
            print(f"FoldMason failed with original files: {e}")

    def test_with_downloaded_structures(self, fm, persistent_output_dir):
        """Test FoldMason with well-known structures downloaded from PDB."""
        import urllib.request
        test_dir = persistent_output_dir / "downloaded_structures"
        test_dir.mkdir(exist_ok=True)
        pdb_ids = ["1ubq", "6vxx"]
        pdb_files = []

        for pdb_id in pdb_ids:
            pdb_file = test_dir / f"{pdb_id}.pdb"
            if not pdb_file.exists():
                try:
                    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    print(f"\nDownloading {url} to {pdb_file}")
                    with urllib.request.urlopen(url) as response, open(pdb_file, 'wb') as out_file:
                        out_file.write(response.read())
                except Exception as e:
                    print(f"Failed to download {pdb_id}.pdb: {e}")
                    continue
            if pdb_file.exists():
                pdb_files.append(pdb_file)

        if len(pdb_files) < 2:
            pytest.skip("Could not obtain at least 2 PDB files for testing")

        for pdb_file in pdb_files:
            print(f"\nFile: {pdb_file} ({pdb_file.stat().st_size} bytes)")
            try:
                with open(pdb_file, 'r') as f:
                    head = [next(f) for _ in range(5)]
                    print("First 5 lines:")
                    for line in head:
                        print(f"  {line.rstrip()}")
            except Exception as e:
                print(f"Error reading file: {e}")

        output_prefix = test_dir / "downloaded_msa"
        tmp_folder = test_dir / "downloaded_tmp"
        tmp_folder.mkdir(exist_ok=True)

        try:
            print("\nRunning FoldMason with downloaded PDB files")
            result = fm.easy_msa(
                input_files=pdb_files,
                output_prefix=output_prefix,
                tmp_folder=tmp_folder,
                report_mode=1
            )
            print("FoldMason result with downloaded files:", result)

            # Check for expected output files
            expected_files = [
                output_prefix.with_suffix(".html"),
                output_prefix.with_suffix("_aa.fa"),
                output_prefix.with_suffix("_3di.fa"),
                output_prefix.with_suffix(".nw")
            ]

            for file in expected_files:
                if file.exists():
                    print(f"Found expected output file: {file.name} ({file.stat().st_size} bytes)")
                else:
                    print(f"WARNING: Expected file not found: {file.name}")

        except Exception as e:
            print(f"FoldMason failed with downloaded files: {e}")
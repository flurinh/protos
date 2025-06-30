"""
Simple test to verify PDB downloads work correctly.
"""

import os
import pytest
from pathlib import Path
from Bio.PDB import PDBList

from protos.loaders.download_structures import download_protein_structures
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.io.paths.path_config import ProtosPaths


def test_simple_pdb_download():
    """Test downloading a single small PDB file."""
    # Use a very small, reliable PDB file
    pdb_id = "1l2y"  # Trp-cage miniprotein - very small
    
    # Use processor's structure directory (managed by ProtosPaths)
    struct_proc = CifBaseProcessor(name="test_download")
    target_dir = struct_proc.path_structure_dir
    
    # Test with PDBList directly first
    pdbl = PDBList()
    print(f"Downloading {pdb_id} to {target_dir}")
    
    # PDBList creates subdirectories - let's see what happens
    file_path = pdbl.retrieve_pdb_file(
        pdb_id, 
        pdir=str(target_dir), 
        file_format="mmCif",
        overwrite=True
    )
    
    print(f"Downloaded to: {file_path}")
    
    # Check what files exist
    all_files = list(target_dir.rglob("*"))
    print(f"Files in target_dir: {all_files}")
    
    # Find the actual downloaded file
    cif_files = list(target_dir.rglob("*.cif"))
    print(f"CIF files found: {cif_files}")
    
    assert len(cif_files) > 0, "Should have downloaded at least one CIF file"
    
    # Now test with our download function
    successful, failed = download_protein_structures(
        [pdb_id],
        target_folder=str(target_dir)
    )
    
    assert len(successful) == 1, f"Should have downloaded {pdb_id}"
    assert len(failed) == 0, f"Should not have any failures"


def test_real_workflow_with_correct_paths():
    """Test complete workflow with correct file paths."""
    # ProtosPaths already configured by conftest
    
    processor = CifBaseProcessor(
        name="test",
        processor_data_dir="structure"
    )
    
    # Download a small protein
    pdb_id = "1l2y"
    target_dir = processor.path_structure_dir
    
    # Use PDBList directly to understand the file structure
    pdbl = PDBList()
    file_path = pdbl.retrieve_pdb_file(
        pdb_id,
        pdir=str(target_dir),
        file_format="mmCif",
        overwrite=True
    )
    
    print(f"File downloaded to: {file_path}")
    
    # PDBList creates a subdirectory structure like mmcif/l2/1l2y.cif
    # We need to move the file to the expected location
    if file_path and os.path.exists(file_path):
        # Move to expected location
        expected_path = target_dir / f"{pdb_id}.cif"
        if not expected_path.exists():
            # Find the actual file
            actual_files = list(target_dir.rglob(f"*{pdb_id}*.cif"))
            if actual_files:
                actual_file = actual_files[0]
                # Copy to expected location
                import shutil
                shutil.copy2(actual_file, expected_path)
                print(f"Copied {actual_file} to {expected_path}")
    
    # Now load with processor
    processor.load_structure(pdb_id)
    
    # Verify
    assert processor.data is not None
    assert len(processor.data) > 0
    assert processor.data['pdb_id'].iloc[0] == pdb_id
    
    # No cleanup needed - conftest handles it
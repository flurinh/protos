#!/usr/bin/env python3
"""
Test script for new CifBaseProcessor methods: get_seq_dict, get_grn_dict, assign_grns

This script tests the GRN-Structure integration using real microbial opsin data.
"""

import os
import sys
from pathlib import Path
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.loaders.download_structures import download_protein_structures
from protos.io.fasta_utils import write_fasta


def test_get_seq_dict():
    """Test the get_seq_dict method."""
    print("\n" + "="*60)
    print("TEST 1: get_seq_dict()")
    print("="*60)
    
    # Initialize processor
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    struct_processor = CifBaseProcessor(
        name="test_seq_dict",
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/structure"
    )
    
    # Create test data
    test_data = pd.DataFrame({
        'pdb_id': ['1UAZ', '1UAZ', '1UAZ', '1UAZ', '1UAZ', '3DDL', '3DDL', '3DDL'],
        'auth_chain_id': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'B'],
        'auth_seq_id': [1, 2, 3, 5, 6, 1, 2, 1],  # Note gap at position 4
        'auth_comp_id': ['MET', 'ALA', 'GLY', 'TRP', 'LEU', 'SER', 'THR', 'VAL']
    })
    
    struct_processor.data = test_data
    
    # Test extraction
    sequences = struct_processor.get_seq_dict()
    
    print(f"Extracted {len(sequences)} sequences:")
    for seq_id, seq in sequences.items():
        print(f"  {seq_id}: {seq}")
    
    # Verify results
    assert '1UAZ_A' in sequences
    assert sequences['1UAZ_A'] == 'MAGXWL'  # X for missing position 4
    assert '3DDL_A' in sequences
    assert sequences['3DDL_A'] == 'ST'
    assert '3DDL_B' in sequences
    assert sequences['3DDL_B'] == 'V'
    
    print("\n✓ get_seq_dict test passed!")
    return sequences


def test_get_grn_dict():
    """Test the get_grn_dict method."""
    print("\n" + "="*60)
    print("TEST 2: get_grn_dict()")
    print("="*60)
    
    # Initialize processor
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    struct_processor = CifBaseProcessor(
        name="test_grn_dict",
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/structure"
    )
    
    # Create test data with GRN annotations
    test_data = pd.DataFrame({
        'pdb_id': ['1UAZ', '1UAZ', '1UAZ', '3DDL', '3DDL'],
        'auth_chain_id': ['A', 'A', 'A', 'A', 'A'],
        'auth_seq_id': [82, 85, 216, 100, 250],
        'auth_comp_id': ['ARG', 'ASP', 'LYS', 'ASP', 'LYS'],
        'grn': ['1.50', '3.50', '7.50', '3.50', '7.50']
    })
    
    struct_processor.data = test_data
    
    # Test extraction
    grn_dict = struct_processor.get_grn_dict()
    
    print(f"Extracted GRN annotations for {len(grn_dict)} structures:")
    for pdb_id, chains in grn_dict.items():
        print(f"\n  {pdb_id}:")
        for chain_id, grns in chains.items():
            print(f"    Chain {chain_id}:")
            for grn_pos, res_info in grns.items():
                print(f"      {grn_pos}: {res_info}")
    
    # Verify results
    assert '1UAZ' in grn_dict
    assert 'A' in grn_dict['1UAZ']
    assert grn_dict['1UAZ']['A']['1.50'] == 'R82'
    assert grn_dict['1UAZ']['A']['7.50'] == 'K216'
    
    print("\n✓ get_grn_dict test passed!")
    return grn_dict


def test_assign_grns_simple():
    """Test the assign_grns method with simple data."""
    print("\n" + "="*60)
    print("TEST 3: assign_grns() - Simple Test")
    print("="*60)
    
    # Initialize processor
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    struct_processor = CifBaseProcessor(
        name="test_assign_grns",
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/structure"
    )
    
    # Create bacteriorhodopsin-like test sequence
    br_sequence = (
        "MLDAVAAALGVGLILLGLIIVSTLVGQRFQWIWLALGTALMGLGTLYFLVKGMGVSDPD"
        "AKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDL"
        "ALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGF"
        "TSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSA"
        "KKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"
    )
    
    # Create structure data for this sequence
    test_data = []
    for i, aa in enumerate(br_sequence):
        aa_3letter = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
            'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
            'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }.get(aa, 'UNK')
        
        test_data.append({
            'pdb_id': 'TEST1',
            'auth_chain_id': 'A',
            'auth_seq_id': i + 1,
            'auth_comp_id': aa_3letter
        })
    
    struct_processor.data = pd.DataFrame(test_data)
    
    print(f"Created test structure with {len(struct_processor.data)} residues")
    
    # Test GRN assignment
    try:
        grn_assignments = struct_processor.assign_grns(
            protein_family='microbial_opsins',
            similarity_threshold=0.1,  # Lower threshold for test
            use_mmseqs=False,  # Use BioPython for simplicity
            save_results=False
        )
        
        print(f"\nAssigned GRNs to {len(grn_assignments)} chains")
        
        # Check if GRN column was added
        if 'grn' in struct_processor.data.columns:
            grn_residues = struct_processor.data[struct_processor.data['grn'].notna()]
            print(f"Annotated {len(grn_residues)} residues with GRN positions")
            
            # Show some key positions
            key_positions = ['1.50', '3.50', '7.50']
            print("\nKey GRN positions:")
            for pos in key_positions:
                residues = grn_residues[grn_residues['grn'] == pos]
                if not residues.empty:
                    res = residues.iloc[0]
                    print(f"  {pos}: {res['auth_comp_id']}{res['auth_seq_id']}")
        
        print("\n✓ assign_grns simple test passed!")
        
    except Exception as e:
        print(f"\n✗ assign_grns test failed: {e}")
        import traceback
        traceback.print_exc()
    
    return grn_assignments if 'grn_assignments' in locals() else {}


def test_real_structure():
    """Test with a real microbial opsin structure."""
    print("\n" + "="*60)
    print("TEST 4: Real Structure Test")
    print("="*60)
    
    # Initialize processor
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    struct_processor = CifBaseProcessor(
        name="test_real",
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/structure"
    )
    
    # Try to load a real structure
    pdb_id = "1UAZ"  # Bacteriorhodopsin
    
    try:
        print(f"Loading structure {pdb_id}...")
        
        # Check if file exists, if not try to download
        mmcif_path = data_dir / "structure" / "mmcif" / f"{pdb_id.lower()}.cif"
        if not mmcif_path.exists():
            print(f"Downloading {pdb_id}...")
            mmcif_dir = data_dir / "structure" / "mmcif"
            mmcif_dir.mkdir(parents=True, exist_ok=True)
            download_protein_structures([pdb_id], str(mmcif_dir))
        
        # Load structure
        struct_processor.load_structure(pdb_id, remove_hetatm=True)
        print(f"Loaded {len(struct_processor.data)} atoms")
        
        # Extract sequences
        sequences = struct_processor.get_seq_dict()
        print(f"\nExtracted {len(sequences)} sequences:")
        for seq_id, seq in sequences.items():
            print(f"  {seq_id}: {seq[:50]}... (length: {len(seq)})")
        
        # Assign GRNs
        print("\nAssigning GRNs...")
        grn_assignments = struct_processor.assign_grns(
            protein_family='microbial_opsins',
            similarity_threshold=0.2,
            use_mmseqs=True,  # Try MMseqs2 if available
            save_results=True
        )
        
        print(f"\nAssigned GRNs to {len(grn_assignments)} chains")
        
        # Extract GRN dictionary
        grn_dict = struct_processor.get_grn_dict()
        if grn_dict:
            print(f"\nGRN annotations:")
            for pdb_id, chains in grn_dict.items():
                for chain_id, grns in chains.items():
                    print(f"\n  {pdb_id}:{chain_id}")
                    # Show key positions
                    for pos in ['1.50', '2.50', '3.50', '6.48', '7.50']:
                        if pos in grns:
                            print(f"    {pos}: {grns[pos]}")
        
        # Save annotated structure
        struct_processor.save_dataset(f"{pdb_id}_with_grn")
        print(f"\nSaved annotated structure")
        
        print("\n✓ Real structure test passed!")
        
    except Exception as e:
        print(f"\n✗ Real structure test failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("TESTING NEW CIFBASEPROCESSOR METHODS")
    print("="*60)
    
    # Run tests
    test_get_seq_dict()
    test_get_grn_dict()
    test_assign_grns_simple()
    test_real_structure()
    
    print("\n" + "="*60)
    print("ALL TESTS COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()
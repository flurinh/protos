#!/usr/bin/env python
"""
Test script for StructureProcessor.assign_grns() method.

This demonstrates how to use the GRN assignment functionality with structure data.
"""

import os
from pathlib import Path
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor

def test_assign_grns_microbial_opsins():
    """Test GRN assignment for microbial opsin structures."""
    
    # Setup paths
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    print("Initializing StructureProcessor...")
    struct_proc = StructureProcessor(name="test_grn_assignment", paths=paths)
    
    # Load microbial opsin dataset
    print("\nLoading microbial opsin structures...")
    struct_proc.load_dataset('test_mo')  # or load individual structures
    
    # Display loaded structures
    print(f"\nLoaded {len(struct_proc.pdb_ids)} structures:")
    for pdb_id in struct_proc.pdb_ids:
        chains = struct_proc.get_chains(pdb_id)
        print(f"  {pdb_id}: chains {', '.join(chains)}")
    
    # Assign GRNs
    print("\nAssigning GRN positions...")
    grn_assignments = struct_proc.assign_grns(
        protein_family='microbial_opsins',
        similarity_threshold=0.2,  # 20% identity threshold
        use_mmseqs=True,           # Use fast MMseqs2 alignment
    )
    
    # Check results
    print(f"\nGRN assignment complete!")
    print(f"Assigned GRNs to {len(grn_assignments)} chains")
    
    # Check if GRN column was added to structure data
    if 'grn' in struct_proc.data.columns:
        grn_residues = struct_proc.data[struct_proc.data['grn'].notna()]
        print(f"Total residues with GRN annotations: {len(grn_residues)}")
        
        # Show key positions for each chain
        for chain_id in grn_assignments:
            print(f"\n{chain_id}:")
            pdb_id, chain = chain_id.rsplit('_', 1)
            
            # Get key positions
            key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
            chain_data = struct_proc.data[
                (struct_proc.data['pdb_id'] == pdb_id) & 
                (struct_proc.data['auth_chain_id'] == chain)
            ]
            
            for pos in key_positions:
                residues = chain_data[chain_data['grn'] == pos]
                if not residues.empty:
                    res = residues.iloc[0]
                    print(f"  {pos}: {res['res_name1l']}{res['auth_seq_id']}")
    
    # Get GRN dictionary
    print("\nTesting get_grn_dict()...")
    grn_dict = struct_proc.get_grn_dict()
    print(f"GRN dictionary contains {len(grn_dict)} entries")
    
    return grn_assignments


def test_assign_grns_gpcr():
    """Test GRN assignment for GPCR structures."""
    
    # Setup paths
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    print("Initializing StructureProcessor...")
    struct_proc = StructureProcessor(name="test_gpcr_grn", paths=paths)
    
    # Load GPCR dataset
    print("\nLoading GPCR structures...")
    struct_proc.load_dataset('gpcr_agonist_inverse_agonist')
    
    # Limit to a few structures for testing
    struct_proc.pdb_ids = struct_proc.pdb_ids[:5]
    print(f"\nProcessing first {len(struct_proc.pdb_ids)} structures")
    
    # Assign GRNs
    print("\nAssigning GRN positions...")
    grn_assignments = struct_proc.assign_grns(
        protein_family='gpcr_a',
        similarity_threshold=0.25,  # 25% identity threshold for GPCRs
        use_mmseqs=True
    )
    
    print(f"\nGRN assignment complete!")
    print(f"Assigned GRNs to {len(grn_assignments)} chains")
    
    return grn_assignments


def test_assign_grns_custom():
    """Test GRN assignment with custom parameters."""
    
    # Setup paths
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    print("Initializing StructureProcessor...")
    struct_proc = StructureProcessor(name="test_custom_grn", paths=paths)
    
    # Load a single structure
    print("\nLoading single structure...")
    struct_proc.load_structure('1UAZ')  # Sensory rhodopsin II
    
    # Assign GRNs with custom reference table
    print("\nAssigning GRN positions with custom settings...")
    grn_assignments = struct_proc.assign_grns(
        protein_family='microbial_opsins',
        similarity_threshold=0.2,     # Lower threshold
        use_mmseqs=False,              # Use BioPython for single structure
        reference_table='mo_ref'       # Explicit reference table
    )
    
    print(f"\nGRN assignment complete!")
    
    # Show detailed results
    for chain_id, grn_mapping in grn_assignments.items():
        print(f"\n{chain_id}:")
        print(f"  Total GRN positions assigned: {len(grn_mapping)}")
        
        # Show some example mappings
        example_positions = list(grn_mapping.items())[:10]
        print("  Example mappings:")
        for res, grn in example_positions:
            print(f"    {res} -> {grn}")
    
    return grn_assignments


if __name__ == "__main__":
    print("=" * 60)
    print("Testing StructureProcessor.assign_grns()")
    print("=" * 60)
    
    # Test 1: Microbial opsins
    try:
        print("\nTest 1: Microbial Opsins")
        test_assign_grns_microbial_opsins()
    except Exception as e:
        print(f"Test 1 failed: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 2: GPCRs
    try:
        print("\n" + "=" * 60)
        print("Test 2: GPCRs")
        test_assign_grns_gpcr()
    except Exception as e:
        print(f"Test 2 failed: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 3: Custom settings
    try:
        print("\n" + "=" * 60)
        print("Test 3: Custom Settings")
        test_assign_grns_custom()
    except Exception as e:
        print(f"Test 3 failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 60)
    print("All tests completed!")
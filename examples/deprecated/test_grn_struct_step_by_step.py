#!/usr/bin/env python3
"""
Step-by-step test of GRN-Structure integration functionality.

This script tests:
1. Download structure dataset
2. Test get_seq_dict function
3. Test GRN assignment to sequences
4. Test mapping GRNs back to structure data
"""

import os
import sys
from pathlib import Path
import pandas as pd
import json

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.io.ingest.download_structures import download_protein_structures
from protos.io.fasta_utils import write_fasta


def step1_download_structures():
    """Step 1: Download microbial opsin structures."""
    print("\n" + "="*80)
    print("STEP 1: Download Structure Dataset")
    print("="*80)
    
    # Define data directory
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    mmcif_dir = data_dir / "structure" / "mmcif"
    mmcif_dir.mkdir(parents=True, exist_ok=True)
    
    # Select a few microbial opsin structures
    mo_pdb_ids = ["1UAZ", "3DDL", "4PXK"]  # Bacteriorhodopsin and others
    
    print(f"Target directory: {mmcif_dir}")
    print(f"PDB IDs to download: {mo_pdb_ids}")
    
    # Check what's already downloaded
    existing_files = []
    for pdb_id in mo_pdb_ids:
        if (mmcif_dir / f"{pdb_id.lower()}.cif").exists():
            existing_files.append(pdb_id)
    
    if existing_files:
        print(f"\nAlready downloaded: {existing_files}")
    
    # Download missing structures
    to_download = [pid for pid in mo_pdb_ids if pid not in existing_files]
    if to_download:
        print(f"\nDownloading: {to_download}")
        try:
            download_protein_structures(to_download, str(mmcif_dir))
            print("✓ Download complete")
        except Exception as e:
            print(f"✗ Download failed: {e}")
            return False
    else:
        print("\n✓ All structures already downloaded")
    
    # Verify downloads
    downloaded_count = 0
    for pdb_id in mo_pdb_ids:
        if (mmcif_dir / f"{pdb_id.lower()}.cif").exists():
            downloaded_count += 1
    
    print(f"\nStructures available: {downloaded_count}/{len(mo_pdb_ids)}")
    
    # Create dataset definition
    dataset_def = {
        "id": "test_mo",
        "name": "Test Microbial Opsins",
        "description": "Test dataset for GRN-Structure integration",
        "type": "structure",
        "pdb_ids": mo_pdb_ids
    }
    
    dataset_path = data_dir / "structure" / "structure_dataset" / "test_mo.json"
    dataset_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(dataset_path, 'w') as f:
        json.dump(dataset_def, f, indent=2)
    
    print(f"\n✓ Created dataset definition: {dataset_path}")
    
    return True


def step2_test_get_seq_dict():
    """Step 2: Test get_seq_dict function."""
    print("\n" + "="*80)
    print("STEP 2: Test get_seq_dict Function")
    print("="*80)
    
    # Initialize processor with correct paths
    data_dir = Path(__file__).parent.parent / "src" / "protos" / "reference_data"
    struct_processor = CifBaseProcessor(
        name="test_seq_extraction",
        data_root=str(data_dir),
        processor_data_dir="structure"
    )
    
    # Load real structure files
    print("\nLoading real microbial opsin structures...")
    try:
        # Load individual structures
        loaded_count = 0
        for pdb_id in ["1UAZ", "3DDL", "4PXK"]:
            try:
                struct_processor.load_structure(pdb_id.lower())
                loaded_count += 1
                print(f"  ✓ Loaded {pdb_id}")
            except Exception as e:
                print(f"  ✗ Failed to load {pdb_id}: {e}")
        
        if loaded_count == 0:
            raise Exception("No structures loaded")
            
        print(f"\n✓ Loaded {loaded_count} structures with {len(struct_processor.data)} total records")
        
        # Show data format
        print("\nData format (first 5 rows):")
        print(struct_processor.data.head())
        
        print("\nColumns in data:")
        print(list(struct_processor.data.columns))
        
        # Show unique PDB IDs and chains
        unique_pdbs = struct_processor.data['pdb_id'].unique()
        print(f"\nUnique PDB IDs: {unique_pdbs}")
        
        for pdb_id in unique_pdbs:
            pdb_data = struct_processor.data[struct_processor.data['pdb_id'] == pdb_id]
            chains = pdb_data['auth_chain_id'].unique()
            print(f"  {pdb_id}: chains {chains}")
        
    except Exception as e:
        print(f"✗ Failed to load structures: {e}")
        # Create minimal test data
        print("\nCreating minimal test data...")
        test_data = pd.DataFrame({
            'pdb_id': ['TEST1'] * 10 + ['TEST2'] * 5,
            'auth_chain_id': ['A'] * 10 + ['B'] * 5,
            'auth_seq_id': list(range(1, 11)) + list(range(1, 6)),
            'auth_comp_id': ['MET', 'ALA', 'GLY', 'VAL', 'LEU', 
                            'ILE', 'SER', 'THR', 'TYR', 'PRO',
                            'ARG', 'ASP', 'GLU', 'LYS', 'HIS']
        })
        struct_processor.data = test_data
        print(f"✓ Created test data with {len(struct_processor.data)} records")
    
    # Extract sequences
    print("\nExtracting sequences...")
    sequences = struct_processor.get_seq_dict()
    
    print(f"\n✓ Extracted {len(sequences)} sequences:")
    for seq_id, seq in sequences.items():
        print(f"  {seq_id}: {seq[:60]}{'...' if len(seq) > 60 else ''} (length: {len(seq)})")
    
    # Save sequences
    fasta_path = data_dir / "sequence" / "fasta" / "test_mo_extracted.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(sequences, str(fasta_path))
    print(f"\n✓ Saved sequences to: {fasta_path}")
    
    return struct_processor, sequences


def step3_test_assign_grns(struct_processor, sequences):
    """Step 3: Test GRN assignment to sequences."""
    print("\n" + "="*80)
    print("STEP 3: Test GRN Assignment")
    print("="*80)
    
    if not sequences:
        print("✗ No sequences available for GRN assignment")
        return None
    
    print(f"\nAssigning GRNs to {len(sequences)} sequences...")
    print("Using protein family: microbial_opsins")
    print("Similarity threshold: 20%")
    
    try:
        # Assign GRNs
        grn_assignments = struct_processor.assign_grns(
            protein_family='microbial_opsins',
            similarity_threshold=0.2,
            use_mmseqs=True,  # Try MMseqs2 first
            save_results=True
        )
        
        if grn_assignments:
            print(f"\n✓ Successfully assigned GRNs to {len(grn_assignments)} chains")
            
            # Show some key positions for each chain
            for chain_id, grn_series in grn_assignments.items():
                print(f"\n{chain_id}:")
                key_positions = ['1.50', '2.50', '3.50', '6.48', '7.50']
                found_positions = []
                
                for pos in key_positions:
                    if pos in grn_series and grn_series[pos] != '-':
                        found_positions.append(f"{pos}: {grn_series[pos]}")
                
                if found_positions:
                    print("  Key positions: " + ", ".join(found_positions))
                else:
                    print("  No key positions found")
        else:
            print("\n✗ No GRN assignments made")
            
    except Exception as e:
        print(f"\n✗ GRN assignment failed: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    return grn_assignments


def step4_test_grn_mapping(struct_processor, grn_assignments):
    """Step 4: Test mapping GRNs back to structure data."""
    print("\n" + "="*80)
    print("STEP 4: Test GRN Mapping to Structure")
    print("="*80)
    
    if not grn_assignments:
        print("✗ No GRN assignments available for mapping")
        return
    
    # Check if GRN column exists
    if 'grn' not in struct_processor.data.columns:
        print("✗ GRN column not found in structure data")
        return
    
    # Get GRN dictionary
    print("\nExtracting GRN annotations from structure...")
    grn_dict = struct_processor.get_grn_dict()
    
    if grn_dict:
        print(f"\n✓ Extracted GRN annotations for {len(grn_dict)} structures:")
        
        for pdb_id, chains in grn_dict.items():
            print(f"\n{pdb_id}:")
            for chain_id, grns in chains.items():
                print(f"  Chain {chain_id}: {len(grns)} GRN positions")
                
                # Show some key positions
                key_positions = ['1.50', '2.50', '3.50', '6.48', '7.50']
                for pos in key_positions:
                    if pos in grns:
                        print(f"    {pos}: {grns[pos]}")
    else:
        print("\n✗ No GRN annotations found in structure data")
    
    # Verify mapping accuracy
    print("\nVerifying GRN mapping accuracy...")
    
    # Count annotated residues
    grn_residues = struct_processor.data[struct_processor.data['grn'].notna()]
    print(f"\nTotal residues with GRN annotations: {len(grn_residues)}")
    
    # Show distribution by PDB and chain
    grn_summary = grn_residues.groupby(['pdb_id', 'auth_chain_id'])['grn'].count()
    print("\nGRN annotations by structure and chain:")
    for (pdb_id, chain_id), count in grn_summary.items():
        print(f"  {pdb_id}:{chain_id} - {count} residues")
    
    # Show unique GRN positions
    unique_grns = grn_residues['grn'].unique()
    print(f"\nUnique GRN positions annotated: {len(unique_grns)}")
    
    # Check specific example
    if len(grn_residues) > 0:
        print("\nExample annotations (first 10):")
        for idx, row in grn_residues.head(10).iterrows():
            res_col = 'auth_comp_id' if 'auth_comp_id' in row.index else 'res_name3l'
            print(f"  {row['pdb_id']}:{row['auth_chain_id']}:{row[res_col]}{row['auth_seq_id']} -> {row['grn']}")
    
    # Save annotated structure
    output_path = Path(__file__).parent.parent / "src" / "protos" / "reference_data" / "structure" / "datasets"
    output_path.mkdir(parents=True, exist_ok=True)
    
    struct_processor.save_data("test_mo_with_grn", struct_processor.data)
    print(f"\n✓ Saved annotated structure dataset")


def main():
    """Run all steps."""
    print("\n" + "="*80)
    print("GRN-STRUCTURE INTEGRATION TEST")
    print("Step-by-step functionality test")
    print("="*80)
    
    # Step 1: Download structures
    success = step1_download_structures()
    if not success:
        print("\nStep 1 failed. Continuing with available data...")
    
    # Step 2: Test sequence extraction
    struct_processor, sequences = step2_test_get_seq_dict()
    
    # Step 3: Test GRN assignment
    grn_assignments = step3_test_assign_grns(struct_processor, sequences)
    
    # Step 4: Test GRN mapping
    step4_test_grn_mapping(struct_processor, grn_assignments)
    
    print("\n" + "="*80)
    print("TEST COMPLETE")
    print("="*80)
    
    # Summary
    print("\nSummary:")
    print(f"  - Sequences extracted: {len(sequences) if 'sequences' in locals() else 0}")
    print(f"  - GRN assignments: {len(grn_assignments) if grn_assignments else 0}")
    
    if struct_processor and struct_processor.data is not None and 'grn' in struct_processor.data.columns:
        grn_count = struct_processor.data['grn'].notna().sum()
        print(f"  - Residues with GRN: {grn_count}")


if __name__ == "__main__":
    main()
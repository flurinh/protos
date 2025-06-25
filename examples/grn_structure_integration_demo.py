#!/usr/bin/env python3
"""
GRN-Structure Integration Demo
==============================

This script demonstrates the integration between Structure and GRN processors:
1. Load microbial opsin structures
2. Extract sequences from structure data
3. Assign GRN numbers to sequences
4. Map GRN annotations back to structure residues

Usage:
    python examples/grn_structure_integration_demo.py
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
# GRNAssignment class doesn't exist - using functions directly
from protos.processing.grn.grn_assignment import get_correctly_aligned_grns
from protos.processing.grn.grn_table_utils import init_row_from_alignment, expand_annotation
from protos.processing.grn.grn_utils import get_seq
from protos.processing.sequence.seq_alignment import init_aligner, align_blosum62, format_alignment
from protos.io.fasta_utils import write_fasta, read_fasta
from protos.loaders.download_structures import download_protein_structures


def extract_sequences_from_structure(struct_df):
    """
    Extract protein sequences from structure dataframe.
    
    Args:
        struct_df: DataFrame with structure data (must have pdb_id, auth_chain_id, auth_seq_id, auth_comp_id)
    
    Returns:
        Dictionary of {pdb_chain: sequence}
    """
    sequences = {}
    
    # Group by PDB ID and chain
    for (pdb_id, chain_id), group in struct_df.groupby(['pdb_id', 'auth_chain_id']):
        # Sort by sequence ID
        group = group.sort_values('auth_seq_id')
        
        # Convert 3-letter codes to 1-letter
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            'MSE': 'M',  # Selenomethionine
            'UNK': 'X'   # Unknown
        }
        
        sequence = []
        prev_seq_id = None
        
        for _, residue in group.iterrows():
            # Skip if not a standard amino acid
            comp_id = residue['auth_comp_id'].upper()
            if comp_id not in aa_map:
                continue
                
            # Check for gaps in sequence
            seq_id = residue['auth_seq_id']
            if prev_seq_id is not None and seq_id > prev_seq_id + 1:
                # Add X for missing residues
                for _ in range(prev_seq_id + 1, seq_id):
                    sequence.append('X')
            
            sequence.append(aa_map[comp_id])
            prev_seq_id = seq_id
        
        if sequence:
            seq_key = f"{pdb_id}_{chain_id}"
            sequences[seq_key] = ''.join(sequence)
    
    return sequences


def add_grn_to_structure(struct_df, grn_assignments):
    """
    Add GRN annotations to structure dataframe.
    
    Args:
        struct_df: Structure dataframe
        grn_assignments: Dictionary of {pdb_chain: {grn_pos: residue_info}}
    
    Returns:
        Structure dataframe with added 'grn' column
    """
    # Add GRN column if not exists
    if 'grn' not in struct_df.columns:
        struct_df['grn'] = None
    
    # Process each chain's GRN assignments
    for chain_key, grn_row in grn_assignments.items():
        pdb_id, chain_id = chain_key.split('_')
        
        for grn_pos, res_info in grn_row.items():
            if res_info and res_info != '-':
                try:
                    # Extract residue number from format like 'K296'
                    res_num = int(res_info[1:])
                    res_type = res_info[0]
                    
                    # Find matching residues in structure
                    mask = (
                        (struct_df['pdb_id'] == pdb_id) &
                        (struct_df['auth_chain_id'] == chain_id) &
                        (struct_df['auth_seq_id'] == res_num)
                    )
                    
                    # Verify residue type matches
                    if mask.any():
                        struct_res = struct_df.loc[mask, 'auth_comp_id'].iloc[0]
                        expected_res = {
                            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
                            'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
                            'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
                            'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
                        }.get(res_type, 'UNK')
                        
                        if struct_res.upper() == expected_res:
                            struct_df.loc[mask, 'grn'] = grn_pos
                        else:
                            print(f"Warning: Residue type mismatch at {pdb_id}:{chain_id}:{res_num} - "
                                  f"expected {expected_res}, found {struct_res}")
                    
                except (ValueError, IndexError) as e:
                    print(f"Error parsing GRN assignment {res_info}: {e}")
                    continue
    
    return struct_df


def main():
    """Run the GRN-Structure integration demo."""
    print("\n" + "="*60)
    print("GRN-STRUCTURE INTEGRATION DEMO")
    print("="*60)
    
    # Set up paths
    project_root = Path(__file__).parent.parent
    data_dir = project_root / "src" / "protos" / "reference_data"
    
    # Initialize processors
    print("\n1. Initializing processors...")
    struct_processor = CifBaseProcessor(
        name="demo_struct",
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/structure"
    )
    
    grn_processor = GRNBaseProcessor(
        name="demo_grn", 
        data_root=str(data_dir.parent.parent),
        processor_data_dir="reference_data/grn"
    )
    
    # Step 1: Load structure dataset
    print("\n2. Loading microbial opsin structures...")
    try:
        # Try to load existing dataset
        struct_processor.load_dataset("microbial_opsins", apply_dtypes=True)
        print(f"Loaded {len(struct_processor.data)} structure records")
    except Exception as e:
        print(f"Could not load dataset: {e}")
        print("Please ensure microbial opsin structures are downloaded")
        
        # Optionally download a few structures for demo
        demo_pdbs = ["1UAZ", "3DDL", "4PXK"]  # Bacteriorhodopsin and others
        print(f"\nDownloading demo structures: {demo_pdbs}")
        
        mmcif_dir = data_dir / "structure" / "mmcif"
        mmcif_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            download_protein_structures(demo_pdbs, str(mmcif_dir))
            
            # Load individual structures
            all_data = []
            for pdb_id in demo_pdbs:
                try:
                    df = struct_processor.load_structure(pdb_id, remove_hetatm=True)
                    if df is not None:
                        all_data.append(df)
                except Exception as e:
                    print(f"Failed to load {pdb_id}: {e}")
            
            if all_data:
                struct_processor.data = pd.concat(all_data, ignore_index=True)
                print(f"Loaded {len(struct_processor.data)} structure records")
            else:
                print("No structures loaded. Exiting.")
                return
        except Exception as e:
            print(f"Failed to download structures: {e}")
            return
    
    # Step 2: Extract sequences
    print("\n3. Extracting sequences from structures...")
    sequences = extract_sequences_from_structure(struct_processor.data)
    print(f"Extracted {len(sequences)} sequences")
    
    # Show sample sequences
    for seq_id, seq in list(sequences.items())[:3]:
        print(f"  {seq_id}: {seq[:60]}... (length: {len(seq)})")
    
    # Save sequences
    fasta_path = data_dir / "sequence" / "fasta" / "mo_structures.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(sequences, str(fasta_path))
    print(f"Saved sequences to: {fasta_path}")
    
    # Step 3: Load GRN reference table
    print("\n4. Loading GRN reference table...")
    try:
        grn_processor.load_grn_table("mo_ref")
        print(f"Loaded reference table with {len(grn_processor.grn_table)} proteins")
    except Exception as e:
        print(f"Failed to load GRN reference: {e}")
        return
    
    # Step 4: Assign GRNs to sequences
    print("\n5. Assigning GRNs to sequences...")
    
    # Initialize aligner
    aligner = init_aligner()
    grn_assignments = {}
    
    for seq_id, sequence in sequences.items():
        print(f"\nProcessing {seq_id}...")
        
        # Find best reference match
        best_ref = None
        best_score = -float('inf')
        best_alignment = None
        
        # Try a few key references
        ref_candidates = ['BR', '1UAZ', '3DDL']  # Bacteriorhodopsin variants
        
        for ref_id in ref_candidates:
            if ref_id in grn_processor.grn_table.index:
                ref_seq = get_seq(ref_id, grn_processor.grn_table)
                if ref_seq:
                    try:
                        alignment = align_blosum62(sequence, ref_seq, aligner)
                        if alignment.score > best_score:
                            best_score = alignment.score
                            best_ref = ref_id
                            best_alignment = alignment
                    except Exception as e:
                        print(f"  Failed to align with {ref_id}: {e}")
        
        if best_ref and best_alignment:
            print(f"  Best match: {best_ref} (score: {best_score:.1f})")
            
            # Format alignment and transfer GRNs
            formatted = format_alignment(best_alignment)
            
            # Get reference GRN mapping
            ref_row = grn_processor.grn_table.loc[best_ref]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
            seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
            
            # Initialize GRN row from alignment
            new_row = init_row_from_alignment(formatted, seq_pos2grn)
            
            # Try to expand annotation
            try:
                grn_list, rn_list, missing = expand_annotation(
                    new_row,
                    sequence,
                    formatted,
                    max_alignment_gap=1,
                    protein_family='microbial_opsins',
                    verbose=0
                )
                
                # Create final GRN row
                final_row = pd.Series(dict(zip(grn_list, rn_list)))
                
                # Add missing columns
                for col in grn_processor.grn_table.columns:
                    if col not in final_row.index:
                        final_row[col] = '-'
                
                # Reorder to match reference
                final_row = final_row[grn_processor.grn_table.columns]
                
                grn_assignments[seq_id] = final_row
                
                # Show key positions
                key_positions = ['1.50', '2.50', '3.50', '6.48', '7.50']
                print("  Key GRN positions:")
                for pos in key_positions:
                    if pos in final_row and final_row[pos] != '-':
                        print(f"    {pos}: {final_row[pos]}")
                        
            except Exception as e:
                print(f"  Failed to expand annotation: {e}")
                grn_assignments[seq_id] = new_row
        else:
            print(f"  No suitable reference found")
    
    # Step 5: Add GRNs to structure data
    print("\n6. Adding GRN annotations to structure data...")
    
    # Convert Series to dict format for add_grn_to_structure
    grn_dict = {}
    for seq_id, grn_series in grn_assignments.items():
        grn_dict[seq_id] = grn_series.to_dict()
    
    annotated_data = add_grn_to_structure(struct_processor.data, grn_dict)
    
    # Count annotated residues
    grn_residues = annotated_data[annotated_data['grn'].notna()]
    print(f"\nAnnotated {len(grn_residues)} residues with GRN positions")
    
    # Show statistics by PDB
    print("\nGRN annotations by structure:")
    for pdb_id in annotated_data['pdb_id'].unique():
        pdb_grns = annotated_data[
            (annotated_data['pdb_id'] == pdb_id) & 
            (annotated_data['grn'].notna())
        ]
        if len(pdb_grns) > 0:
            unique_grns = pdb_grns['grn'].nunique()
            print(f"  {pdb_id}: {len(pdb_grns)} residues, {unique_grns} unique GRN positions")
    
    # Step 6: Save results
    print("\n7. Saving results...")
    
    # Save GRN table
    grn_table = pd.DataFrame(grn_assignments).T
    grn_output = data_dir / "grn" / "tables" / "mo_structures_grn.csv"
    grn_output.parent.mkdir(parents=True, exist_ok=True)
    grn_table.to_csv(grn_output)
    print(f"Saved GRN table to: {grn_output}")
    
    # Save annotated structure
    struct_processor.data = annotated_data
    struct_processor.save_dataset("mo_with_grn")
    print(f"Saved annotated structure dataset")
    
    # Step 7: Demonstrate usage - find all lysines at position 7.50
    print("\n8. Example usage - Finding Schiff base lysines (7.50):")
    schiff_base = annotated_data[annotated_data['grn'] == '7.50']
    
    if not schiff_base.empty:
        for (pdb_id, chain_id), group in schiff_base.groupby(['pdb_id', 'auth_chain_id']):
            residue = group.iloc[0]
            print(f"  {pdb_id}:{chain_id} - {residue['auth_comp_id']}{residue['auth_seq_id']} at position 7.50")
            
            # Calculate distance to retinal if coordinates available
            if all(coord in residue for coord in ['x', 'y', 'z']):
                coords = np.array([residue['x'], residue['y'], residue['z']])
                print(f"    Coordinates: ({coords[0]:.1f}, {coords[1]:.1f}, {coords[2]:.1f})")
    
    print("\n" + "="*60)
    print("DEMO COMPLETE")
    print("="*60)
    print("\nSummary:")
    print(f"  - Structures processed: {len(annotated_data['pdb_id'].unique())}")
    print(f"  - Sequences extracted: {len(sequences)}")
    print(f"  - GRN assignments: {len(grn_assignments)}")
    print(f"  - Residues annotated: {len(grn_residues)}")
    print(f"  - Unique GRN positions: {grn_residues['grn'].nunique()}")


if __name__ == "__main__":
    main()
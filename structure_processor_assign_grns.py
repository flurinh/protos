#!/usr/bin/env python
"""
Revised implementation of assign_grns method for StructureProcessor.

This follows the exact pattern from test_grn_assignment.py and annotate_microbial_opsins_structures.py
but stores results in the 'grn' column of structure data.
"""

def assign_grns(self, protein_family='microbial_opsins', reference_table=None, 
                similarity_threshold=0.2, use_mmseqs=True, verbose=True):
    """
    Assign GRN annotations to structure data.
    
    This method follows the same process as test_grn_assignment.py:
    1. Extracts sequences from the loaded structures
    2. Loads the specified reference GRN table
    3. Aligns sequences against reference
    4. Assigns GRN positions using init_row_from_alignment and expand_annotation
    5. Stores results in 'grn' column of self.data
    
    Args:
        protein_family: Protein family for GRN config ('microbial_opsins', 'gpcr_a', etc.)
        reference_table: Name of reference GRN table (e.g., 'mo_ref', 'gpcrdb_ref')
                        If None, must be specified based on protein_family
        similarity_threshold: Minimum sequence identity for alignment
        use_mmseqs: Whether to use MMseqs2 for alignment
        verbose: Whether to print progress messages
        
    Returns:
        Dictionary mapping chain IDs to pandas Series of GRN assignments
    """
    import pandas as pd
    from pathlib import Path
    import shutil
    import platform
    
    # Protos imports (same as test scripts)
    from protos.processing.grn.grn_processor import GRNProcessor
    from protos.processing.grn.grn_table_utils import (
        init_row_from_alignment,
        expand_annotation,
    )
    from protos.processing.grn.grn_utils import (
        get_seq, sort_grns_str, GRNConfigManager, parse_grn_str2float,
        parse_grn_float2str, get_grn_interval
    )
    from protos.processing.sequence.seq_alignment import (
        init_aligner,
        align_blosum62,
        mmseqs2_align2,
        format_alignment
    )
    
    if self.data is None:
        raise ValueError("No structure data loaded")
    
    # Initialize GRN processor
    grn_proc = GRNProcessor(paths=self.paths)
    
    # Determine reference table
    if reference_table is None:
        if protein_family == 'microbial_opsins':
            reference_table = 'mo_ref'
        elif protein_family == 'gpcr_a':
            reference_table = 'gpcrdb_ref'
        else:
            raise ValueError(f"No default reference table for protein_family '{protein_family}'. "
                           "Please specify reference_table parameter.")
    
    # Extract sequences from structure (same as test_grn_assignment.py line 60)
    if verbose:
        print("\nExtracting sequences from structures...")
    sequences = self.get_seq_dict()
    
    if not sequences:
        if verbose:
            print("No sequences found in structure data")
        return {}
    
    # Load reference GRN table (same as test_grn_assignment.py line 72)
    if verbose:
        print(f"\nLoading reference GRN table: {reference_table}")
    
    ref_file = grn_proc.path_ref_dir / f"{reference_table}.csv"
    if not ref_file.exists():
        raise FileNotFoundError(f"Reference GRN table not found: {ref_file}")
    
    grn_proc.data = grn_proc.load_reference_table(reference_table)
    grn_proc.ids = grn_proc.data.index.tolist()
    
    if verbose:
        print(f"Loaded {len(grn_proc.data)} sequences from reference table")
        print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Get reference sequences
    ref_sequences = grn_proc.get_seq_dict()
    
    # Clean up temp directory before alignment (same as test_grn_assignment.py line 82)
    mmseqs_tmp_dir = Path("temp/mmseqs_tmp")
    if mmseqs_tmp_dir.exists():
        if verbose:
            print("Cleaning up existing MMseqs2 temp directory...")
        if platform.system() == "Windows":
            try:
                for item in mmseqs_tmp_dir.iterdir():
                    if item.is_symlink():
                        item.unlink()
                shutil.rmtree(mmseqs_tmp_dir)
            except Exception as e:
                if verbose:
                    print(f"Warning: Could not fully clean temp directory: {e}")
                import subprocess
                try:
                    subprocess.run(['wsl', 'rm', '-rf', 'temp/mmseqs_tmp'], check=True)
                    if verbose:
                        print("Cleaned up using WSL")
                except:
                    if verbose:
                        print("Could not clean up temp directory, continuing anyway...")
        else:
            shutil.rmtree(mmseqs_tmp_dir)
    
    # Perform alignment (same as test_grn_assignment.py line 111)
    if verbose:
        print("\nPerforming MMseqs2 alignment...")
    alignment_df = mmseqs2_align2(sequences, ref_sequences)
    
    if alignment_df is None:
        if verbose:
            print("ERROR: MMseqs2 alignment failed")
        return {}
    
    # Filter by sequence identity
    alignment_df = alignment_df[alignment_df['sequence_identity'] > similarity_threshold]
    if verbose:
        print(f"Found {len(alignment_df)} alignments with > {similarity_threshold:.0%} sequence identity")
    
    # Get unique queries
    queries = alignment_df['query_id'].unique().tolist()
    if verbose:
        print(f"\nProcessing {len(queries)} unique queries")
    
    # Initialize aligner and config (same as test_grn_assignment.py line 127)
    aligner = init_aligner()
    config = GRNConfigManager(paths=self.paths)
    grn_config_str = config.get_config(protein_family=protein_family, strict=True)
    
    grns_str_str = []
    if grn_config_str:
        for region_name, (start_grn, end_grn) in grn_config_str.items():
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_str.extend(region_grns)
    
    # Remove duplicates and sort
    grns_str_str = list(set(grns_str_str))
    grns_str_str = sort_grns_str(grns_str_str)
    
    # Initialize GRN column if not exists
    if 'grn' not in self.data.columns:
        self.data['grn'] = None
    
    # Process each query (same pattern as test_grn_assignment.py line 147)
    rows = {}
    if verbose:
        print("\nAnnotating sequences...")
    
    for i, q in enumerate(queries):
        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing {i+1}/{len(queries)}: {q}")
            print(f"{'='*60}")
        
        # Get best sequence match
        query_alignments = alignment_df[alignment_df['query_id'] == q]
        best_alignment = query_alignments.loc[query_alignments['sequence_identity'].idxmax()]
        ref_id = best_alignment['target_id']
        
        if verbose:
            print(f"Best match: {ref_id} (identity: {best_alignment['sequence_identity']:.1%})")
        
        # Get sequences
        test_seq = sequences[q]
        ref_seq = ref_sequences[ref_id]
        
        # Align sequences
        alignment = align_blosum62(test_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        
        # Create initial annotation (same as test_grn_assignment.py line 172)
        ref_row = grn_proc.data.loc[ref_id]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        grns = [grn for grn in grns_str_str if grn in new_row.index]
        new_row = new_row[grns]
        
        if verbose:
            print(f"\nInitial annotation has {len(new_row)} positions")
        
        # Expand annotation (same as test_grn_assignment.py line 189)
        try:
            grn_list, rn_list, missing = expand_annotation(
                new_row,
                test_seq.replace('-', ''),
                formatted,
                max_alignment_gap=1,
                protein_family=protein_family,
                verbose=1 if verbose else 0
            )
            
            # Create final row
            final_row = pd.Series(dict(zip(grn_list, rn_list)))
            rows[q] = final_row
            
            # Now assign GRN values to structure data
            # Parse chain ID from query (format: PDBID_CHAIN)
            pdb_id, chain_id = q.rsplit('_', 1)
            
            # Create mapping from residue position to GRN
            # rn_list contains entries like 'M1', 'L2', etc.
            # We need to extract position number and map to GRN
            residue_grn_map = {}
            for rn, grn in zip(rn_list, grn_list):
                # Extract position number from residue notation (e.g., 'M1' -> 1)
                import re
                match = re.match(r'([A-Z])(\d+)', rn)
                if match:
                    res_type = match.group(1)
                    res_pos = int(match.group(2))
                    residue_grn_map[res_pos] = grn
            
            # Assign GRN values to matching residues in structure
            chain_mask = (self.data['pdb_id'] == pdb_id) & (self.data['auth_chain_id'] == chain_id)
            
            for idx in self.data[chain_mask].index:
                seq_id = self.data.loc[idx, 'auth_seq_id']
                if seq_id in residue_grn_map:
                    self.data.loc[idx, 'grn'] = residue_grn_map[seq_id]
            
            if verbose:
                assigned_count = self.data[chain_mask & self.data['grn'].notna()].shape[0]
                print(f"\nAssigned GRN to {assigned_count} residues in structure")
                
        except Exception as e:
            if verbose:
                print(f"Error in annotation: {e}")
            # Store partial results
            rows[q] = new_row
    
    # Summary
    if verbose:
        total_grn_residues = self.data[self.data['grn'].notna()].shape[0]
        print(f"\n\nTotal residues with GRN annotations: {total_grn_residues}")
        
        # Show some key positions
        key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
        print("\nKey positions found in structure:")
        for pos in key_positions:
            residues = self.data[self.data['grn'] == pos]
            if not residues.empty:
                for _, res in residues.iterrows():
                    print(f"  {pos}: {res['pdb_id']}_{res['auth_chain_id']} "
                          f"{res['auth_comp_id']}{res['auth_seq_id']}")
    
    return rows
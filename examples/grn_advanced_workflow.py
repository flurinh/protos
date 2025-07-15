#!/usr/bin/env python3
"""
Advanced GRN Assignment Workflow
================================

This script demonstrates an advanced 5-phase workflow for comprehensive
GRN (Generic Residue Numbering) assignment to protein sequences.

Workflow Phases:
1. Initial Alignment - Find best reference and transfer GRNs
2. Terminal Extensions - Extend GRNs to N and C terminals  
3. Missing Standard GRNs - Fill in conserved positions
4. Gap and Loop Annotation - Handle insertions/deletions
5. Validation and Output - Quality control and save results

Key Features:
- Comprehensive GRN assignment including terminals and loops
- Handles gaps and insertions in sequences
- Validates assignments and checks coverage
- Detailed logging of each phase

Usage:
    python examples/grn_advanced_workflow.py
    
Output:
    - Complete GRN assignment saved to: data/grn/tables/complete_workflow_result.csv
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.grn import GRNProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.sequence import (
    init_aligner, align_blosum62, format_alignment, mmseqs2_align2
)
from protos.processing.grn.grn_utils import (
    get_seq, sort_grns_str, GRNConfigManager, parse_grn_str2float, parse_grn_float2str
)
from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment, expand_annotation, init_grn_intervals
)
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)8s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class CompleteGRNWorkflow:
    """Complete GRN assignment workflow with real data."""
    
    def __init__(self):
        # Set up ProtosPaths to use the project's data directory
        from protos.io.paths.path_config import ProtosPaths
        
        # Use the project's data directory (non-deprecated approach)
        data_dir = Path(__file__).parent.parent / "data"  # Go up to project root
        paths = ProtosPaths(data_root=str(data_dir))
        
        # Initialize processors
        self.grn_proc = GRNProcessor(name="complete_grn_workflow", paths=paths)
        self.seq_proc = SequenceProcessor(name="complete_seq_workflow", paths=paths)
        self.aligner = init_aligner()
        
        # GRN configuration for microbial opsins
        # For this demo, we'll use standard positions
        self.standard_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
        
    def load_real_data(self) -> Dict[str, str]:
        """Load real opsin sequences."""
        sequence_name = "opsin_sequences_from_yaml"
        
        try:
            # Use the sequence processor to load the data
            # No hardcoded paths - processor knows where to look!
            if self.seq_proc.entity_exists(sequence_name):
                sequences = self.seq_proc.load_entity(sequence_name)
                if isinstance(sequences, str):
                    sequences = {sequence_name: sequences}
            else:
                sequences = self.seq_proc.load_sequences(sequence_name)
            
            logger.info(f"Loaded {len(sequences)} opsin sequences")
            return sequences
        except FileNotFoundError:
            logger.error(f"Sequence file '{sequence_name}' not found")
            logger.error(f"Expected location: {self.seq_proc.path_fasta_dir}/{sequence_name}.fasta")
            return {}
    
    def load_reference_data(self) -> pd.DataFrame:
        """Load or create reference GRN table."""
        try:
            # Try to load existing reference
            self.grn_proc.load_grn_table("ref/mo_ref")
            logger.info(f"Loaded reference table with {len(self.grn_proc.data)} proteins")
            return self.grn_proc.data
        except:
            # Create minimal reference
            logger.info("Creating reference GRN table")
            return self.create_reference_table()
    
    def create_reference_table(self) -> pd.DataFrame:
        """Create reference GRN table for microbial opsins."""
        # Bacteriorhodopsin reference
        br_grns = {}
        
        # TM1 (positions 36-62)
        tm1_seq = "WIWLALGTALMGLGTLYFLVKG"
        for i, aa in enumerate(tm1_seq):
            grn = f"1.{36+i}"
            br_grns[grn] = aa
        
        # TM2 (positions 39-65)
        tm2_seq = "DPAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVP"
        for i, aa in enumerate(tm2_seq[:27]):
            grn = f"2.{39+i}"
            br_grns[grn] = aa
        
        # TM3 (positions 39-65)
        tm3_seq = "NQDPIYWARYADWLFTTPLLLLDLALLAK"
        for i, aa in enumerate(tm3_seq[:27]):
            grn = f"3.{39+i}"
            br_grns[grn] = aa
        
        # TM4-7 conserved positions
        conserved = {
            '4.50': 'F', '5.50': 'G', '6.48': 'W', '6.50': 'P', '7.50': 'K'
        }
        br_grns.update(conserved)
        
        # Create DataFrame
        df = pd.DataFrame([br_grns], index=['BR_reference'])
        df = df.fillna('-')
        
        # Sort columns - don't use sort_grns_str since it converts to x notation
        grn_cols = [col for col in df.columns if '.' in col]
        def sort_key(grn):
            try:
                parts = grn.split('.')
                if len(parts) == 2:
                    return (int(parts[0]), int(parts[1]))
                return (999, 999)
            except:
                return (999, 999)
        sorted_cols = sorted(grn_cols, key=sort_key)
        df = df[sorted_cols]
        
        return df
    
    def demonstrate_workflow(self, sequence_name: str, query_seq: str):
        """Demonstrate the complete GRN assignment workflow."""
        logger.info("="*80)
        logger.info(f"COMPLETE GRN WORKFLOW FOR: {sequence_name}")
        logger.info("="*80)
        
        # Phase 1: Initial Alignment
        aligned_positions = self.phase1_initial_alignment(query_seq)
        
        # Phase 2: Terminal Extensions
        extended_positions = self.phase2_terminal_extensions(query_seq, aligned_positions)
        
        # Phase 3: Missing Standard GRNs
        filled_positions = self.phase3_missing_std_grns(query_seq, extended_positions)
        
        # Phase 4: Gap and Loop Annotation
        complete_positions = self.phase4_gaps_and_loops(query_seq, filled_positions)
        
        # Phase 5: Validation and Output
        final_grn = self.phase5_validation(sequence_name, query_seq, complete_positions)
        
        return final_grn
    
    def phase1_initial_alignment(self, query_seq: str) -> pd.Series:
        """Phase 1: Process initial alignment results."""
        logger.info("\n" + "="*60)
        logger.info("PHASE 1: INITIAL ALIGNMENT")
        logger.info("="*60)
        
        # Get reference sequences
        ref_sequences = {}
        for protein_id in self.grn_proc.ids:
            seq = get_seq(protein_id, self.grn_proc.data)
            if seq and len(seq) > 50:
                ref_sequences[protein_id] = seq
        
        logger.info(f"Using {len(ref_sequences)} reference sequences")
        
        # Find best match
        query_dict = {'query': query_seq}
        best_ref = list(ref_sequences.keys())[0]  # Default
        
        try:
            hits = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_sequences)
            if hits is not None and not hits.empty:
                best_ref = hits.iloc[0]['target_id']
                logger.info(f"Best match: {best_ref}")
        except:
            logger.info(f"Using default reference: {best_ref}")
        
        # Perform alignment
        ref_seq = ref_sequences[best_ref]
        alignment = align_blosum62(query_seq, ref_seq, self.aligner)
        formatted = format_alignment(alignment)
        
        logger.info(f"Alignment score: {alignment.score}")
        
        # Transfer GRN annotations
        ref_row = self.grn_proc.data.loc[best_ref]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        # Initialize GRN row from alignment
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        logger.info(f"Transferred {len(new_row)} GRN positions")
        
        # Show some assignments
        logger.info("\nExample GRN assignments:")
        for i, (grn, (res, pos)) in enumerate(new_row.items()):
            if i < 5 and res != '-':
                logger.info(f"  {res}{pos} -> {grn}")
        
        return new_row
    
    def phase2_terminal_extensions(self, query_seq: str, aligned_positions: pd.Series) -> pd.Series:
        """Phase 2: Extend GRN assignments to terminals."""
        logger.info("\n" + "="*60)
        logger.info("PHASE 2: TERMINAL EXTENSIONS")
        logger.info("="*60)
        
        # Get first and last aligned positions
        aligned_pos = [(pos, grn) for grn, (res, pos) in aligned_positions.items() if res != '-']
        if not aligned_pos:
            return aligned_positions
            
        aligned_pos.sort(key=lambda x: x[0])
        first_pos, first_grn = aligned_pos[0]
        last_pos, last_grn = aligned_pos[-1]
        
        logger.info(f"First aligned: position {first_pos} -> {first_grn}")
        logger.info(f"Last aligned: position {last_pos} -> {last_grn}")
        
        # N-terminal extension
        logger.info(f"\nN-terminal extension (positions 1-{first_pos-1}):")
        
        # Parse first GRN
        tm, pos = first_grn.split('.')
        start_pos = int(pos)
        
        # Extend backwards
        for i in range(first_pos - 1, 0, -1):
            new_pos = start_pos - (first_pos - i)
            if new_pos > 0:
                new_grn = f"{tm}.{new_pos}"
                res = query_seq[i-1]
                aligned_positions[new_grn] = (res, i)
                if i >= first_pos - 5:  # Show first few
                    logger.info(f"  {res}{i} -> {new_grn}")
        
        # C-terminal extension
        logger.info(f"\nC-terminal extension (positions {last_pos+1}-{len(query_seq)}):")
        
        # Parse last GRN
        tm, pos = last_grn.split('.')
        end_pos = int(pos)
        
        # Extend forwards
        for i in range(last_pos + 1, len(query_seq) + 1):
            new_pos = end_pos + (i - last_pos)
            if new_pos <= 99:  # Keep within reasonable range
                new_grn = f"{tm}.{new_pos}"
                res = query_seq[i-1]
                aligned_positions[new_grn] = (res, i)
                if i <= last_pos + 5:  # Show first few
                    logger.info(f"  {res}{i} -> {new_grn}")
        
        logger.info(f"\nTotal positions after extension: {len(aligned_positions)}")
        return aligned_positions
    
    def phase3_missing_std_grns(self, query_seq: str, aligned_positions: pd.Series) -> pd.Series:
        """Phase 3: Assign missing standard GRN positions."""
        logger.info("\n" + "="*60)
        logger.info("PHASE 3: MISSING STANDARD GRNs")
        logger.info("="*60)
        
        # Standard GRNs that should be present
        standard_grns = self.standard_grns
        assigned_grns = set(aligned_positions.index)
        missing_std = [g for g in standard_grns if g not in assigned_grns]
        
        logger.info(f"Standard positions: {standard_grns}")
        logger.info(f"Currently assigned: {sorted([g for g in standard_grns if g in assigned_grns])}")
        logger.info(f"Missing: {missing_std}")
        
        # Try to assign missing positions
        for missing_grn in missing_std:
            tm, pos = missing_grn.split('.')
            tm_num = int(tm)
            pos_num = int(pos)
            
            # Find nearest assigned GRN in same TM
            same_tm_grns = [g for g in assigned_grns if g.startswith(f"{tm}.")]
            if not same_tm_grns:
                logger.info(f"  {missing_grn}: No reference in TM{tm}")
                continue
            
            # Find closest position
            distances = []
            for grn in same_tm_grns:
                _, p = grn.split('.')
                distances.append((abs(int(p) - pos_num), grn))
            
            distances.sort()
            closest_dist, closest_grn = distances[0]
            
            # Get position of closest GRN
            closest_res, closest_pos = aligned_positions[closest_grn]
            
            # Calculate expected position
            _, closest_p = closest_grn.split('.')
            offset = pos_num - int(closest_p)
            expected_pos = closest_pos + offset
            
            # Check if position is valid
            if 1 <= expected_pos <= len(query_seq):
                res = query_seq[expected_pos - 1]
                aligned_positions[missing_grn] = (res, expected_pos)
                logger.info(f"  {missing_grn}: Assigned to {res}{expected_pos} (based on {closest_grn})")
            else:
                logger.info(f"  {missing_grn}: Expected position {expected_pos} out of range")
        
        return aligned_positions
    
    def phase4_gaps_and_loops(self, query_seq: str, aligned_positions: pd.Series) -> pd.Series:
        """Phase 4: Annotate gaps and loops."""
        logger.info("\n" + "="*60)
        logger.info("PHASE 4: GAP AND LOOP ANNOTATION")
        logger.info("="*60)
        
        # Find all assigned positions
        assigned_pos = set()
        for grn, (res, pos) in aligned_positions.items():
            if res != '-':
                assigned_pos.add(pos)
        
        # Find missing positions
        all_positions = set(range(1, len(query_seq) + 1))
        missing_positions = sorted(all_positions - assigned_pos)
        
        logger.info(f"Total positions: {len(query_seq)}")
        logger.info(f"Assigned positions: {len(assigned_pos)}")
        logger.info(f"Missing positions: {len(missing_positions)}")
        
        if not missing_positions:
            return aligned_positions
        
        # Group consecutive missing positions
        intervals = []
        current = [missing_positions[0]]
        
        for pos in missing_positions[1:]:
            if pos == current[-1] + 1:
                current.append(pos)
            else:
                intervals.append(current)
                current = [pos]
        intervals.append(current)
        
        logger.info(f"\nGrouped into {len(intervals)} intervals:")
        for i, interval in enumerate(intervals[:5]):  # Show first 5
            logger.info(f"  Interval {i+1}: positions {interval[0]}-{interval[-1]} ({len(interval)} residues)")
        
        # Annotate each interval
        for interval in intervals:
            start_pos = interval[0]
            end_pos = interval[-1]
            
            # Find flanking GRNs
            before_grn = None
            after_grn = None
            
            for grn, (res, pos) in aligned_positions.items():
                if pos == start_pos - 1:
                    before_grn = grn
                elif pos == end_pos + 1:
                    after_grn = grn
            
            if before_grn and after_grn:
                # Determine if gap or loop
                before_tm = before_grn.split('.')[0]
                after_tm = after_grn.split('.')[0]
                
                if before_tm == after_tm:
                    # Gap within TM
                    logger.info(f"  Gap in TM{before_tm}: positions {start_pos}-{end_pos}")
                    
                    # Annotate gap positions
                    for i, pos in enumerate(interval):
                        gap_grn = f"{before_grn}.{i+1:03d}"
                        res = query_seq[pos-1]
                        aligned_positions[gap_grn] = (res, pos)
                else:
                    # Loop between TMs
                    logger.info(f"  Loop between TM{before_tm} and TM{after_tm}: positions {start_pos}-{end_pos}")
                    
                    # Annotate loop positions
                    for i, pos in enumerate(interval):
                        loop_grn = f"{before_tm}{after_tm}.{i+1:02d}"
                        res = query_seq[pos-1]
                        aligned_positions[loop_grn] = (res, pos)
        
        return aligned_positions
    
    def phase5_validation(self, sequence_name: str, query_seq: str, 
                         aligned_positions: pd.Series) -> pd.Series:
        """Phase 5: Validate and output results."""
        logger.info("\n" + "="*60)
        logger.info("PHASE 5: VALIDATION AND OUTPUT")
        logger.info("="*60)
        
        # Check for duplicate positions
        position_counts = {}
        for grn, (res, pos) in aligned_positions.items():
            if res != '-':
                if pos in position_counts:
                    position_counts[pos].append(grn)
                else:
                    position_counts[pos] = [grn]
        
        duplicates = {pos: grns for pos, grns in position_counts.items() if len(grns) > 1}
        if duplicates:
            logger.warning(f"Found {len(duplicates)} positions with multiple GRNs")
            for pos, grns in list(duplicates.items())[:3]:
                logger.warning(f"  Position {pos}: {grns}")
        
        # Check key positions
        key_positions = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
        assigned_keys = [k for k in key_positions if k in aligned_positions]
        logger.info(f"\nKey positions assigned: {len(assigned_keys)}/{len(key_positions)}")
        
        for pos in key_positions:
            if pos in aligned_positions:
                res, seq_pos = aligned_positions[pos]
                logger.info(f"  {pos}: {res}{seq_pos}")
        
        # Calculate coverage
        assigned_count = len([1 for res, pos in aligned_positions.values() if res != '-'])
        coverage = assigned_count / len(query_seq) * 100
        
        logger.info(f"\nSummary statistics:")
        logger.info(f"  - Sequence length: {len(query_seq)}")
        logger.info(f"  - Assigned positions: {assigned_count}")
        logger.info(f"  - Coverage: {coverage:.1f}%")
        
        # Create final GRN series
        grn_dict = {}
        for grn, (res, pos) in aligned_positions.items():
            if res != '-':
                grn_dict[grn] = f"{res}{pos}"
        
        final_grn = pd.Series(grn_dict, name=sequence_name)
        
        # Add missing columns
        for col in self.grn_proc.grns:
            if col not in final_grn.index:
                final_grn[col] = '-'
        
        # Sort - don't use sort_grns_str since it converts to x notation
        grn_indices = [idx for idx in final_grn.index if '.' in idx]
        def sort_key(grn):
            try:
                parts = grn.split('.')
                if len(parts) == 2:
                    return (int(parts[0]), int(parts[1]))
                elif len(parts) == 3:  # Loop positions like "12.003"
                    return (int(parts[0]), int(parts[1]) + float(f"0.{parts[2]}"))
                return (999, 999)
            except:
                return (999, 999)
        sorted_indices = sorted(grn_indices, key=sort_key)
        final_grn = final_grn[sorted_indices]
        
        logger.info("\n" + "="*60)
        logger.info("WORKFLOW COMPLETE")
        logger.info("="*60)
        
        return final_grn


def main():
    """Run the complete workflow demonstration with real data."""
    workflow = CompleteGRNWorkflow()
    
    # Load real opsin sequences
    sequences = workflow.load_real_data()
    if not sequences:
        logger.error("No sequences loaded")
        return
    
    # Load reference data
    ref_data = workflow.load_reference_data()
    workflow.grn_proc.data = ref_data
    workflow.grn_proc.ids = ref_data.index.tolist()
    workflow.grn_proc.grns = ref_data.columns.tolist()
    
    # Process first sequence as demonstration
    seq_name, query_seq = list(sequences.items())[0]
    logger.info(f"\nProcessing: {seq_name} ({len(query_seq)} AA)")
    
    # Run complete workflow
    final_grn = workflow.demonstrate_workflow(seq_name, query_seq)
    
    # Save result
    result_df = pd.DataFrame([final_grn])
    workflow.grn_proc.data = result_df
    workflow.grn_proc.ids = result_df.index.tolist()
    workflow.grn_proc.grns = result_df.columns.tolist()
    
    saved_path = workflow.grn_proc.save_grn_table("complete_workflow_result")
    logger.info(f"\nSaved result to: {saved_path}")


if __name__ == "__main__":
    main()
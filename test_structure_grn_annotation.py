#!/usr/bin/env python3
"""Annotate GPCR structures with GRN (Generic Residue Numbering) positions."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos

from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.io.ingest.structure_loader import StructureLoader


def ensure_data_root() -> Path:
    """Set up data root directory."""
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


# Test with a few well-known GPCR structures
STRUCTURE_IDS = ["3sn6", "5d5a", "6b73", "4daj"]  # Beta2, A2A, mGluR5, M3 muscarinic
GPCR_DATASET = "gpcr_structures_grn"
SEQUENCE_DATASET = "gpcr_chains_for_grn"
GRN_TABLE_NAME = "gpcr_structure_grn_annotations"
REFERENCE_TABLE = "gpcrdb_ref"
PROTEIN_FAMILY = "gpcr_a"


def ensure_structures(structure_ids: List[str]) -> None:
    """Download structures if not available."""
    loader = StructureLoader()
    try:
        loader.download_batch(structure_ids, dataset_name=GPCR_DATASET)
        print(f"✓ Ensured {len(structure_ids)} GPCR structures")
    except Exception as exc:
        print(f"Warning: could not refresh structures: {exc}")


def extract_and_register_sequences(
    struct_proc: StructureProcessor,
    structure_ids: List[str]
) -> Dict[str, str]:
    """Extract sequences from structures and register them."""
    # Ensure structures are loaded
    struct_proc.load_dataset(GPCR_DATASET, return_format="dict")
    
    # Extract chain sequences and register them
    struct_proc.register_chain_sequences(
        structure_ids,
        dataset_prefix="gpcr_chain_dataset",
        create_dataset=True,
        overwrite=False,
    )
    
    # Collect all chain sequences
    all_sequences = {}
    for struct_id in structure_ids:
        chain_seqs = struct_proc.collect_chain_sequences([struct_id])
        for chains in chain_seqs.values():
            for chain_data in chains.values():
                seq_name = chain_data["entity_name"]
                sequence = chain_data["sequence"]
                all_sequences[seq_name] = sequence
    
    print(f"✓ Extracted {len(all_sequences)} total chain sequences from {len(structure_ids)} structures")
    return all_sequences


def filter_gpcr_chains(
    seq_proc: SequenceProcessor,
    all_sequences: Dict[str, str],
    reference_gpcr_sequence: Optional[str] = None,
    similarity_threshold: float = 0.35
) -> Dict[str, str]:
    """Filter sequences to identify GPCR chains using sequence similarity."""
    from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine
    
    # If no reference provided, use a known GPCR sequence (e.g., from beta2-adrenergic receptor)
    if reference_gpcr_sequence is None:
        # This is human beta2-adrenergic receptor sequence
        reference_gpcr_sequence = "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL"
    
    engine = SequenceAlignmentEngine()
    gpcr_sequences = {}
    
    print(f"\nFiltering for GPCR chains (threshold: {similarity_threshold})...")
    
    for seq_name, sequence in all_sequences.items():
        # Skip very short sequences
        if len(sequence) < 100:
            continue
            
        # Align against reference GPCR
        result = engine.align_pairwise(
            seq_name, sequence,
            "reference_gpcr", reference_gpcr_sequence
        )
        
        # Calculate normalized score
        normalized_score = result.score / max(len(sequence), len(reference_gpcr_sequence))
        
        if normalized_score >= similarity_threshold:
            gpcr_sequences[seq_name] = sequence
            print(f"  ✓ {seq_name}: score={normalized_score:.3f} (GPCR-like)")
        else:
            print(f"  ✗ {seq_name}: score={normalized_score:.3f} (not GPCR)")
    
    print(f"\nIdentified {len(gpcr_sequences)} GPCR chains out of {len(all_sequences)} total chains")
    return gpcr_sequences


def create_sequence_dataset(
    seq_proc: SequenceProcessor,
    sequences: Dict[str, str]
) -> None:
    """Create a sequence dataset for GRN annotation."""
    # Save sequences as a dataset
    seq_proc.save_sequences(
        sequences,
        output_file=SEQUENCE_DATASET,
        dataset_name=SEQUENCE_DATASET,
        metadata={
            "description": "GPCR chain sequences for GRN annotation",
            "source": "structure extraction"
        },
        materialize_entities=True,
    )
    print(f"✓ Created sequence dataset '{SEQUENCE_DATASET}' with {len(sequences)} sequences")


def annotate_sequences_with_grn(seq_proc: SequenceProcessor) -> pd.DataFrame:
    """Annotate sequences with GRN positions."""
    print("\nAnnotating sequences with GRN...")
    
    grn_table, summary = seq_proc.annotate_with_grn(
        dataset_name=SEQUENCE_DATASET,
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        output_table=GRN_TABLE_NAME,
        return_summary=True,
        allow_create=True,
    )
    
    # Print annotation summary
    print("\nGRN Annotation Summary:")
    print(f"Total sequences: {summary['global']['total']}")
    print(f"Successfully annotated: {summary['global']['annotated']}")
    
    print("\nPer-sequence coverage:")
    for seq_id, info in summary["per_sequence"].items():
        coverage = info.get("coverage", 0)
        status = info.get("status", "unknown")
        assigned = info.get("assigned_positions", 0)
        print(f"  {seq_id}: coverage={coverage:.1%}, positions={assigned}, status={status}")
    
    return grn_table


def map_grn_to_structures(
    struct_proc: StructureProcessor,
    grn_proc: GRNProcessor,
    sequences: Dict[str, str]
) -> Dict[str, int]:
    """Map GRN annotations back to structure residues."""
    # Load the GRN table
    grn_table = grn_proc.load_table(GRN_TABLE_NAME)
    
    # Track how many residues get GRN annotations per structure
    annotation_counts = {}
    
    # For each sequence, map GRN back to its source structure
    for seq_name in sequences:
        # Parse structure and chain from sequence name (e.g., "3sn6_chain_A")
        parts = seq_name.split("_chain_")
        if len(parts) != 2:
            continue
            
        struct_id, chain_id = parts
        
        # Get GRN data for this sequence
        if seq_name not in grn_table.index:
            print(f"⚠ No GRN data for {seq_name}")
            continue
            
        grn_data = grn_table.loc[seq_name]
        
        # Load structure data
        struct_df = struct_proc.load_entity(struct_id)
        if struct_df is None:
            continue
            
        # Filter to specific chain
        chain_df = struct_df[struct_df['auth_chain_id'] == chain_id].copy()
        
        # Create a mapping from sequence position to GRN
        seq_pos_to_grn = {}
        for grn_pos, value in grn_data.items():
            if value != '-' and isinstance(value, str) and len(value) > 1:
                # Extract position number from value (e.g., "N51" -> 51)
                try:
                    pos = int(value[1:])
                    seq_pos_to_grn[pos] = grn_pos
                except (ValueError, IndexError):
                    continue
        
        # Add GRN column to structure dataframe
        grn_annotations = []
        count = 0
        
        # Get unique residues in order
        residue_positions = chain_df.groupby(['auth_seq_id', 'insertion']).first().reset_index()
        
        for idx, row in residue_positions.iterrows():
            seq_id = row['auth_seq_id']
            # Map auth_seq_id to sequence position (1-indexed)
            # This assumes auth_seq_id corresponds to position in sequence
            grn = seq_pos_to_grn.get(seq_id, '-')
            if grn != '-':
                count += 1
            
            # Apply to all atoms of this residue
            mask = (chain_df['auth_seq_id'] == seq_id) & \
                   (chain_df['insertion'] == row['insertion'])
            chain_df.loc[mask, 'grn'] = grn
        
        # Update the structure dataframe with GRN annotations
        struct_df.loc[chain_df.index, 'grn'] = chain_df['grn']
        
        # Save the annotated structure
        struct_proc.frames[struct_id] = struct_df
        annotation_counts[f"{struct_id}_{chain_id}"] = count
        
        print(f"✓ Mapped {count} GRN positions to {struct_id} chain {chain_id}")
    
    return annotation_counts


def demonstrate_grn_mapping(
    struct_proc: StructureProcessor,
    struct_id: str = "5d5a",
    chain_id: str = "A"
) -> None:
    """Demonstrate GRN mapping by showing key conserved positions."""
    print(f"\nDemonstrating GRN mapping for {struct_id} chain {chain_id}:")
    
    # Load annotated structure
    struct_df = struct_proc.load_entity(struct_id)
    if struct_df is None:
        print(f"Structure {struct_id} not found")
        return
        
    # Filter to chain and CA atoms
    ca_atoms = struct_df[
        (struct_df['auth_chain_id'] == chain_id) & 
        (struct_df['atom_name'] == 'CA') &
        (struct_df.get('grn', '-') != '-')
    ].copy()
    
    if ca_atoms.empty:
        print("No GRN annotations found")
        return
    
    # Show some key conserved GPCR positions
    key_positions = {
        "1.50": "Most conserved position in TM1 (usually N)",
        "2.50": "Most conserved position in TM2 (usually D)", 
        "3.50": "DRY motif position (usually R)",
        "4.50": "Most conserved position in TM4 (usually W)",
        "5.50": "Most conserved position in TM5 (usually P)",
        "6.50": "Most conserved position in TM6 (usually P)",
        "7.50": "NPxxY motif position (usually P)",
    }
    
    print("\nKey conserved GPCR positions:")
    for grn_pos, description in key_positions.items():
        residues = ca_atoms[ca_atoms['grn'] == grn_pos]
        if not residues.empty:
            row = residues.iloc[0]
            print(f"  GRN {grn_pos}: {row['res_name']} {row['auth_seq_id']} - {description}")
    
    # Summary statistics
    grn_coverage = ca_atoms['grn'].value_counts()
    print(f"\nTotal CA atoms with GRN: {len(ca_atoms)}")
    print(f"Unique GRN positions: {len(grn_coverage)}")
    
    # Show GRN distribution by TM region
    tm_regions = {}
    for grn in ca_atoms['grn'].unique():
        if '.' in grn and grn[0].isdigit():
            tm = grn.split('.')[0]
            tm_regions[tm] = tm_regions.get(tm, 0) + 1
    
    print("\nGRN positions by TM helix:")
    for tm in sorted(tm_regions.keys(), key=int):
        print(f"  TM{tm}: {tm_regions[tm]} positions")


def save_annotated_structures(
    struct_proc: StructureProcessor,
    structure_ids: List[str],
    output_prefix: str = "grn_annotated"
) -> None:
    """Save structures with GRN annotations."""
    output_dir = Path("data/structure/annotated")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for struct_id in structure_ids:
        struct_df = struct_proc.frames.get(struct_id)
        if struct_df is None or 'grn' not in struct_df.columns:
            continue
            
        # Save as pickle for easy loading
        output_path = output_dir / f"{output_prefix}_{struct_id}.pkl"
        struct_df.to_pickle(output_path)
        print(f"✓ Saved GRN-annotated structure to {output_path}")


def main() -> None:
    """Main workflow for structure GRN annotation."""
    # Initialize
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")
    
    # Initialize processors
    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    
    # Step 1: Ensure structures are available
    ensure_structures(STRUCTURE_IDS)
    
    # Step 2: Extract sequences from structures
    all_sequences = extract_and_register_sequences(struct_proc, STRUCTURE_IDS)
    
    # Step 3: Filter for GPCR chains using sequence similarity
    gpcr_sequences = filter_gpcr_chains(seq_proc, all_sequences)
    
    if not gpcr_sequences:
        print("❌ No GPCR chains identified!")
        return
    
    # Step 4: Create sequence dataset with only GPCR chains
    create_sequence_dataset(seq_proc, gpcr_sequences)
    
    # Step 5: Annotate sequences with GRN
    grn_table = annotate_sequences_with_grn(seq_proc)
    
    # Step 6: Map GRN annotations back to structures
    print("\nMapping GRN annotations to structures...")
    annotation_counts = map_grn_to_structures(struct_proc, grn_proc, gpcr_sequences)
    
    # Step 6: Demonstrate the mapping
    demonstrate_grn_mapping(struct_proc, "5d5a", "A")
    
    # Step 7: Save annotated structures
    save_annotated_structures(struct_proc, STRUCTURE_IDS)
    
    # Summary
    print("\n=== Summary ===")
    print(f"Processed {len(STRUCTURE_IDS)} structures")
    print(f"Identified {len(gpcr_sequences)} GPCR chains from {len(all_sequences)} total chains")
    print(f"Total residues with GRN annotations: {sum(annotation_counts.values())}")
    
    # Show which structures have the most GRN coverage
    if annotation_counts:
        sorted_counts = sorted(annotation_counts.items(), key=lambda x: x[1], reverse=True)
        print("\nTop structures by GRN coverage:")
        for struct_chain, count in sorted_counts[:5]:
            print(f"  {struct_chain}: {count} positions")


if __name__ == "__main__":
    main()
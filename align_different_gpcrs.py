#!/usr/bin/env python3
"""
Align different GPCR structures using sequence-independent alignment.

Since different GPCRs have different sequences and residue numbering,
we need to use a more sophisticated alignment approach.
"""

import os
import sys
from pathlib import Path
import numpy as np

# Add project root to path
project_root = Path(__file__).parent.absolute()
sys.path.insert(0, str(project_root))


def main():
    print("=== Different GPCR Structure Alignment ===\n")
    
    # 1. Setup ProtosPaths
    print("1. Setting up ProtosPaths...")
    os.environ["PROTOS_DATA_ROOT"] = str(project_root / "data")
    
    # 2. Initialize processor
    from protos.processing.structure import StructureProcessor
    processor = StructureProcessor()
    
    # 3. Check available GPCRs
    print("\n2. Checking available GPCR structures...")
    entities = processor.list_entities()
    
    # Look for GPCR structures
    gpcr_ids = []
    for entity in entities:
        # Common GPCR PDB IDs
        if entity in ['1u19', '3sn6', '2rh1', '3pbl', '4djh', '5d5a']:
            gpcr_ids.append(entity)
    
    print(f"   Found {len(gpcr_ids)} GPCR structures: {', '.join(gpcr_ids)}")
    
    if len(gpcr_ids) < 2:
        print("\n   Downloading additional GPCRs...")
        from protos.io.ingest.structure_loader import StructureLoader
        loader = StructureLoader()
        
        # Download some well-known GPCRs
        additional_gpcrs = {
            "2rh1": "Beta2-adrenergic (inactive)",
            "5d5a": "M2 muscarinic receptor"
        }
        
        for pdb_id, desc in additional_gpcrs.items():
            if pdb_id not in entities:
                print(f"   Downloading {pdb_id} - {desc}")
                result = loader.download_and_register(pdb_id)
                if result:
                    # Convert to PKL
                    cif_path = Path(processor.paths.get_subdir_path('structure', 'structure_dir')) / f"{pdb_id}.cif"
                    if cif_path.exists():
                        from protos.io.formats.cif_utils import load_structure_from_cif
                        df = load_structure_from_cif(str(cif_path))
                        processor.save_entity(pdb_id, df, format='pkl')
                        gpcr_ids.append(pdb_id)
    
    if len(gpcr_ids) < 2:
        print("✗ Need at least 2 GPCRs. Exiting.")
        return
    
    # Use first two GPCRs
    ref_id = gpcr_ids[0]
    mob_id = gpcr_ids[1]
    
    # 4. Load and prepare structures
    print(f"\n3. Loading structures...")
    print(f"   Reference: {ref_id}")
    print(f"   Mobile: {mob_id}")
    
    # Load structures
    ref_df = processor.load_entity(ref_id)
    mob_df = processor.load_entity(mob_id)
    
    if ref_df is None or mob_df is None:
        print("✗ Failed to load structures")
        return
    
    # Clean structures
    processor.remove_hetatm(ref_id)
    processor.remove_hetatm(mob_id)
    
    # Get main chains
    def get_main_chain(df):
        df_reset = df.reset_index()
        ca_counts = df_reset[df_reset['atom_name'] == 'CA'].groupby('auth_chain_id').size()
        return ca_counts.idxmax() if len(ca_counts) > 0 else None
    
    ref_chain = get_main_chain(processor.frames[ref_id])
    mob_chain = get_main_chain(processor.frames[mob_id])
    
    print(f"\n   Main chains: {ref_id}:{ref_chain}, {mob_id}:{mob_chain}")
    
    # Filter to main chains
    processor.filter_by_chain(ref_id, [ref_chain], new_id=f"{ref_id}_main")
    processor.filter_by_chain(mob_id, [mob_chain], new_id=f"{mob_id}_main")
    
    # 5. Use CE-align or similar structural alignment
    print("\n4. Performing structural alignment...")
    
    # Since residue IDs don't match, we'll use coordinate-based alignment
    from protos.analysis.structure import extract_sequence
    
    # Extract sequences to show they're different
    ref_seq = extract_sequence(processor.frames[f"{ref_id}_main"], one_letter=True)
    mob_seq = extract_sequence(processor.frames[f"{mob_id}_main"], one_letter=True)
    
    print(f"\n   Reference sequence ({len(ref_seq)} residues):")
    print(f"   {ref_seq[:60]}...")
    print(f"\n   Mobile sequence ({len(mob_seq)} residues):")
    print(f"   {mob_seq[:60]}...")
    
    # Try structural alignment using external tools if available
    print("\n5. Attempting structural superposition...")
    
    # For now, we'll do a simple center-based alignment
    from protos.analysis.structure import (
        calculate_center_of_mass,
        calculate_radius_of_gyration
    )
    
    ref_com = calculate_center_of_mass(processor.frames[f"{ref_id}_main"])
    mob_com = calculate_center_of_mass(processor.frames[f"{mob_id}_main"])
    
    print(f"\n   Reference COM: [{ref_com[0]:.1f}, {ref_com[1]:.1f}, {ref_com[2]:.1f}]")
    print(f"   Mobile COM: [{mob_com[0]:.1f}, {mob_com[1]:.1f}, {mob_com[2]:.1f}]")
    
    # Calculate translation to align centers
    translation = ref_com - mob_com
    
    # Apply translation
    translated_id = f"{mob_id}_translated"
    processor._set_frame(translated_id, processor.frames[f"{mob_id}_main"].copy())
    processor.apply_transformation(translated_id, translation=translation)
    
    print(f"\n   Applied translation: [{translation[0]:.1f}, {translation[1]:.1f}, {translation[2]:.1f}]")
    
    # For better alignment, we would need:
    # 1. Secondary structure matching
    # 2. Transmembrane helix detection
    # 3. Sequence-independent structural alignment (CE, TM-align, etc.)
    
    print("\n   Note: For proper alignment of different GPCRs, consider using:")
    print("   - CE-align or TM-align for sequence-independent alignment")
    print("   - GPCR-specific alignment using conserved motifs (DRY, NPxxY)")
    print("   - Transmembrane helix-based alignment")
    
    # 6. Export aligned structures
    print("\n6. Exporting structures...")
    
    output_dir = Path("aligned_different_gpcrs")
    output_dir.mkdir(exist_ok=True)
    
    # Export structures
    try:
        ref_path = processor.export_entity(
            f"{ref_id}_main",
            output_dir / f"{ref_id}_reference.cif",
            format='cif'
        )
        print(f"   ✓ Reference: {ref_path}")
        
        mob_path = processor.export_entity(
            translated_id,
            output_dir / f"{mob_id}_translated.cif",
            format='cif'
        )
        print(f"   ✓ Translated: {mob_path}")
    except Exception as e:
        print(f"   ✗ Export failed: {e}")
    
    # 7. Analyze conserved regions
    print("\n7. GPCR-specific analysis...")
    
    # Look for conserved GPCR motifs
    conserved_motifs = {
        'DRY': 'DRY',  # End of TM3
        'CWxP': 'CW.P', # TM6
        'NPxxY': 'NP..Y'  # TM7
    }
    
    print("\n   Conserved GPCR motifs:")
    for name, pattern in conserved_motifs.items():
        ref_pos = ref_seq.find(pattern.replace('.', ''))
        mob_pos = mob_seq.find(pattern.replace('.', ''))
        
        if ref_pos >= 0:
            print(f"   {name} in {ref_id}: position {ref_pos}")
        if mob_pos >= 0:
            print(f"   {name} in {mob_id}: position {mob_pos}")
    
    print("\n✓ Basic alignment complete!")
    print(f"\nStructures exported to '{output_dir}' directory.")
    print("\nFor better alignment of different GPCRs:")
    print("1. Use specialized tools like GPCRdb structure superposition")
    print("2. Align based on transmembrane helices")
    print("3. Use conserved motifs as anchor points")


if __name__ == "__main__":
    main()
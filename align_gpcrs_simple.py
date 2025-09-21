#!/usr/bin/env python3
"""
Simple GPCR alignment using the correct workflow.
"""

import os
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.absolute()
sys.path.insert(0, str(project_root))


def main():
    print("=== Simple GPCR Alignment ===\n")
    
    # 1. Setup
    os.environ["PROTOS_DATA_ROOT"] = str(project_root / "data")

    from protos.io.ingest.structure_loader import StructureLoader
    from protos.processing.structure import StructureProcessor
    
    # 2. Download structures
    print("1. Downloading structures...")
    loader = StructureLoader()
    
    # Clear any bad registrations and re-download
    structures_to_align = ["2rh1", "3pbl"]  # Two beta-adrenergic receptors
    
    for pdb_id in structures_to_align:
        print(f"\n   Downloading {pdb_id}...")
        result = loader.download_and_register(pdb_id)
        if result:
            print(f"   ✓ Downloaded {result}")
        else:
            print(f"   ✗ Failed to download {pdb_id}")
    
    # 3. Process with StructureProcessor
    print("\n2. Loading and processing structures...")
    processor = StructureProcessor()
    
    # Canonical PKLs already created by loader; nothing else needed here
    
    # 4. Load structures
    print("\n3. Loading structures for alignment...")
    
    ref_id = structures_to_align[0]
    mob_id = structures_to_align[1]
    
    ref_df = processor.load_entity(ref_id)
    mob_df = processor.load_entity(mob_id)
    
    if ref_df is None or mob_df is None:
        print("✗ Failed to load structures")
        return
        
    print(f"   Reference ({ref_id}): {len(ref_df)} atoms")
    print(f"   Mobile ({mob_id}): {len(mob_df)} atoms")
    
    # 5. Clean structures
    print("\n4. Cleaning structures...")
    processor.remove_hetatm(ref_id)
    processor.remove_hetatm(mob_id)
    
    # 6. Simple alignment
    print("\n5. Performing alignment...")
    
    # Get cleaned structures
    ref_clean = processor.frames[ref_id]
    mob_clean = processor.frames[mob_id]
    
    print(f"   Clean reference: {len(ref_clean)} atoms")
    print(f"   Clean mobile: {len(mob_clean)} atoms")
    
    # Try alignment
    try:
        results = processor.align_structures(
            structure_ids=[mob_id],
            reference_id=ref_id,
            atom_selection='CA',
            save_aligned=True
        )
        
        for struct_id, result in results.items():
            if 'error' not in result:
                print(f"   ✓ Aligned {struct_id}: RMSD = {result['rmsd']:.2f} Å")
            else:
                print(f"   ✗ Failed {struct_id}: {result['error']}")
                
    except Exception as e:
        print(f"   ✗ Alignment error: {e}")
        
        # Try manual alignment
        print("\n   Attempting manual CA alignment...")
        from protos.analysis.structure.alignment import simple_align_structures
        
        try:
            rotation, translation, rmsd = simple_align_structures(
                ref_clean, mob_clean,
                atom_selection='CA'
            )
            print(f"   ✓ Manual alignment RMSD: {rmsd:.2f} Å")
            
            # Apply transformation
            mob_aligned = processor.apply_transformation(mob_id, rotation=rotation, translation=translation)
            processor._set_frame(f"{mob_id}_aligned", mob_aligned)
            
        except Exception as e2:
            print(f"   ✗ Manual alignment also failed: {e2}")
    
    # 7. Export results
    print("\n6. Exporting results...")
    
    output_dir = Path("gpcr_alignment_output")
    output_dir.mkdir(exist_ok=True)
    
    # Export reference
    try:
        ref_path = processor.export_entity(
            ref_id,
            output_dir / f"{ref_id}_reference.cif"
        )
        print(f"   ✓ Exported reference: {ref_path}")
    except Exception as e:
        print(f"   ✗ Export failed: {e}")
    
    # Export aligned if exists
    aligned_id = f"{mob_id}_aligned"
    if aligned_id in processor.frames:
        try:
            aln_path = processor.export_entity(
                aligned_id,
                output_dir / f"{mob_id}_aligned.cif"
            )
            print(f"   ✓ Exported aligned: {aln_path}")
        except Exception as e:
            print(f"   ✗ Export failed: {e}")
    
    print("\n✓ Done!")


if __name__ == "__main__":
    main()

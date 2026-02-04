#!/usr/bin/env python3
"""
Complete BoltzGen binder design workflow using Protos.

This script demonstrates a real-world workflow:
1. Download a target structure from PDB
2. Register it in StructureProcessor
3. Extract sequence and identify binding site
4. Prepare a BoltzGen job to design a binder
5. Run via ModelManager with GPU
6. Register designed structures back into Protos

Target: MDM2 (1YCR) - classic peptide binder design target
        MDM2 binds p53 transactivation domain
"""
import sys
from pathlib import Path
import pandas as pd

# Ensure protos is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor
from protos.processing.structure.structure_utils import load_structure
from protos.models import ModelManager, JobStatus


def main():
    print("=" * 60)
    print("BoltzGen Binder Design Workflow")
    print("=" * 60)

    # =========================================================================
    # Step 1: Download target structure from PDB
    # =========================================================================
    print("\n[1/6] Downloading target structure...")

    # Using 1YCR: MDM2-p53 peptide complex
    # Chain A = MDM2 (target), Chain B = p53 peptide (we'll design a replacement)
    pdb_id = "1ycr"

    loader = StructureLoader()
    sp = StructureProcessor()

    # Download and register
    registered_name = loader.download_and_register(
        identifier=pdb_id,
        name=pdb_id,
        metadata={"description": "MDM2-p53 peptide complex", "target_chain": "A"}
    )

    if not registered_name:
        print(f"  ✗ Failed to download {pdb_id}")
        return

    print(f"  ✓ Downloaded and registered: {registered_name}")

    # =========================================================================
    # Step 2: Load and analyze the structure
    # =========================================================================
    print("\n[2/6] Analyzing target structure...")

    df = sp.load_entity(registered_name)
    if df is None:
        print(f"  ✗ Failed to load structure")
        return

    # Reset index for easier manipulation
    df = df.reset_index()

    # Get chain info
    chains = df['auth_chain_id'].unique()
    print(f"  Chains found: {list(chains)}")

    # Analyze each chain
    for chain in chains:
        chain_df = df[df['auth_chain_id'] == chain]
        ca_atoms = chain_df[chain_df['atom_name'] == 'CA']
        n_residues = len(ca_atoms)
        seq = ''.join(ca_atoms.drop_duplicates('gen_seq_id').sort_values('gen_seq_id')['res_name1l'].tolist())
        print(f"  Chain {chain}: {n_residues} residues")
        print(f"    Sequence: {seq[:50]}..." if len(seq) > 50 else f"    Sequence: {seq}")

    # MDM2 is chain A (target), p53 peptide is chain B
    target_chain = "A"

    # =========================================================================
    # Step 3: Locate original CIF for BoltzGen
    # =========================================================================
    print("\n[3/6] Preparing target for BoltzGen...")

    job_name = "mdm2_binder_design"

    # Use the ORIGINAL downloaded CIF (has full metadata including _entity_poly_seq)
    # The StructureLoader downloads to: {data_root}/structure/mmcif/{pdb_id}.cif
    original_cif = Path(sp.paths.data_root) / "structure" / "mmcif" / f"{pdb_id}.cif"

    if not original_cif.exists():
        print(f"  ✗ Original CIF not found at {original_cif}")
        return

    print(f"  ✓ Using original CIF: {original_cif}")

    # =========================================================================
    # Step 4: Prepare BoltzGen job via ModelManager
    # =========================================================================
    print("\n[4/6] Preparing BoltzGen job...")

    mm = ModelManager()
    print(f"  GPU enabled: {mm._executor.use_gpu}")

    # BoltzGen config for peptide binder design against MDM2
    # We'll design a peptide (15-25 residues) to bind to MDM2's p53 binding site
    config = {
        "job_name": job_name,
        "entities": [
            # Designed peptide binder
            {
                "protein": {
                    "id": "B",  # Chain B for designed binder
                    "sequence": "15..25",  # Design 15-25 residue peptide
                }
            },
            # Target structure (MDM2) - use original CIF with full metadata
            {
                "file": {
                    "path": str(original_cif),
                    "include": [
                        {"chain": {"id": "A"}}  # Only use chain A (MDM2)
                    ],
                    # Specify binding site (p53 binding groove residues)
                    "binding_types": [
                        {
                            "chain": {
                                "id": "A",
                                "binding": "25..109"  # MDM2 p53-binding domain
                            }
                        }
                    ]
                }
            }
        ],
        "protocol": "peptide-anything",  # Peptide design protocol
        "num_designs": 4,  # Generate 4 designs for testing
    }

    print(f"  Config: peptide binder (15-25 aa) against MDM2 chain A")
    print(f"  Protocol: {config['protocol']}")
    print(f"  Num designs: {config['num_designs']}")

    invocation = mm.prepare("boltzgen", config=config)
    print(f"  ✓ Job prepared")
    print(f"    Working dir: {invocation.job.working_dir}")
    print(f"    Command: {' '.join(invocation.job.command[:10])}...")

    # =========================================================================
    # Step 5: Submit and run the job
    # =========================================================================
    print("\n[5/6] Running BoltzGen (this may take several minutes)...")

    state = mm.submit_job(invocation)
    print(f"  Job ID: {state.job_id}")
    print(f"  Status: {state.status.value}")

    # Wait for completion
    print("  Waiting for completion...")
    final_state = mm.wait_for_job(
        state.job_id,
        timeout_seconds=3600,  # 1 hour max
        poll_interval=15.0
    )

    print(f"\n  Final status: {final_state.status.value}")
    if final_state.result:
        print(f"  Duration: {final_state.result.duration_seconds:.1f}s")
        print(f"  Exit code: {final_state.result.exit_code}")
        print(f"  Output files: {len(final_state.result.output_files)}")

    if final_state.error:
        print(f"  Error: {final_state.error}")

    # =========================================================================
    # Step 6: Register designed structures
    # =========================================================================
    print("\n[6/6] Registering designed structures...")

    working_dir = Path(invocation.job.working_dir)

    # Check multiple possible output locations
    output_dirs = [
        working_dir / "predictions" / "intermediate_designs_inverse_folded" / "refold_cif",
        working_dir / "predictions" / "final_ranked_designs",
        working_dir / "predictions" / "intermediate_designs",
    ]

    designs_registered = 0
    for output_dir in output_dirs:
        if not output_dir.exists():
            continue

        for cif_file in sorted(output_dir.glob("*.cif")):
            structure_id = f"boltzgen_{job_name}_{cif_file.stem}"

            try:
                df = load_structure(cif_file, structure_id=structure_id)
                df_reset = df.reset_index() if 'atom_id' in df.index.names else df

                # BoltzGen outputs: Chain A = designed binder, Chain B = target
                # Extract the designed binder sequence (chain A)
                if 'auth_chain_id' in df_reset.columns:
                    binder_df = df_reset[df_reset['auth_chain_id'] == 'A']
                    target_df = df_reset[df_reset['auth_chain_id'] == 'B']
                else:
                    binder_df = df_reset
                    target_df = pd.DataFrame()

                # Get binder sequence
                binder_ca = binder_df[binder_df['atom_name'] == 'CA']
                if len(binder_ca) > 0:
                    binder_seq = ''.join(
                        binder_ca.drop_duplicates('gen_seq_id')
                        .sort_values('gen_seq_id')['res_name1l'].tolist()
                    )
                else:
                    binder_seq = "unknown"

                # Get target sequence for reference
                target_ca = target_df[target_df['atom_name'] == 'CA'] if len(target_df) > 0 else pd.DataFrame()
                target_seq = ''.join(
                    target_ca.drop_duplicates('gen_seq_id')
                    .sort_values('gen_seq_id')['res_name1l'].tolist()
                ) if len(target_ca) > 0 else ""

                # Register the full complex
                sp.save_entity(
                    structure_id,
                    df,
                    metadata={
                        "source": "boltzgen",
                        "job_id": state.job_id,
                        "job_name": job_name,
                        "target_pdb": pdb_id,
                        "binder_chain": "A",
                        "target_chain": "B",
                        "binder_sequence": binder_seq,
                        "binder_length": len(binder_seq),
                        "target_sequence": target_seq,
                        "target_length": len(target_seq),
                    }
                )
                designs_registered += 1
                print(f"  ✓ {structure_id}")
                print(f"      Binder ({len(binder_seq)} aa): {binder_seq}")
                if target_seq:
                    print(f"      Target ({len(target_seq)} aa): {target_seq[:30]}...")

            except Exception as e:
                print(f"  ✗ {cif_file.name}: {e}")

    print(f"\n  Total designs registered: {designs_registered}")

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("WORKFLOW COMPLETE")
    print("=" * 60)
    print(f"Target: {pdb_id} (MDM2)")
    print(f"Job: {job_name}")
    print(f"Status: {final_state.status.value}")
    print(f"Designs registered: {designs_registered}")

    if designs_registered > 0:
        print("\nTo load a designed structure:")
        print(f"  sp = StructureProcessor()")
        print(f"  df = sp.load_entity('boltzgen_{job_name}_config_0')")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Minimal BoltzGen workflow example using Protos ModelManager.

This demonstrates:
1. Creating a BoltzGen job via ModelManager
2. Running it in Docker with GPU support
3. Registering results into StructureProcessor

Usage:
    python scripts/boltzgen_minimal_example.py
"""
from pathlib import Path
from protos.models import ModelManager, JobStatus
from protos.processing.structure import StructureProcessor
from protos.processing.structure.structure_utils import load_structure


def main():
    # 1. Initialize ModelManager
    print("=== Initializing ModelManager ===")
    mm = ModelManager()
    print(f"GPU enabled: {mm._executor.use_gpu}")

    # 2. Prepare a BoltzGen job
    # De novo protein design (40-60 residues)
    config = {
        "job_name": "minimal_example",
        "entities": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "40..60",  # Design 40-60 residue protein
                }
            }
        ],
        "num_designs": 2,  # Generate 2 designs
        # Optional: skip analysis/filtering for de novo (no target)
        # These steps require a target structure for SASA calculation
    }

    print("\n=== Preparing Job ===")
    invocation = mm.prepare("boltzgen", config=config)
    print(f"Working dir: {invocation.job.working_dir}")
    print(f"Command: {' '.join(invocation.job.command)}")

    # 3. Submit and wait for completion
    print("\n=== Submitting Job ===")
    state = mm.submit_job(invocation)
    print(f"Job ID: {state.job_id}")

    print("\n=== Waiting for completion (this may take several minutes) ===")
    final_state = mm.wait_for_job(state.job_id, timeout_seconds=3600, poll_interval=10.0)

    print(f"\nStatus: {final_state.status.value}")
    if final_state.result:
        print(f"Duration: {final_state.result.duration_seconds:.1f}s")
        print(f"Output files: {len(final_state.result.output_files)}")

    # 4. Register results into StructureProcessor
    # Note: For de novo designs, analysis step may fail (no target for SASA)
    # but design/folding outputs are still valid
    working_dir = Path(invocation.job.working_dir)
    refold_dir = working_dir / "predictions" / "intermediate_designs_inverse_folded" / "refold_cif"

    if refold_dir.exists():
        print("\n=== Registering Structures ===")
        sp = StructureProcessor()

        for cif_file in sorted(refold_dir.glob("*.cif")):
            structure_id = f"boltzgen_{config['job_name']}_{cif_file.stem}"
            try:
                df = load_structure(cif_file, structure_id=structure_id)

                # Extract sequence
                ca_atoms = df[df['atom_name'] == 'CA'].drop_duplicates('gen_seq_id').sort_values('gen_seq_id')
                seq = ''.join(ca_atoms['res_name1l'].tolist())

                # Register
                sp.save_entity(
                    structure_id,
                    df,
                    metadata={
                        "source": "boltzgen",
                        "job_id": state.job_id,
                        "sequence": seq,
                        "sequence_length": len(seq),
                    }
                )
                print(f"  ✓ {structure_id}: {len(seq)} residues, {len(df)} atoms")

            except Exception as e:
                print(f"  ✗ {cif_file.name}: {e}")

        # Verify
        print("\n=== Verification ===")
        for cif_file in sorted(refold_dir.glob("*.cif"))[:2]:
            sid = f"boltzgen_{config['job_name']}_{cif_file.stem}"
            entity = sp.load_entity(sid)
            if entity is not None:
                print(f"Loaded {sid}: {len(entity)} atoms")
    else:
        print(f"\nNo refold outputs found at {refold_dir}")
        print("Check job logs for errors.")

    print("\n=== Done ===")


if __name__ == "__main__":
    main()

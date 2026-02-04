#!/usr/bin/env python3
"""Test RFdiffusion2 job submission for rhodozyme design.

This script tests the RFdiffusion2 integration with a single small run
using placement_00 from the BoltzGen rhodozyme outputs.
"""

from pathlib import Path
from protos.models.model_manager import ModelManager
from protos.models.rfdiffusion2_configs import (
    RFD2ConfigBuilder,
    partial_diffusion_config,
)

# Paths
PLACEMENT_DIR = Path("/data/fast/projects/protos/data/models/boltzgen/20260202_rhodozyme")
OUTPUT_BASE = Path("/data/fast/projects/protos/data/models/rfdiffusion2")

# Triad positions per placement (from config template)
TRIAD_POSITIONS = {
    "placement_00": {"SER": 230, "HIS": 245, "ASP": 250},
    "placement_01": {"SER": 134, "HIS": 138, "ASP": 226},
    "placement_02": {"SER": 230, "HIS": 242, "ASP": 246},
    "placement_03": {"SER": 249, "HIS": 253, "ASP": 311},
    "placement_04": {"SER": 226, "HIS": 230, "ASP": 246},
    "placement_05": {"SER": 229, "HIS": 230, "ASP": 250},
    "placement_06": {"SER": 139, "HIS": 227, "ASP": 230},
    "placement_07": {"SER": 243, "HIS": 247, "ASP": 249},
}


def build_test_config(placement_name: str = "placement_00") -> dict:
    """Build a minimal test config for RFdiffusion2.

    For the test, we do partial diffusion on the full structure
    with the catalytic triad atoms constrained.
    """
    placement_cif = PLACEMENT_DIR / placement_name / f"{placement_name}.cif"
    output_dir = OUTPUT_BASE / f"test_{placement_name}"

    triad = TRIAD_POSITIONS[placement_name]

    # Build atomic constraints for catalytic triad
    atom_constraints = {
        f"A{triad['SER']}": ["OG", "CB", "CA"],           # SER hydroxyl
        f"A{triad['HIS']}": ["NE2", "ND1", "CE1", "CD2", "CG", "CB", "CA"],  # HIS imidazole
        f"A{triad['ASP']}": ["OD1", "OD2", "CG", "CB", "CA"],  # ASP carboxyl
    }

    # Use partial_diffusion_config helper for a simple test
    config = partial_diffusion_config(
        input_pdb=placement_cif,  # Will be converted to PDB by adapter
        output_dir=output_dir,
        atom_constraints=atom_constraints,
        partial_T=20,        # Conservative - less change
        num_designs=2,       # Just 2 designs for test
        ligands=["RET"],     # Include retinal
        chain="A",
        chain_length=326,    # Approximate length of rhodopsin
        job_name=f"test_{placement_name}",
        stop_step="sweep",   # Stop after design, skip LigandMPNN for test
    )

    return config


def main():
    """Run test RFdiffusion2 job."""
    print("=" * 60)
    print("RFdiffusion2 Test Job Submission")
    print("=" * 60)

    # Build test config
    placement = "placement_00"
    config = build_test_config(placement)

    print(f"\nPlacement: {placement}")
    print(f"Input: {config['input_pdb']}")
    print(f"Output: {config['output_dir']}")
    print(f"Num designs: {config['num_designs']}")
    print(f"Partial T: {config['partial_T']}")
    print(f"Triad constraints: {list(config['contig_atoms'].keys())}")

    # Initialize ModelManager
    print("\nInitializing ModelManager...")
    manager = ModelManager()

    # Prepare job
    print("\nPreparing RFdiffusion2 job...")
    try:
        invocation = manager.prepare(
            "rfdiffusion2",
            config=config,
        )
        print(f"Job prepared: {invocation.model}")
        print(f"Has external job: {invocation.is_external()}")

        # Submit job
        print("\nSubmitting job...")
        result = manager.submit_job(invocation)
        print(f"\nJob submitted successfully!")
        print(f"Job ID: {result.job_id}")
        print(f"Status: {result.status}")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

#!/usr/bin/env python3
"""RFdiffusion2 runner for Protos.

Provides a wrapper around RFdiffusion2 for protein backbone generation
with support for partial diffusion and motif scaffolding.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field


# RFdiffusion2 installation path
RFDIFFUSION2_ROOT = Path(__file__).parent / "RFdiffusion2"
SIF_PATH = RFDIFFUSION2_ROOT / "rf_diffusion" / "exec" / "bakerlab_rf_diffusion_aa.sif"


@dataclass
class RFD2Config:
    """Configuration for RFdiffusion2 run."""

    # Input structure (must be PDB format with ORI HETATM for centering)
    input_pdb: Path

    # Output settings
    output_dir: Path
    output_prefix: str = "design"
    num_designs: int = 10

    # Diffusion settings
    diffusion_steps: int = 100  # T parameter

    # Ligand settings (e.g., "RET" for retinal)
    ligand: Optional[str] = None

    # Config name (aa, guidepost, etc.)
    config_name: str = "aa"

    # Additional inference parameters as hydra overrides
    extra_params: Dict[str, Any] = field(default_factory=dict)

    def to_command_args(self) -> List[str]:
        """Convert config to Hydra command line arguments."""
        args = []

        # Input PDB - path will be mapped in container
        args.append(f"inference.input_pdb={self.input_pdb}")

        # Output
        args.append(f"inference.num_designs={self.num_designs}")
        args.append(f"inference.output_prefix={self.output_dir}/{self.output_prefix}")

        # Diffusion steps
        args.append(f"diffuser.T={self.diffusion_steps}")

        # Ligand
        if self.ligand:
            args.append(f"inference.ligand={self.ligand}")

        # Extra params
        for key, value in self.extra_params.items():
            if isinstance(value, bool):
                args.append(f"{key}={str(value)}")
            elif isinstance(value, str):
                args.append(f"{key}={value}")
            else:
                args.append(f"{key}={value}")

        return args


def run_rfdiffusion2(
    config: RFD2Config,
    use_gpu: bool = True,
    verbose: bool = True,
) -> Path:
    """Run RFdiffusion2 with given configuration.

    Args:
        config: RFD2Config with run parameters
        use_gpu: Whether to use GPU (--nv flag)
        verbose: Print command output

    Returns:
        Path to output directory
    """
    if not SIF_PATH.exists():
        raise FileNotFoundError(
            f"RFdiffusion2 SIF not found at {SIF_PATH}. "
            f"Download from Baker Lab or build container."
        )

    if not config.input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {config.input_pdb}")

    # Ensure output directory exists
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Build apptainer command
    cmd = ["apptainer", "exec"]

    if use_gpu:
        cmd.append("--nv")

    # Environment variables
    cmd.extend([
        "--env", "HOME=/tmp",
        "--env", "DGLBACKEND=pytorch",
        "--env", "PYTHONPATH=/app:/app/rf2aa",
    ])

    # Bind mounts
    cmd.extend([
        "--bind", f"{RFDIFFUSION2_ROOT}:/app",
        "--bind", f"{config.input_pdb.parent}:{config.input_pdb.parent}",
        "--bind", f"{config.output_dir}:{config.output_dir}",
    ])

    # Working directory and container
    cmd.extend([
        "--pwd", "/app/rf_diffusion",
        str(SIF_PATH),
    ])

    # Python command
    cmd.extend([
        "python", "run_inference.py",
        f"--config-name={config.config_name}",
    ])

    # Add config arguments
    cmd.extend(config.to_command_args())

    if verbose:
        print(f"Running RFdiffusion2:")
        print(f"  Input: {config.input_pdb}")
        print(f"  Output: {config.output_dir}")
        print(f"  Designs: {config.num_designs}")
        print(f"  Steps: {config.diffusion_steps}")

    # Run
    result = subprocess.run(
        cmd,
        capture_output=not verbose,
        text=True,
    )

    if result.returncode != 0:
        error_msg = result.stderr if result.stderr else "Unknown error"
        raise RuntimeError(f"RFdiffusion2 failed: {error_msg}")

    return config.output_dir


def prepare_input_from_cif(
    cif_path: Path,
    output_dir: Path,
    chains_to_keep: Optional[List[str]] = None,
) -> Path:
    """Convert CIF file to PDB format suitable for RFdiffusion2.

    Args:
        cif_path: Path to input CIF file
        output_dir: Directory for output PDB
        chains_to_keep: Optional list of chain IDs to include

    Returns:
        Path to output PDB file with ORI token
    """
    try:
        from .rfdiffusion2_utils import cif_to_pdb_with_ori
    except ImportError:
        from rfdiffusion2_utils import cif_to_pdb_with_ori

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = output_dir / f"{cif_path.stem}.pdb"
    return cif_to_pdb_with_ori(cif_path, pdb_path, chains_to_keep=chains_to_keep)


def run_rhodozyme_design(
    placement_cif: Path,
    output_dir: Path,
    num_designs: int = 1,
    diffusion_steps: int = 100,
    chains_to_keep: Optional[List[str]] = None,
    verbose: bool = True,
) -> Path:
    """Run RFdiffusion2 on a rhodozyme placement.

    This is a high-level function for rhodozyme design that handles
    CIF->PDB conversion and runs the diffusion.

    Args:
        placement_cif: Path to placement CIF file
        output_dir: Output directory
        num_designs: Number of designs to generate
        diffusion_steps: Number of diffusion steps (100 = full, lower = more like input)
        chains_to_keep: Which chains to include (default: all)
        verbose: Print progress

    Returns:
        Path to output directory containing designs
    """
    output_dir = Path(output_dir)

    # Convert CIF to PDB
    if verbose:
        print(f"Converting {placement_cif.name} to PDB format...")

    pdb_path = prepare_input_from_cif(
        placement_cif,
        output_dir / "input",
        chains_to_keep=chains_to_keep,
    )

    # Create config
    config = RFD2Config(
        input_pdb=pdb_path,
        output_dir=output_dir / "designs",
        output_prefix="rhodozyme",
        num_designs=num_designs,
        diffusion_steps=diffusion_steps,
        config_name="aa",
    )

    # Run
    return run_rfdiffusion2(config, verbose=verbose)


if __name__ == "__main__":
    import sys
    # Add parent directory to path for standalone testing
    sys.path.insert(0, str(Path(__file__).parent))

    if len(sys.argv) > 1:
        placement = Path(sys.argv[1])
    else:
        placement = Path("/data/fast/projects/protos/data/models/boltzgen/20260202_rhodozyme/placement_00/placement_00.cif")

    if not placement.exists():
        print(f"Placement file not found: {placement}")
        sys.exit(1)

    output = Path("/tmp/rfd2_test_runner")

    print(f"Running rhodozyme design test...")
    print(f"  Input: {placement}")
    print(f"  Output: {output}")

    result_dir = run_rhodozyme_design(
        placement_cif=placement,
        output_dir=output,
        num_designs=1,
        diffusion_steps=100,
        chains_to_keep=["A"],  # Just opsin chain for now
    )

    print(f"\nDesigns saved to: {result_dir}")

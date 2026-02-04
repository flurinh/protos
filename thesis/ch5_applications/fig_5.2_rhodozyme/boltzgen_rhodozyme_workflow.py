#!/usr/bin/env python3
"""BoltzGen Rhodozyme Design Workflow - Substrate-Specific Binding Pocket Design.

Uses BoltzGen to design the intracellular domain of rhodopsin to create
catalytic binding pockets for different enzyme substrates.

Design Strategy:
- Use active rhodopsin (3PQR) as the scaffold/target
- Design a peptide binder (the ICL region) that creates a binding pocket
- Include substrate molecules as ligands during design
- BoltzGen will optimize the designed region to accommodate the substrate

Substrates tested:
1. Benzamidine (trypsin substrate) - small, charged
2. Phenylalanine (papain substrate) - amino acid
3. Leucine (subtilisin substrate) - hydrophobic amino acid

This creates light-activated enzyme designs where the ICL conformational
change upon light activation could expose/hide the catalytic pocket.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime

import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.structure.structure_utils import load_structure
from protos.io.ingest.structure_loader import StructureLoader
from protos.models import ModelManager, JobStatus


# =============================================================================
# Configuration
# =============================================================================

# Rhodopsin scaffold - active state with open intracellular domain
RHODOPSIN_PDB = "3pqr"  # Metarhodopsin II (active state)
RHODOPSIN_CHAIN = "A"

# Substrates for different enzyme types
SUBSTRATES = {
    "trypsin": {
        "name": "benzamidine",
        "smiles": "NC(=N)c1ccccc1",
        "description": "Trypsin substrate - small charged molecule",
        "enzyme_type": "serine_protease",
    },
    "papain": {
        "name": "phenylalanine",
        "smiles": "NC(Cc1ccccc1)C(=O)O",
        "description": "Papain substrate - aromatic amino acid",
        "enzyme_type": "cysteine_protease",
    },
    "subtilisin": {
        "name": "leucine",
        "smiles": "CC(C)CC(N)C(=O)O",
        "description": "Subtilisin substrate - hydrophobic amino acid",
        "enzyme_type": "serine_protease",
    },
}

# Design parameters
NUM_DESIGNS_PER_SUBSTRATE = 4  # Number of designs to generate per substrate
DESIGNED_PEPTIDE_LENGTH = "20..40"  # Length range for designed ICL region


def download_and_register_rhodopsin(
    loader: StructureLoader,
    sp: StructureProcessor,
) -> tuple[str, Path]:
    """Download rhodopsin structure and return the original CIF path."""
    print(f"\n[1] Downloading rhodopsin scaffold ({RHODOPSIN_PDB})...")

    registered_name = loader.download_and_register(
        identifier=RHODOPSIN_PDB,
        name=RHODOPSIN_PDB,
        metadata={
            "description": "Metarhodopsin II - active state",
            "use": "rhodozyme_scaffold",
        }
    )

    if not registered_name:
        raise RuntimeError(f"Failed to download {RHODOPSIN_PDB}")

    # Get the original CIF path (has full metadata for BoltzGen)
    original_cif = Path(sp.paths.data_root) / "structure" / "mmcif" / f"{RHODOPSIN_PDB}.cif"

    if not original_cif.exists():
        raise FileNotFoundError(f"Original CIF not found at {original_cif}")

    print(f"  ✓ Downloaded and registered: {registered_name}")
    print(f"  ✓ Original CIF: {original_cif}")

    return registered_name, original_cif


def analyze_rhodopsin_structure(sp: StructureProcessor, name: str) -> Dict[str, Any]:
    """Analyze the rhodopsin structure to understand the scaffold."""
    print(f"\n[2] Analyzing rhodopsin structure...")

    df = sp.load_entity(name)
    if df is None:
        raise RuntimeError(f"Failed to load {name}")

    df = df.reset_index()

    # Get chain info
    chain_df = df[df['auth_chain_id'] == RHODOPSIN_CHAIN]

    # Extract sequence
    ca_atoms = chain_df[chain_df['atom_name'] == 'CA'].drop_duplicates('gen_seq_id').sort_values('gen_seq_id')
    sequence = ''.join(ca_atoms['res_name1l'].tolist())
    n_residues = len(ca_atoms)

    print(f"  Chain {RHODOPSIN_CHAIN}: {n_residues} residues")
    print(f"  Sequence: {sequence[:50]}...")

    # Get all chains for reference
    chains = df['auth_chain_id'].unique()
    print(f"  All chains: {list(chains)}")

    return {
        "name": name,
        "chain": RHODOPSIN_CHAIN,
        "n_residues": n_residues,
        "sequence": sequence,
        "all_chains": list(chains),
    }


def prepare_boltzgen_jobs(
    mm: ModelManager,
    rhodopsin_cif: Path,
    substrates: Dict[str, Dict],
) -> List[Dict[str, Any]]:
    """Prepare BoltzGen jobs for each substrate."""
    print(f"\n[3] Preparing BoltzGen jobs for {len(substrates)} substrates...")

    jobs = []

    for substrate_key, substrate_config in substrates.items():
        job_name = f"rhodozyme_{substrate_key}"

        print(f"\n  Substrate: {substrate_config['name']} ({substrate_key})")
        print(f"    SMILES: {substrate_config['smiles']}")
        print(f"    Description: {substrate_config['description']}")

        # BoltzGen config:
        # - Design a peptide that will form the catalytic binding pocket
        # - Use rhodopsin as the target/scaffold structure
        # - Include the substrate as a ligand to design around
        config = {
            "job_name": job_name,
            "entities": [
                # Designed peptide - this will be the catalytic binding domain
                {
                    "protein": {
                        "id": "B",  # Designed chain
                        "sequence": DESIGNED_PEPTIDE_LENGTH,
                    }
                },
                # Rhodopsin scaffold - target structure
                {
                    "file": {
                        "path": str(rhodopsin_cif),
                        "include": [
                            {"chain": {"id": RHODOPSIN_CHAIN}}
                        ],
                        # The designed peptide should bind to the ICL region
                        # TM5-TM6 cytoplasmic ends move during activation
                        "binding_types": [
                            {
                                "chain": {
                                    "id": RHODOPSIN_CHAIN,
                                    # ICL3 region (between TM5-TM6) - most dynamic
                                    "binding": "225..265"
                                }
                            }
                        ]
                    }
                },
                # Substrate ligand - design the pocket around this
                {
                    "ligand": {
                        "id": "SUB",
                        "smiles": substrate_config["smiles"],
                    }
                },
            ],
            "protocol": "peptide-anything",
            "num_designs": NUM_DESIGNS_PER_SUBSTRATE,
        }

        try:
            invocation = mm.prepare("boltzgen", config=config)

            jobs.append({
                "substrate_key": substrate_key,
                "substrate_config": substrate_config,
                "job_name": job_name,
                "invocation": invocation,
                "config": config,
            })

            print(f"    ✓ Job prepared: {invocation.job.working_dir}")
            print(f"    Command: {' '.join(invocation.job.command[:8])}...")

        except Exception as e:
            print(f"    ✗ Failed to prepare job: {e}")

    return jobs


def run_boltzgen_jobs(
    mm: ModelManager,
    jobs: List[Dict[str, Any]],
    timeout_seconds: int = 3600,
) -> List[Dict[str, Any]]:
    """Submit and run all BoltzGen jobs."""
    print(f"\n[4] Running {len(jobs)} BoltzGen jobs...")

    results = []

    for job_info in jobs:
        substrate_key = job_info["substrate_key"]
        invocation = job_info["invocation"]

        print(f"\n  Running: {job_info['job_name']}...")

        try:
            # Submit job
            state = mm.submit_job(invocation)
            print(f"    Job ID: {state.job_id}")
            print(f"    Status: {state.status.value}")

            # Wait for completion
            print(f"    Waiting for completion...")
            final_state = mm.wait_for_job(
                state.job_id,
                timeout_seconds=timeout_seconds,
                poll_interval=15.0,
            )

            print(f"    Final status: {final_state.status.value}")

            if final_state.result:
                print(f"    Duration: {final_state.result.duration_seconds:.1f}s")
                print(f"    Exit code: {final_state.result.exit_code}")
                print(f"    Output files: {len(final_state.result.output_files)}")

            results.append({
                **job_info,
                "state": final_state,
                "job_id": state.job_id,
                "success": final_state.status == JobStatus.COMPLETED,
            })

        except Exception as e:
            print(f"    ✗ Job failed: {e}")
            results.append({
                **job_info,
                "state": None,
                "job_id": None,
                "success": False,
                "error": str(e),
            })

    return results


def register_designs(
    sp: StructureProcessor,
    results: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Register designed structures into StructureProcessor."""
    print(f"\n[5] Registering designed structures...")

    registered = {}

    for result in results:
        if not result.get("success"):
            print(f"\n  Skipping {result['job_name']} (failed)")
            continue

        substrate_key = result["substrate_key"]
        job_name = result["job_name"]
        invocation = result["invocation"]

        print(f"\n  Processing: {job_name}")

        working_dir = Path(invocation.job.working_dir)

        # Check multiple output locations
        output_dirs = [
            working_dir / "predictions" / "intermediate_designs_inverse_folded" / "refold_cif",
            working_dir / "predictions" / "final_ranked_designs",
            working_dir / "predictions" / "intermediate_designs",
        ]

        designs = []

        for output_dir in output_dirs:
            if not output_dir.exists():
                continue

            for cif_file in sorted(output_dir.glob("*.cif")):
                structure_id = f"{job_name}_{cif_file.stem}"

                try:
                    df = load_structure(cif_file, structure_id=structure_id)
                    df_reset = df.reset_index() if 'atom_id' in df.index.names else df

                    # Extract designed binder (chain A in BoltzGen output)
                    if 'auth_chain_id' in df_reset.columns:
                        binder_df = df_reset[df_reset['auth_chain_id'] == 'A']
                    else:
                        binder_df = df_reset

                    # Get binder sequence
                    binder_ca = binder_df[binder_df['atom_name'] == 'CA']
                    if len(binder_ca) > 0:
                        binder_seq = ''.join(
                            binder_ca.drop_duplicates('gen_seq_id')
                            .sort_values('gen_seq_id')['res_name1l'].tolist()
                        )
                    else:
                        binder_seq = "unknown"

                    # Register
                    sp.save_entity(
                        structure_id,
                        df,
                        metadata={
                            "source": "boltzgen",
                            "workflow": "rhodozyme_design",
                            "substrate": substrate_key,
                            "substrate_name": result["substrate_config"]["name"],
                            "substrate_smiles": result["substrate_config"]["smiles"],
                            "designed_sequence": binder_seq,
                            "designed_length": len(binder_seq),
                            "job_id": result.get("job_id"),
                        }
                    )

                    design_info = {
                        "structure_id": structure_id,
                        "sequence": binder_seq,
                        "length": len(binder_seq),
                        "cif_path": str(cif_file),
                    }
                    designs.append(design_info)

                    print(f"    ✓ {structure_id}")
                    print(f"        Designed ({len(binder_seq)} aa): {binder_seq[:40]}{'...' if len(binder_seq) > 40 else ''}")

                except Exception as e:
                    print(f"    ✗ {cif_file.name}: {e}")

        registered[substrate_key] = designs

    return registered


def save_summary(
    results: List[Dict[str, Any]],
    registered: Dict[str, List[Dict[str, Any]]],
    rhodopsin_info: Dict[str, Any],
) -> Path:
    """Save workflow summary."""
    print(f"\n[6] Saving workflow summary...")

    summary = {
        "workflow": "boltzgen_rhodozyme_design",
        "timestamp": datetime.now().isoformat(),
        "scaffold": {
            "pdb": RHODOPSIN_PDB,
            "chain": RHODOPSIN_CHAIN,
            "n_residues": rhodopsin_info["n_residues"],
        },
        "design_parameters": {
            "peptide_length_range": DESIGNED_PEPTIDE_LENGTH,
            "num_designs_per_substrate": NUM_DESIGNS_PER_SUBSTRATE,
            "protocol": "peptide-anything",
        },
        "substrates": {},
    }

    for substrate_key, substrate_config in SUBSTRATES.items():
        designs = registered.get(substrate_key, [])

        # Find the corresponding result
        result = next((r for r in results if r["substrate_key"] == substrate_key), None)

        summary["substrates"][substrate_key] = {
            "name": substrate_config["name"],
            "smiles": substrate_config["smiles"],
            "enzyme_type": substrate_config["enzyme_type"],
            "job_success": result.get("success", False) if result else False,
            "num_designs": len(designs),
            "designs": [
                {
                    "id": d["structure_id"],
                    "sequence": d["sequence"],
                    "length": d["length"],
                }
                for d in designs
            ],
        }

    summary_path = OUTPUT_DIR / "boltzgen_rhodozyme_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  ✓ Saved: {summary_path}")

    return summary_path


def main() -> int:
    """Run the BoltzGen rhodozyme design workflow."""
    print("=" * 70)
    print("BOLTZGEN RHODOZYME DESIGN WORKFLOW")
    print("Substrate-Specific Catalytic Binding Pocket Design")
    print("=" * 70)
    print()
    print("CONCEPT:")
    print("  Use BoltzGen to design peptides that form catalytic binding pockets")
    print("  around specific enzyme substrates, using rhodopsin as the scaffold.")
    print("  The designed regions correspond to the ICL domain which undergoes")
    print("  conformational changes upon light activation.")
    print()
    print(f"  Scaffold: {RHODOPSIN_PDB} (active rhodopsin)")
    print(f"  Substrates: {', '.join(SUBSTRATES.keys())}")
    print(f"  Designs per substrate: {NUM_DESIGNS_PER_SUBSTRATE}")
    print("=" * 70)

    # Initialize Protos - use project data directory
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    print(f"  Data root: {data_root}")

    # Initialize processors
    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)
    mm = ModelManager()

    print(f"\nGPU enabled: {mm._executor.use_gpu}")

    try:
        # Step 1: Download rhodopsin
        rhodopsin_name, rhodopsin_cif = download_and_register_rhodopsin(loader, sp)

        # Step 2: Analyze structure
        rhodopsin_info = analyze_rhodopsin_structure(sp, rhodopsin_name)

        # Step 3: Prepare BoltzGen jobs
        jobs = prepare_boltzgen_jobs(mm, rhodopsin_cif, SUBSTRATES)

        if not jobs:
            print("\nERROR: No jobs were prepared")
            return 1

        # Step 4: Run jobs
        results = run_boltzgen_jobs(mm, jobs)

        # Step 5: Register designs
        registered = register_designs(sp, results)

        # Step 6: Save summary
        summary_path = save_summary(results, registered, rhodopsin_info)

        # Final summary
        print("\n" + "=" * 70)
        print("WORKFLOW COMPLETE")
        print("=" * 70)

        total_designs = sum(len(designs) for designs in registered.values())
        successful_substrates = sum(1 for r in results if r.get("success"))

        print(f"Scaffold: {RHODOPSIN_PDB}")
        print(f"Substrates processed: {successful_substrates}/{len(SUBSTRATES)}")
        print(f"Total designs: {total_designs}")
        print()

        print("Designs per substrate:")
        for substrate_key, designs in registered.items():
            substrate_name = SUBSTRATES[substrate_key]["name"]
            print(f"  {substrate_key} ({substrate_name}): {len(designs)} designs")
            for d in designs[:2]:  # Show first 2
                print(f"    - {d['structure_id']}: {d['sequence'][:30]}...")

        print()
        print(f"Output directory: {OUTPUT_DIR}")
        print(f"Summary: {summary_path}")

        return 0

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""BoltzGen Inpainting Workflow for Rhodozyme Design.

Uses BoltzGen's inpainting capabilities to design the cytoplasmic loops
of rhodopsin while keeping specific positions fixed for catalytic residues.

Design Strategy:
1. FIXED: TM helix cores (membrane-spanning regions)
2. FIXED: Catalytic Cα positions (from triplet geometry matching)
3. DESIGNABLE: Cytoplasmic loops connecting the helices
4. PRESENT: Substrate at center of helix bundle

BoltzGen Configuration:
- file entity: rhodopsin structure
- design: cytoplasmic regions (58-76, 134-157, 222-256, 306-326)
- not_design: catalytic residue positions (keeps their Cα positions)
- structure_groups visibility=0: hide designable regions for inpainting
- ligand: substrate SMILES

The model will:
1. See the fixed TM cores
2. See the catalytic residue positions (but not their identities)
3. Generate the connecting loops to create a functional binding pocket
"""
from __future__ import annotations

import json
import sys
import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
from dataclasses import dataclass

import numpy as np
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
RHODOPSIN_PDB = "3pqr"
RHODOPSIN_CHAIN = "A"

# Designable cytoplasmic regions (from GRN analysis)
# These are the residue ranges that will be designed/inpainted
DESIGNABLE_REGIONS = [
    (58, 76),    # TM1_cyto + ICL1 + TM2_cyto
    (134, 157),  # TM3_cyto + ICL2 + TM4_cyto
    (222, 256),  # TM5_cyto + ICL3 + TM6_cyto (most dynamic!)
    (306, 326),  # TM7_cyto + H8 (excluding gap at 322)
]

# Alternative: use residue list (for structures with gaps)
DESIGNABLE_RESIDUES = "58..76,134..157,222..256,306..321,323..326"

# Reference enzymes with catalytic triads
REFERENCE_ENZYMES = {
    "trypsin": {
        "id": "TRP",
        "pdb": "1S81",
        "chain": "A",
        "description": "Bovine trypsin (serine protease)",
        "triad": {
            "nucleophile": {"res_num": 195, "res_name": "SER"},
            "base": {"res_num": 57, "res_name": "HIS"},
            "electrostatic": {"res_num": 102, "res_name": "ASP"},
        },
        "substrate_smiles": "NC(=N)c1ccccc1",  # Benzamidine
        "substrate_name": "benzamidine",
    },
    "papain": {
        "id": "PAP",
        "pdb": "9PAP",
        "chain": "A",
        "description": "Papain (cysteine protease)",
        "triad": {
            "nucleophile": {"res_num": 25, "res_name": "CYS"},
            "base": {"res_num": 159, "res_name": "HIS"},
            "electrostatic": {"res_num": 175, "res_name": "ASN"},
        },
        "substrate_smiles": "NC(Cc1ccccc1)C(=O)O",  # Phenylalanine
        "substrate_name": "phenylalanine",
    },
}

# Design parameters
NUM_DESIGNS = 4


@dataclass
class CatalyticSite:
    """Configuration for a catalytic site placement."""
    enzyme_key: str
    enzyme_id: str
    positions: Tuple[int, int, int]  # (nucleophile_pos, base_pos, electrostatic_pos)
    grn_positions: Tuple[str, str, str]  # GRN annotations
    score: float
    has_dynamic_helix: bool
    substrate_smiles: str
    substrate_name: str


def format_residue_ranges(ranges: List[Tuple[int, int]]) -> str:
    """Convert list of (start, end) tuples to BoltzGen range string."""
    parts = []
    for start, end in ranges:
        parts.append(f"{start}..{end}")
    return ",".join(parts)


def create_inpainting_yaml_config(
    structure_path: Path,
    chain_id: str,
    designable_regions: str,
    fixed_catalytic_positions: List[int],
    substrate_smiles: str,
    job_name: str,
    visibility_mode: str = "hide_designable",
    secondary_structure: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Create a BoltzGen YAML configuration for inpainting.

    Args:
        structure_path: Path to the rhodopsin CIF file
        chain_id: Chain ID in the structure
        designable_regions: Comma-separated residue ranges for design
        fixed_catalytic_positions: Residue numbers to keep fixed (Cα positions)
        substrate_smiles: SMILES string for the substrate
        job_name: Name for the job
        visibility_mode: How to set structure visibility
            - "hide_designable": Hide designable regions (standard inpainting)
            - "show_all": Show all regions (design with context)
            - "hide_loops_only": Hide only the connecting loops
        secondary_structure: Optional SS conditioning (e.g., {"helix": "58..65", "loop": "66..76"})

    Returns:
        Dict that can be converted to YAML
    """
    # Build the file entity with design/not_design specifications
    file_entity = {
        "file": {
            "path": str(structure_path),
            "include": [
                {"chain": {"id": chain_id}}
            ],
            # Mark cytoplasmic regions as designable
            "design": [
                {"chain": {"id": chain_id, "res_index": designable_regions}}
            ],
        }
    }

    # Add not_design for fixed catalytic positions
    if fixed_catalytic_positions:
        fixed_str = ",".join(str(p) for p in fixed_catalytic_positions)
        file_entity["file"]["not_design"] = [
            {"chain": {"id": chain_id, "res_index": fixed_str}}
        ]

    # Add structure_groups for visibility control
    if visibility_mode == "hide_designable":
        # Standard inpainting: hide designable regions, show scaffold
        file_entity["file"]["structure_groups"] = [
            # Show all by default
            {"group": {"id": "all", "visibility": 1}},
            # Hide designable regions
            {"group": {"id": chain_id, "res_index": designable_regions, "visibility": 0}},
        ]
    elif visibility_mode == "show_all":
        # Show everything (design with full context)
        file_entity["file"]["structure_groups"] = [
            {"group": {"id": "all", "visibility": 1}},
        ]
    elif visibility_mode == "hide_loops_only":
        # Only hide the ICL loops, show helix extensions
        file_entity["file"]["structure_groups"] = [
            {"group": {"id": "all", "visibility": 1}},
            # Hide only ICL regions
            {"group": {"id": chain_id, "res_index": "66..69,142..148,238..240", "visibility": 0}},
        ]

    # Add secondary structure conditioning if provided
    if secondary_structure:
        ss_list = []
        for ss_type, res_index in secondary_structure.items():
            ss_list.append({"chain": {"id": chain_id, ss_type: res_index}})
        file_entity["file"]["secondary_structure"] = ss_list

    # Add binding site specification (substrate should bind near catalytic residues)
    if fixed_catalytic_positions:
        binding_str = ",".join(str(p) for p in fixed_catalytic_positions)
        # Expand binding region around catalytic residues
        file_entity["file"]["binding_types"] = [
            {"chain": {"id": chain_id, "binding": binding_str}}
        ]

    # Build full config
    config = {
        "job_name": job_name,
        "entities": [
            file_entity,
            # Substrate ligand
            {
                "ligand": {
                    "id": "SUB",
                    "smiles": substrate_smiles,
                }
            },
        ],
        "num_designs": NUM_DESIGNS,
    }

    return config


def create_grid_search_configs(
    structure_path: Path,
    chain_id: str,
    catalytic_sites: List[CatalyticSite],
    base_output_dir: Path,
) -> List[Tuple[str, Dict[str, Any], Path]]:
    """Create grid of BoltzGen configurations for all design variants.

    Grid dimensions:
    1. Catalytic site placements (from triplet matching)
    2. Visibility modes
    3. Secondary structure conditioning (optional)

    Returns:
        List of (job_name, config_dict, yaml_path) tuples
    """
    configs = []

    # Visibility modes to try
    visibility_modes = [
        ("inpaint", "hide_designable"),
        ("context", "show_all"),
        ("loopsonly", "hide_loops_only"),
    ]

    # Secondary structure conditioning variants
    ss_variants = [
        ("noss", None),
        ("helixext", {
            # Extend helices slightly into cytoplasmic region
            "helix": "58..62,70..76,134..138,153..157,222..228,250..256,306..309"
        }),
    ]

    for site in catalytic_sites:
        for vis_name, vis_mode in visibility_modes:
            for ss_name, ss_config in ss_variants:
                job_name = f"rhodo_{site.enzyme_id}_{site.positions[0]}_{site.positions[1]}_{site.positions[2]}_{vis_name}_{ss_name}"

                config = create_inpainting_yaml_config(
                    structure_path=structure_path,
                    chain_id=chain_id,
                    designable_regions=DESIGNABLE_RESIDUES,
                    fixed_catalytic_positions=list(site.positions),
                    substrate_smiles=site.substrate_smiles,
                    job_name=job_name,
                    visibility_mode=vis_mode,
                    secondary_structure=ss_config,
                )

                # Save YAML config
                yaml_dir = base_output_dir / "configs"
                yaml_dir.mkdir(parents=True, exist_ok=True)
                yaml_path = yaml_dir / f"{job_name}.yaml"

                with open(yaml_path, "w") as f:
                    yaml.dump({"entities": config["entities"]}, f, default_flow_style=False)

                configs.append((job_name, config, yaml_path))

    return configs


def find_catalytic_triplets(
    rhodopsin_df: pd.DataFrame,
    grn_df: pd.DataFrame,
    enzyme_key: str,
    enzyme_config: Dict,
    top_n: int = 3,
) -> List[CatalyticSite]:
    """Find matching triplet positions in rhodopsin for a given enzyme.

    This is a simplified version - in practice, use the full geometry matching
    from rhodozyme_design_workflow.py.
    """
    # For now, return some example positions based on manual analysis
    # In a full implementation, this would call find_top_triplets()

    enzyme_id = enzyme_config["id"]
    substrate_smiles = enzyme_config["substrate_smiles"]
    substrate_name = enzyme_config["substrate_name"]

    # Example triplet placements (these should come from geometry matching)
    example_sites = [
        # Site 1: ICL2-TM5-TM6 configuration
        CatalyticSite(
            enzyme_key=enzyme_key,
            enzyme_id=enzyme_id,
            positions=(143, 235, 250),  # nucleophile, base, electrostatic
            grn_positions=("34.51", "5.68", "6.33"),
            score=0.85,
            has_dynamic_helix=True,
            substrate_smiles=substrate_smiles,
            substrate_name=substrate_name,
        ),
        # Site 2: TM5-TM6-H8 configuration
        CatalyticSite(
            enzyme_key=enzyme_key,
            enzyme_id=enzyme_id,
            positions=(229, 248, 311),
            grn_positions=("5.62", "6.31", "8.50"),
            score=0.72,
            has_dynamic_helix=True,
            substrate_smiles=substrate_smiles,
            substrate_name=substrate_name,
        ),
        # Site 3: ICL1-TM5-TM6 configuration
        CatalyticSite(
            enzyme_key=enzyme_key,
            enzyme_id=enzyme_id,
            positions=(67, 234, 245),
            grn_positions=("12.49", "5.67", "6.28"),
            score=0.68,
            has_dynamic_helix=True,
            substrate_smiles=substrate_smiles,
            substrate_name=substrate_name,
        ),
    ]

    return example_sites[:top_n]


def download_rhodopsin(
    loader: StructureLoader,
    sp: StructureProcessor,
) -> Tuple[str, Path, pd.DataFrame]:
    """Download rhodopsin and return structure info."""
    print(f"\n[1] Downloading rhodopsin scaffold ({RHODOPSIN_PDB})...")

    registered_name = loader.download_and_register(
        identifier=RHODOPSIN_PDB,
        name=RHODOPSIN_PDB,
        metadata={
            "description": "Metarhodopsin II - active state",
            "use": "rhodozyme_inpainting_scaffold",
        }
    )

    if not registered_name:
        raise RuntimeError(f"Failed to download {RHODOPSIN_PDB}")

    # Get paths
    original_cif = Path(sp.paths.data_root) / "structure" / "mmcif" / f"{RHODOPSIN_PDB}.cif"

    if not original_cif.exists():
        raise FileNotFoundError(f"Original CIF not found at {original_cif}")

    # Load structure
    df = sp.load_entity(registered_name)
    if df is None:
        raise RuntimeError(f"Failed to load structure {registered_name}")

    print(f"  ✓ Downloaded: {registered_name}")
    print(f"  ✓ CIF path: {original_cif}")

    return registered_name, original_cif, df.reset_index()


def main() -> int:
    """Run the BoltzGen inpainting workflow for rhodozyme design."""
    print("=" * 70)
    print("BOLTZGEN INPAINTING WORKFLOW - RHODOZYME DESIGN")
    print("Constrained Binding Pocket Generation")
    print("=" * 70)
    print()
    print("DESIGN STRATEGY:")
    print("  1. FIXED: TM helix cores (membrane-spanning)")
    print("  2. FIXED: Catalytic Cα positions (from triplet matching)")
    print("  3. DESIGN: Cytoplasmic loops connecting helices")
    print("  4. PRESENT: Substrate at binding site")
    print()
    print("DESIGNABLE REGIONS:")
    for start, end in DESIGNABLE_REGIONS:
        print(f"    {start}-{end}")
    print()
    print(f"  Reference enzymes: {list(REFERENCE_ENZYMES.keys())}")
    print("=" * 70)

    # Initialize Protos - use project data directory
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    print(f"  Data root: {data_root}")

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)
    mm = ModelManager()

    print(f"\nGPU enabled: {mm._executor.use_gpu}")

    try:
        # Step 1: Download rhodopsin
        rhodopsin_name, rhodopsin_cif, rhodopsin_df = download_rhodopsin(loader, sp)

        # Step 2: Annotate with GRN
        print(f"\n[2] Annotating with GRN...")
        grn_df = sp.annotate_with_grn(rhodopsin_name, chains=[RHODOPSIN_CHAIN])

        if grn_df is not None:
            print(f"  ✓ Annotated {len(grn_df)} atoms with GRN")
        else:
            print("  ⚠ GRN annotation not available")
            grn_df = pd.DataFrame()

        # Step 3: Find catalytic triplet placements for each enzyme
        print(f"\n[3] Finding catalytic triplet placements...")

        all_sites = []
        for enzyme_key, enzyme_config in REFERENCE_ENZYMES.items():
            print(f"\n  {enzyme_config['description']}:")
            sites = find_catalytic_triplets(
                rhodopsin_df, grn_df, enzyme_key, enzyme_config, top_n=3
            )
            for site in sites:
                print(f"    {site.positions} | score={site.score:.2f} | dynamic={site.has_dynamic_helix}")
            all_sites.extend(sites)

        print(f"\n  Total catalytic site configurations: {len(all_sites)}")

        # Step 4: Create grid search configurations
        print(f"\n[4] Creating BoltzGen configurations...")

        configs = create_grid_search_configs(
            structure_path=rhodopsin_cif,
            chain_id=RHODOPSIN_CHAIN,
            catalytic_sites=all_sites,
            base_output_dir=OUTPUT_DIR,
        )

        print(f"  ✓ Created {len(configs)} configurations")
        print(f"  ✓ Config files saved to: {OUTPUT_DIR}/configs/")

        # Show example config
        if configs:
            print(f"\n  Example configuration ({configs[0][0]}):")
            example_config = configs[0][1]
            print(f"    Entities: {len(example_config['entities'])}")

            file_entity = example_config["entities"][0]["file"]
            print(f"    Design regions: {file_entity.get('design', [{}])[0].get('chain', {}).get('res_index', 'N/A')}")
            print(f"    Not design: {file_entity.get('not_design', [{}])[0].get('chain', {}).get('res_index', 'N/A')}")

            # Show full YAML
            print(f"\n  Full YAML for example config:")
            print("  " + "-" * 60)
            yaml_content = yaml.dump({"entities": example_config["entities"]}, default_flow_style=False)
            for line in yaml_content.split("\n"):
                print(f"  {line}")
            print("  " + "-" * 60)

        # Step 5: Prepare and run jobs (limited for testing)
        print(f"\n[5] Preparing BoltzGen jobs...")

        # For testing, only run first config
        test_configs = configs[:1]

        jobs_prepared = []
        for job_name, config, yaml_path in test_configs:
            try:
                invocation = mm.prepare("boltzgen", config=config)
                jobs_prepared.append({
                    "job_name": job_name,
                    "config": config,
                    "yaml_path": yaml_path,
                    "invocation": invocation,
                })
                print(f"  ✓ {job_name}")
                print(f"      Working dir: {invocation.job.working_dir}")
            except Exception as e:
                print(f"  ✗ {job_name}: {e}")

        # Save summary of all configurations
        summary = {
            "workflow": "boltzgen_inpainting_rhodozyme",
            "timestamp": datetime.now().isoformat(),
            "scaffold": {
                "pdb": RHODOPSIN_PDB,
                "chain": RHODOPSIN_CHAIN,
            },
            "designable_regions": DESIGNABLE_RESIDUES,
            "num_configs": len(configs),
            "catalytic_sites": [
                {
                    "enzyme": site.enzyme_key,
                    "positions": site.positions,
                    "grn": site.grn_positions,
                    "score": site.score,
                    "substrate": site.substrate_name,
                }
                for site in all_sites
            ],
            "configs": [
                {
                    "job_name": name,
                    "yaml_path": str(yaml_path),
                }
                for name, _, yaml_path in configs
            ],
        }

        summary_path = OUTPUT_DIR / "inpainting_workflow_summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        print(f"\n  ✓ Saved summary: {summary_path}")

        # Final output
        print("\n" + "=" * 70)
        print("WORKFLOW PREPARED")
        print("=" * 70)
        print(f"Total configurations: {len(configs)}")
        print(f"Jobs prepared: {len(jobs_prepared)}")
        print()
        print("GRID DIMENSIONS:")
        print(f"  Catalytic sites: {len(all_sites)}")
        print(f"  Visibility modes: 3 (inpaint, context, loopsonly)")
        print(f"  SS conditioning: 2 (none, helix-extended)")
        print()
        print("NEXT STEPS:")
        print("  1. Review YAML configs in: {OUTPUT_DIR}/configs/")
        print("  2. Run full grid search with: mm.submit_job()")
        print("  3. Compare results by catalytic geometry preservation")
        print()
        print(f"Output: {OUTPUT_DIR}")

        return 0

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""GPCR Binding Pocket Workflow: Comprehensive Mechanism Analysis.

This workflow demonstrates ProtOS capabilities for GPCR structure analysis:
- Loading and annotating structures with GRN (Generic Residue Numbers)
- Ligand interaction analysis using ProtOS utilities
- Hypothesis testing for activation mechanisms
- PyMOL visualization generation

Structures analyzed (8 adrenergic receptors):
  ADRB2 (full agonists, active):
    - 3SN6: BI-167107 + Gs protein
    - 4LDO: Adrenaline + Nb6B9 nanobody
  ADRB1 (full agonist, active-like):
    - 2Y02: Isoprenaline
  ADRB1 (partial agonists, intermediate):
    - 2Y04: Salbutamol
    - 2Y00: Dobutamine
  ADRB2 (inverse agonists, inactive):
    - 2RH1: Carazolol
    - 3NY9: ICI 118,551
  ADRB1 (antagonist, inactive):
    - 2VT4: Cyanopindolol

Hypotheses tested:
  H1: Agonists bind CLOSER to S5.43 than inverse agonists
  H2: Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists
  H3: Water at N6.55 is EXCLUSIVE to agonist-bound active structures

Processors used: Structure -> GRN annotation -> Ligand Analysis -> Property
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from collections import defaultdict

import pandas as pd
import numpy as np

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
SRC_DIR = REPO_ROOT / "src"
LOG_DIR = THESIS_DIR / "logs"
OUTPUT_DIR = THESIS_DIR / "outputs" / "binding_pocket"
FIGURES_DIR = THESIS_DIR / "workflows" / "figures"

if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

LOG_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Configure logging
log_file = LOG_DIR / f"gpcr_binding_pocket_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.structure.ligand_interactions import LigandInteractionAnalyzer
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Structure Configuration
# =============================================================================
STRUCTURES: Dict[str, Dict[str, Any]] = {
    # ADRB2 - Full Agonists (Active)
    "3SN6": {
        "chain": "R",
        "ligand": "P0G",
        "ligand_name": "BI-167107",
        "ligand_type": "full_agonist",
        "state": "active",
        "receptor": "ADRB2",
        "color": "forest",
    },
    "4LDO": {
        "chain": "A",
        "ligand": "ALE",
        "ligand_name": "Adrenaline",
        "ligand_type": "full_agonist",
        "state": "active",
        "receptor": "ADRB2",
        "color": "lime",
    },
    # ADRB1 - Full Agonist (Active-like)
    "2Y02": {
        "chain": "A",
        "ligand": "WHJ",
        "ligand_name": "Isoprenaline",
        "ligand_type": "full_agonist",
        "state": "active_like",
        "receptor": "ADRB1",
        "color": "cyan",
    },
    # ADRB1 - Partial Agonists (Intermediate)
    "2Y04": {
        "chain": "A",
        "ligand": "68H",
        "ligand_name": "Salbutamol",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
        "color": "orange",
    },
    "2Y00": {
        "chain": "A",
        "ligand": "Y00",
        "ligand_name": "Dobutamine",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
        "color": "tv_orange",
    },
    # ADRB2 - Inverse Agonists (Inactive)
    "2RH1": {
        "chain": "A",
        "ligand": "CAU",
        "ligand_name": "Carazolol",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
        "color": "firebrick",
    },
    "3NY9": {
        "chain": "A",
        "ligand": "JSZ",
        "ligand_name": "ICI 118,551",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
        "color": "salmon",
    },
    # ADRB1 - Antagonist (Inactive)
    "2VT4": {
        "chain": "A",
        "ligand": "P32",
        "ligand_name": "Cyanopindolol",
        "ligand_type": "antagonist",
        "state": "inactive",
        "receptor": "ADRB1",
        "color": "slate",
    },
}

# Hypothesis-specific GRN positions
H1_GRN = "5.43"  # Agonists bind closer to S5.43
H2_GRN = "6.48"  # Inverse agonists bind closer to W6.48 (toggle switch)
H3_GRN = "6.55"  # Water at N6.55 exclusive to active structures

BINDING_CUTOFF = 5.0
WATER_CUTOFF = 3.5


# =============================================================================
# Step 1: Load Structures using ProtOS
# =============================================================================

def load_structures(struct_proc: StructureProcessor) -> Dict[str, pd.DataFrame]:
    """Load all 8 GPCR structures using ProtOS StructureProcessor."""
    logger.info("=" * 60)
    logger.info("Step 1: Loading structures via ProtOS")
    logger.info("=" * 60)

    loader = StructureLoader(processor=struct_proc)
    structures = {}

    for pdb_id, info in STRUCTURES.items():
        logger.info(f"  {pdb_id}: {info['receptor']} + {info['ligand_name']} ({info['ligand_type']})")

        try:
            # Try cache first
            cached = struct_proc.load_entity(pdb_id.lower())
            if cached is not None and isinstance(cached, pd.DataFrame) and not cached.empty:
                structures[pdb_id] = cached
                logger.info(f"    Loaded from cache: {len(cached)} atoms")
                continue

            # Download via ProtOS
            result = loader.download_and_register(
                pdb_id,
                name=pdb_id.lower(),
                materialize_entities=True,
            )

            if result:
                df = struct_proc.load_entity(pdb_id.lower())
                if df is not None and isinstance(df, pd.DataFrame) and not df.empty:
                    structures[pdb_id] = df
                    logger.info(f"    Downloaded: {len(df)} atoms")

        except Exception as e:
            logger.warning(f"    Failed: {e}")

    logger.info(f"\nLoaded {len(structures)}/8 structures")
    return structures


# =============================================================================
# Step 2: Annotate Structures with GRN using ProtOS
# =============================================================================

def annotate_structures_with_grn(
    struct_proc: StructureProcessor,
    structures: Dict[str, pd.DataFrame],
) -> Dict[str, pd.DataFrame]:
    """Annotate each structure with GRN using ProtOS annotate_with_grn."""
    logger.info("=" * 60)
    logger.info("Step 2: GRN annotation via ProtOS")
    logger.info("=" * 60)

    annotated = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        chain = info["chain"]

        logger.info(f"  {pdb_id} chain {chain}...")

        try:
            # Store in processor for annotation
            struct_proc.frames[pdb_id.lower()] = df

            # Use ProtOS annotate_with_grn - handles sequence extraction,
            # GRN alignment, and mapping back to structure residues
            annotated_df = struct_proc.annotate_with_grn(
                pdb_id.lower(),
                reference_table="gpcrdb_ref",
                protein_family="gpcr_a",
                chains=[chain],
                save=False,  # Don't overwrite original
            )

            annotated[pdb_id] = annotated_df

            # Count GRN annotations
            df_reset = annotated_df.reset_index() if hasattr(annotated_df.index, 'names') else annotated_df
            grn_count = (df_reset['grn'] != '').sum() if 'grn' in df_reset.columns else 0
            logger.info(f"    Annotated {grn_count} atoms with GRN")

        except Exception as e:
            logger.warning(f"    GRN annotation failed: {e}")
            annotated[pdb_id] = df  # Keep unannotated

    return annotated


# =============================================================================
# Step 3: Ligand Interaction Analysis using ProtOS
# =============================================================================

def analyze_ligand_interactions(
    annotated_structures: Dict[str, pd.DataFrame],
) -> Dict[str, Dict[str, Any]]:
    """Analyze ligand-protein interactions using ProtOS LigandInteractionAnalyzer."""
    logger.info("=" * 60)
    logger.info("Step 3: Ligand interaction analysis via ProtOS")
    logger.info("=" * 60)

    results = {}

    for pdb_id, df in annotated_structures.items():
        info = STRUCTURES[pdb_id]
        target_ligand = info["ligand"]
        chain = info["chain"]

        logger.info(f"  {pdb_id} - {info['ligand_name']}...")

        try:
            df_reset = df.reset_index() if hasattr(df.index, 'names') else df

            # Use ProtOS LigandInteractionAnalyzer
            analyzer = LigandInteractionAnalyzer(df_reset)

            # Extract the target ligand
            ligands = analyzer.extract_ligands(
                exclude_common=True,
                allowed_res_names={target_ligand},
            )

            if not ligands:
                logger.warning(f"    Ligand {target_ligand} not found")
                continue

            ligand_info = ligands[0]
            ligand_atoms = ligand_info["atoms"]

            # Get binding site residues via ProtOS
            binding_residues = analyzer.get_binding_site_residues(
                ligand_atoms,
                cutoff=BINDING_CUTOFF,
            )

            if binding_residues.empty:
                logger.warning(f"    No binding residues found")
                continue

            # Map binding residues to GRN using the annotated structure
            grn_mapping = {}
            if 'grn' in df_reset.columns:
                chain_atoms = df_reset[df_reset['auth_chain_id'] == chain]
                for _, row in binding_residues.iterrows():
                    res_id = row['res_id']
                    res_atoms = chain_atoms[chain_atoms['auth_seq_id'] == res_id]
                    if not res_atoms.empty:
                        grn = res_atoms['grn'].iloc[0]
                        if grn and grn != '':
                            grn_mapping[res_id] = grn

            binding_residues['grn'] = binding_residues['res_id'].map(grn_mapping)

            logger.info(f"    Ligand: {ligand_info['num_atoms']} atoms")
            logger.info(f"    Binding residues: {len(binding_residues)}")
            logger.info(f"    GRN-annotated: {sum(binding_residues['grn'].notna())}")

            results[pdb_id] = {
                "ligand": ligand_info,
                "binding_residues": binding_residues,
                "grn_mapping": grn_mapping,
                "annotated_structure": df_reset,
                "chain": chain,
            }

        except Exception as e:
            logger.warning(f"    Analysis failed: {e}")
            import traceback
            traceback.print_exc()

    return results


# =============================================================================
# Step 4: Build GRN to auth_seq_id mapping for hypothesis testing
# =============================================================================

def build_grn_residue_mapping(
    interaction_results: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, int]]:
    """Build mapping from GRN position to auth_seq_id for each structure."""
    logger.info("=" * 60)
    logger.info("Step 4: Building GRN -> residue mappings")
    logger.info("=" * 60)

    grn_mappings = {}

    for pdb_id, data in interaction_results.items():
        df = data["annotated_structure"]
        chain = data["chain"]

        chain_atoms = df[df['auth_chain_id'] == chain]
        if 'grn' not in chain_atoms.columns:
            logger.warning(f"  {pdb_id}: No GRN column")
            continue

        # Build GRN -> auth_seq_id mapping (take first occurrence per residue)
        mapping = {}
        for _, atom in chain_atoms.iterrows():
            grn = atom.get('grn', '')
            if grn and grn != '':
                res_id = atom['auth_seq_id']
                if grn not in mapping:
                    mapping[grn] = int(res_id)

        grn_mappings[pdb_id] = mapping

        # Log key GRN positions
        key_grns = ["5.43", "6.48", "6.55", "2.50", "3.32"]
        found = [f"{g}={mapping.get(g, '?')}" for g in key_grns if g in mapping]
        logger.info(f"  {pdb_id}: {len(mapping)} GRNs ({', '.join(found[:4])})")

    return grn_mappings


# =============================================================================
# Step 5: Hypothesis Testing
# =============================================================================

def calculate_ligand_to_grn_distance(
    df: pd.DataFrame,
    chain: str,
    ligand_code: str,
    grn_pos: str,
    grn_to_resid: Dict[str, int],
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    """Calculate minimum distance from ligand to a GRN-annotated residue."""
    # Get ligand atoms (non-hydrogen)
    ligand_atoms = df[(df['res_name3l'] == ligand_code) & (~df['element'].isin(["H"]))]
    if ligand_atoms.empty:
        return None, None, None

    # Get residue ID for target GRN
    target_res_id = grn_to_resid.get(grn_pos)
    if target_res_id is None:
        return None, None, None

    # Get residue atoms
    res_atoms = df[
        (df['auth_chain_id'] == chain) &
        (df['auth_seq_id'] == target_res_id) &
        (~df['element'].isin(["H"]))
    ]
    if res_atoms.empty:
        return None, None, None

    # Calculate distances
    ligand_coords = ligand_atoms[["x", "y", "z"]].values
    ligand_names = ligand_atoms["atom_name"].values
    res_coords = res_atoms[["x", "y", "z"]].values
    res_names = res_atoms["atom_name"].values

    min_dist = float("inf")
    closest_lig = None
    closest_res = None

    for i, l_coord in enumerate(ligand_coords):
        for j, r_coord in enumerate(res_coords):
            dist = np.linalg.norm(l_coord - r_coord)
            if dist < min_dist:
                min_dist = dist
                closest_lig = ligand_names[i]
                closest_res = res_names[j]

    return min_dist, closest_lig, closest_res


def test_hypothesis_1(
    interaction_results: Dict[str, Dict[str, Any]],
    grn_mappings: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, Any]]:
    """H1: Agonists bind CLOSER to S5.43 than inverse agonists."""
    logger.info("=" * 60)
    logger.info("H1: Ligand distance to S5.43 (agonist binding marker)")
    logger.info("=" * 60)

    results = {}

    for pdb_id, data in interaction_results.items():
        info = STRUCTURES[pdb_id]
        df = data["annotated_structure"]
        chain = info["chain"]
        ligand = info["ligand"]
        grn_to_resid = grn_mappings.get(pdb_id, {})

        dist, lig_atom, res_atom = calculate_ligand_to_grn_distance(
            df, chain, ligand, H1_GRN, grn_to_resid
        )

        if dist is not None:
            results[pdb_id] = {
                "ligand_type": info["ligand_type"],
                "state": info["state"],
                "receptor": info["receptor"],
                "distance": dist,
                "ligand_atom": lig_atom,
                "residue_atom": res_atom,
                "grn_resid": grn_to_resid.get(H1_GRN),
            }
            logger.info(f"  {pdb_id} ({info['ligand_type']:15s}): {dist:.2f} A (resi {grn_to_resid.get(H1_GRN)})")
        else:
            logger.warning(f"  {pdb_id}: Could not calculate distance")

    return results


def test_hypothesis_2(
    interaction_results: Dict[str, Dict[str, Any]],
    grn_mappings: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, Any]]:
    """H2: Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists."""
    logger.info("=" * 60)
    logger.info("H2: Ligand distance to W6.48 (toggle switch)")
    logger.info("=" * 60)

    results = {}

    for pdb_id, data in interaction_results.items():
        info = STRUCTURES[pdb_id]
        df = data["annotated_structure"]
        chain = info["chain"]
        ligand = info["ligand"]
        grn_to_resid = grn_mappings.get(pdb_id, {})

        dist, lig_atom, res_atom = calculate_ligand_to_grn_distance(
            df, chain, ligand, H2_GRN, grn_to_resid
        )

        if dist is not None:
            results[pdb_id] = {
                "ligand_type": info["ligand_type"],
                "state": info["state"],
                "receptor": info["receptor"],
                "distance": dist,
                "ligand_atom": lig_atom,
                "residue_atom": res_atom,
                "grn_resid": grn_to_resid.get(H2_GRN),
            }
            logger.info(f"  {pdb_id} ({info['ligand_type']:15s}): {dist:.2f} A (resi {grn_to_resid.get(H2_GRN)})")
        else:
            logger.warning(f"  {pdb_id}: Could not calculate distance")

    return results


def test_hypothesis_3(
    interaction_results: Dict[str, Dict[str, Any]],
    grn_mappings: Dict[str, Dict[str, int]],
) -> Dict[str, Dict[str, Any]]:
    """H3: Water at N6.55 is EXCLUSIVE to agonist-bound active structures."""
    logger.info("=" * 60)
    logger.info("H3: Water at N6.55 (activation marker)")
    logger.info("=" * 60)
    logger.info("  NOTE: Crystallographic waters depend on resolution")

    results = {}

    for pdb_id, data in interaction_results.items():
        info = STRUCTURES[pdb_id]
        df = data["annotated_structure"]
        chain = info["chain"]
        grn_to_resid = grn_mappings.get(pdb_id, {})

        # Find residue with GRN 6.55
        target_res_id = grn_to_resid.get(H3_GRN)
        if target_res_id is None:
            logger.warning(f"  {pdb_id}: GRN {H3_GRN} not found")
            continue

        # Get N6.55 atoms
        n655_atoms = df[
            (df['auth_chain_id'] == chain) &
            (df['auth_seq_id'] == target_res_id) &
            (~df['element'].isin(["H"]))
        ]

        if n655_atoms.empty:
            logger.warning(f"  {pdb_id}: No atoms for residue {target_res_id}")
            continue

        # Get water oxygens
        waters = df[df['res_name3l'] == "HOH"]

        if waters.empty:
            results[pdb_id] = {
                "ligand_type": info["ligand_type"],
                "state": info["state"],
                "receptor": info["receptor"],
                "has_water": False,
                "min_distance": None,
                "note": "No crystallographic waters",
            }
            logger.info(f"  {pdb_id} ({info['ligand_type']:15s}): NO WATERS")
            continue

        water_oxygens = waters[waters["atom_name"] == "O"][["auth_seq_id", "x", "y", "z"]]
        n655_coords = n655_atoms[["x", "y", "z"]].values

        min_dist = float("inf")
        for _, water in water_oxygens.iterrows():
            water_coord = np.array([water["x"], water["y"], water["z"]])
            for n655_coord in n655_coords:
                dist = np.linalg.norm(water_coord - n655_coord)
                if dist < min_dist:
                    min_dist = dist

        has_contact = min_dist <= WATER_CUTOFF

        results[pdb_id] = {
            "ligand_type": info["ligand_type"],
            "state": info["state"],
            "receptor": info["receptor"],
            "has_water": has_contact,
            "min_distance": min_dist,
        }

        status = "YES" if has_contact else "NO"
        logger.info(f"  {pdb_id} ({info['ligand_type']:15s}, {info['state']:12s}): {status} ({min_dist:.2f} A)")

    return results


# =============================================================================
# Step 6: Visualization
# =============================================================================

def create_hypothesis_visualization(
    h1_results: Dict[str, Dict[str, Any]],
    h2_results: Dict[str, Dict[str, Any]],
    h3_results: Dict[str, Dict[str, Any]],
    output_path: Path,
) -> None:
    """Create visualization of hypothesis testing results."""
    logger.info("=" * 60)
    logger.info("Step 6a: Creating hypothesis visualization")
    logger.info("=" * 60)

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    pdb_ids = list(STRUCTURES.keys())
    ligand_types = [STRUCTURES[p]["ligand_type"] for p in pdb_ids]

    type_colors = {
        "full_agonist": "#2ca02c",
        "partial_agonist": "#ff7f0e",
        "antagonist": "#9467bd",
        "inverse_agonist": "#d62728",
    }
    colors = [type_colors.get(t, "#888") for t in ligand_types]

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "H1: Distance to S5.43 (lower = agonist-like)",
            "H2: Distance to W6.48 (lower = inverse agonist-like)",
            "H3: Water at N6.55 (active state marker)",
            "Summary: Mean Distances by Ligand Type",
        ),
        specs=[
            [{"type": "bar"}, {"type": "bar"}],
            [{"type": "bar"}, {"type": "table"}],
        ],
        vertical_spacing=0.15,
        horizontal_spacing=0.1,
    )

    # H1: S5.43 distance
    h1_distances = [h1_results.get(p, {}).get("distance", 0) for p in pdb_ids]
    fig.add_trace(
        go.Bar(
            x=pdb_ids,
            y=h1_distances,
            marker_color=colors,
            text=[f"{d:.1f}A" if d else "" for d in h1_distances],
            textposition="outside",
        ),
        row=1, col=1,
    )

    # H2: W6.48 distance
    h2_distances = [h2_results.get(p, {}).get("distance", 0) for p in pdb_ids]
    fig.add_trace(
        go.Bar(
            x=pdb_ids,
            y=h2_distances,
            marker_color=colors,
            text=[f"{d:.1f}A" if d else "" for d in h2_distances],
            textposition="outside",
        ),
        row=1, col=2,
    )

    # H3: Water at N6.55
    h3_water = []
    h3_labels = []
    for p in pdb_ids:
        r = h3_results.get(p, {})
        if r.get("has_water"):
            h3_water.append(1)
            h3_labels.append("YES")
        else:
            h3_water.append(0)
            h3_labels.append("NO")

    fig.add_trace(
        go.Bar(
            x=pdb_ids,
            y=h3_water,
            marker_color=colors,
            text=h3_labels,
            textposition="outside",
        ),
        row=2, col=1,
    )

    # Summary table
    summary_data = []
    for lt in ["full_agonist", "partial_agonist", "antagonist", "inverse_agonist"]:
        h1_vals = [r["distance"] for p, r in h1_results.items() if r.get("ligand_type") == lt and r.get("distance")]
        h2_vals = [r["distance"] for p, r in h2_results.items() if r.get("ligand_type") == lt and r.get("distance")]
        h1_mean = f"{np.mean(h1_vals):.2f}A" if h1_vals else "-"
        h2_mean = f"{np.mean(h2_vals):.2f}A" if h2_vals else "-"
        summary_data.append([lt.replace("_", " ").title(), h1_mean, h2_mean])

    fig.add_trace(
        go.Table(
            header=dict(
                values=["Ligand Type", "S5.43 Mean", "W6.48 Mean"],
                fill_color="#1f77b4",
                font=dict(color="white", size=11),
            ),
            cells=dict(
                values=list(zip(*summary_data)) if summary_data else [[], [], []],
                fill_color=[["#f0f0f0", "white"] * 2],
            ),
        ),
        row=2, col=2,
    )

    fig.update_layout(
        title="GPCR Activation Mechanism: Hypothesis Testing",
        height=800,
        width=1200,
        showlegend=False,
        template="plotly_white",
    )

    fig.update_yaxes(title_text="Distance (A)", row=1, col=1)
    fig.update_yaxes(title_text="Distance (A)", row=1, col=2)
    fig.update_yaxes(title_text="Water Present", tickvals=[0, 1], ticktext=["No", "Yes"], row=2, col=1)

    html_path = output_path / "hypothesis_analysis.html"
    fig.write_html(str(html_path))
    logger.info(f"  Saved: {html_path}")

    try:
        png_path = output_path / "hypothesis_analysis.png"
        fig.write_image(str(png_path), scale=2)
        logger.info(f"  Saved: {png_path}")
    except Exception as e:
        logger.warning(f"  PNG save failed: {e}")


def generate_pymol_script(
    grn_mappings: Dict[str, Dict[str, int]],
    output_path: Path,
) -> None:
    """Generate PyMOL visualization script using GRN mappings."""
    logger.info("=" * 60)
    logger.info("Step 6b: Generating PyMOL script")
    logger.info("=" * 60)

    lines = [
        "# GPCR Mechanism Analysis - PyMOL Visualization",
        "# Generated by ProtOS GPCR Binding Pocket Workflow",
        "",
        "# Fetch structures",
    ]

    for pdb_id in STRUCTURES:
        lines.append(f"fetch {pdb_id.lower()}")

    lines.extend([
        "",
        "# Basic setup",
        "bg_color white",
        "set ray_shadows, 0",
        "hide everything",
        "",
        "# Extract GPCR chains and show as cartoon",
    ])

    for pdb_id, info in STRUCTURES.items():
        chain = info["chain"]
        color = info["color"]
        obj_name = f"{pdb_id}_gpcr"
        lines.append(f"create {obj_name}, {pdb_id.lower()} and chain {chain}")
        lines.append(f"show cartoon, {obj_name}")
        lines.append(f"color {color}, {obj_name}")
        lines.append(f"set cartoon_transparency, 0.7, {obj_name}")

    lines.extend(["", "# Show ligands"])

    for pdb_id, info in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        ligand = info["ligand"]
        lines.append(f"show sticks, {obj_name} and resn {ligand}")
        if info["ligand_type"] in ["full_agonist", "partial_agonist"]:
            lines.append(f"color green, {obj_name} and resn {ligand}")
        else:
            lines.append(f"color red, {obj_name} and resn {ligand}")

    lines.extend([
        "",
        "# Key GRN positions for hypothesis testing",
        "# Residue numbers from ProtOS GRN annotation:",
    ])

    # Add GRN residue selections
    for grn_pos in ["5.43", "6.48", "6.55"]:
        grn_safe = grn_pos.replace(".", "_")
        sel_parts = []

        for pdb_id, mapping in grn_mappings.items():
            obj_name = f"{pdb_id}_gpcr"
            chain = STRUCTURES[pdb_id]["chain"]
            res_id = mapping.get(grn_pos)
            if res_id:
                sel_parts.append(f"({obj_name} and chain {chain} and resi {res_id})")
                lines.append(f"# {pdb_id} {grn_pos} = residue {res_id}")

        if sel_parts:
            lines.append(f"select grn_{grn_safe}, {' or '.join(sel_parts)}")
            lines.append(f"show sticks, grn_{grn_safe} and sidechain")

    lines.extend([
        "",
        "# Align all to reference (2RH1)",
        "super 3sn6_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 4ldo_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 2y02_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 2y04_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 2y00_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 3ny9_gpcr and name CA, 2rh1_gpcr and name CA",
        "super 2vt4_gpcr and name CA, 2rh1_gpcr and name CA",
        "",
        "# Create scenes",
        "orient",
        "zoom all, 5",
        "scene overview, store",
        "",
        "# Color GRN positions",
        "color magenta, grn_5_43 and sidechain",
        "color cyan, grn_6_48 and sidechain",
        "color yellow, grn_6_55 and sidechain",
        "",
        "# Create groups",
        "group agonists, 3sn6_gpcr 4ldo_gpcr 2y02_gpcr",
        "group partial_agonists, 2y04_gpcr 2y00_gpcr",
        "group inverse_agonists, 2rh1_gpcr 3ny9_gpcr",
        "group antagonists, 2vt4_gpcr",
        "",
        "# Usage:",
        "# scene overview, recall",
        "# enable agonists; disable inverse_agonists",
        "",
    ])

    pml_path = output_path / "gpcr_mechanism_analysis.pml"
    with open(pml_path, "w") as f:
        f.write("\n".join(lines))

    logger.info(f"  Saved: {pml_path}")


def save_results(
    h1_results: Dict[str, Dict[str, Any]],
    h2_results: Dict[str, Dict[str, Any]],
    h3_results: Dict[str, Dict[str, Any]],
    grn_mappings: Dict[str, Dict[str, int]],
    output_path: Path,
) -> None:
    """Save all analysis results."""
    logger.info("=" * 60)
    logger.info("Step 6c: Saving results")
    logger.info("=" * 60)

    # Hypothesis results as JSON
    results_data = {
        "metadata": {
            "structures": {k: {kk: vv for kk, vv in v.items()} for k, v in STRUCTURES.items()},
            "hypotheses": {
                "H1": "Agonists bind CLOSER to S5.43 than inverse agonists",
                "H2": "Inverse agonists bind CLOSER to W6.48 than agonists",
                "H3": "Water at N6.55 is EXCLUSIVE to agonist-bound active structures",
            },
        },
        "grn_mappings": grn_mappings,
        "H1_S543_distance": h1_results,
        "H2_W648_distance": h2_results,
        "H3_N655_water": h3_results,
    }

    json_path = output_path / "gpcr_mechanism_results.json"
    with open(json_path, "w") as f:
        json.dump(results_data, f, indent=2, default=str)
    logger.info(f"  Saved: {json_path}")

    # Summary CSV
    rows = []
    for pdb_id, info in STRUCTURES.items():
        h1 = h1_results.get(pdb_id, {})
        h2 = h2_results.get(pdb_id, {})
        h3 = h3_results.get(pdb_id, {})
        grn = grn_mappings.get(pdb_id, {})

        rows.append({
            "pdb_id": pdb_id,
            "receptor": info["receptor"],
            "ligand": info["ligand_name"],
            "ligand_type": info["ligand_type"],
            "state": info["state"],
            "grn_5.43_resid": grn.get("5.43"),
            "grn_6.48_resid": grn.get("6.48"),
            "grn_6.55_resid": grn.get("6.55"),
            "H1_S543_dist": round(h1.get("distance", 0), 2) if h1.get("distance") else None,
            "H2_W648_dist": round(h2.get("distance", 0), 2) if h2.get("distance") else None,
            "H3_water": h3.get("has_water", False),
            "H3_water_dist": round(h3.get("min_distance", 0), 2) if h3.get("min_distance") else None,
        })

    summary_df = pd.DataFrame(rows)
    csv_path = output_path / "gpcr_mechanism_summary.csv"
    summary_df.to_csv(csv_path, index=False)
    logger.info(f"  Saved: {csv_path}")

    # Print summary
    logger.info("\n  Summary table:")
    logger.info(f"  {'PDB':5s} {'Type':15s} {'S5.43':>8s} {'W6.48':>8s} {'Water':>6s}")
    logger.info("  " + "-" * 45)
    for row in rows:
        h1_dist = f"{row['H1_S543_dist']:.2f}A" if row['H1_S543_dist'] else "N/A"
        h2_dist = f"{row['H2_W648_dist']:.2f}A" if row['H2_W648_dist'] else "N/A"
        water = "Yes" if row['H3_water'] else "No"
        logger.info(f"  {row['pdb_id']:5s} {row['ligand_type']:15s} {h1_dist:>8s} {h2_dist:>8s} {water:>6s}")


# =============================================================================
# Main Workflow
# =============================================================================

def main() -> int:
    """Run the GPCR binding pocket workflow using ProtOS utilities."""
    logger.info("=" * 70)
    logger.info("GPCR BINDING POCKET WORKFLOW")
    logger.info("Using ProtOS for structure loading, GRN annotation, and ligand analysis")
    logger.info("=" * 70)

    logger.info("\n8 Adrenergic Receptor Structures:")
    for pdb_id, meta in STRUCTURES.items():
        logger.info(
            f"  {pdb_id}: {meta['receptor']} + {meta['ligand_name']} "
            f"({meta['ligand_type']}, {meta['state']})"
        )

    logger.info("\nHypotheses to test:")
    logger.info("  H1: Agonists bind CLOSER to S5.43 than inverse agonists")
    logger.info("  H2: Inverse agonists bind CLOSER to W6.48 (toggle switch)")
    logger.info("  H3: Water at N6.55 is EXCLUSIVE to active structures")

    # Initialize ProtOS
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))

    # Initialize StructureProcessor
    struct_proc = StructureProcessor()

    # Step 1: Load structures
    structures = load_structures(struct_proc)
    if len(structures) < 4:
        logger.error("Not enough structures loaded")
        return 1

    # Step 2: Annotate with GRN
    annotated = annotate_structures_with_grn(struct_proc, structures)

    # Step 3: Ligand interaction analysis
    interaction_results = analyze_ligand_interactions(annotated)
    if not interaction_results:
        logger.error("No interaction results")
        return 1

    # Step 4: Build GRN mappings
    grn_mappings = build_grn_residue_mapping(interaction_results)

    # Step 5: Hypothesis testing
    logger.info("\n" + "=" * 70)
    logger.info("HYPOTHESIS TESTING")
    logger.info("=" * 70)

    h1_results = test_hypothesis_1(interaction_results, grn_mappings)
    h2_results = test_hypothesis_2(interaction_results, grn_mappings)
    h3_results = test_hypothesis_3(interaction_results, grn_mappings)

    # Step 6: Visualization and output
    logger.info("\n" + "=" * 70)
    logger.info("GENERATING OUTPUTS")
    logger.info("=" * 70)

    create_hypothesis_visualization(h1_results, h2_results, h3_results, OUTPUT_DIR)
    generate_pymol_script(grn_mappings, FIGURES_DIR)
    save_results(h1_results, h2_results, h3_results, grn_mappings, OUTPUT_DIR)

    # Final summary
    logger.info("\n" + "=" * 70)
    logger.info("WORKFLOW COMPLETE")
    logger.info("=" * 70)

    # Calculate summary statistics
    agonist_h1 = [r["distance"] for r in h1_results.values() if r.get("ligand_type") == "full_agonist" and r.get("distance")]
    inverse_h1 = [r["distance"] for r in h1_results.values() if r.get("ligand_type") == "inverse_agonist" and r.get("distance")]
    agonist_h2 = [r["distance"] for r in h2_results.values() if r.get("ligand_type") == "full_agonist" and r.get("distance")]
    inverse_h2 = [r["distance"] for r in h2_results.values() if r.get("ligand_type") == "inverse_agonist" and r.get("distance")]

    if agonist_h1 and inverse_h1:
        logger.info(f"\nH1 Summary: Full agonists S5.43 mean={np.mean(agonist_h1):.2f}A vs inverse agonists {np.mean(inverse_h1):.2f}A")
        if np.mean(agonist_h1) < np.mean(inverse_h1):
            logger.info("  -> SUPPORTS H1: Agonists bind closer to S5.43")
        else:
            logger.info("  -> DOES NOT support H1")

    if agonist_h2 and inverse_h2:
        logger.info(f"H2 Summary: Full agonists W6.48 mean={np.mean(agonist_h2):.2f}A vs inverse agonists {np.mean(inverse_h2):.2f}A")
        if np.mean(inverse_h2) < np.mean(agonist_h2):
            logger.info("  -> SUPPORTS H2: Inverse agonists bind closer to W6.48")
        else:
            logger.info("  -> DOES NOT support H2")

    active_water = sum(1 for p, r in h3_results.items() if r.get("has_water") and STRUCTURES[p]["state"] in ["active", "active_like"])
    inactive_water = sum(1 for p, r in h3_results.items() if r.get("has_water") and STRUCTURES[p]["state"] in ["inactive", "intermediate"])
    logger.info(f"H3 Summary: Water at N6.55 - Active states: {active_water}, Inactive states: {inactive_water}")

    logger.info(f"\nOutputs: {OUTPUT_DIR}")
    logger.info(f"PyMOL script: {FIGURES_DIR / 'gpcr_mechanism_analysis.pml'}")
    logger.info(f"Log: {log_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

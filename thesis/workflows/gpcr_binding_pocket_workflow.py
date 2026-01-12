#!/usr/bin/env python3
"""GPCR Binding Pocket Workflow: Comprehensive Mechanism Analysis.

This workflow demonstrates:
- Analyzing 8 adrenergic receptor structures (ADRB1 + ADRB2)
- GRN-annotated ligand-protein interactions
- Agonist vs antagonist vs inverse agonist binding patterns
- Hypothesis testing for activation mechanisms

Structures analyzed:
  ADRB2:
    - 3SN6: full_agonist (BI-167107 + Gs), active
    - 4LDO: full_agonist (adrenaline + Nb6B9), active
    - 2RH1: inverse_agonist (carazolol), inactive
    - 3NY9: inverse_agonist (ICI 118,551), inactive

  ADRB1:
    - 2Y02: full_agonist (isoprenaline), active_like
    - 2Y04: partial_agonist (salbutamol), intermediate
    - 2Y00: partial_agonist (dobutamine), intermediate
    - 2VT4: antagonist (cyanopindolol), inactive

Hypotheses tested:
  H1: Agonists bind CLOSER to S5.43 (serine) than inverse agonists
  H2: Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists
  H3: Water at N6.55 is EXCLUSIVE to agonist-bound active structures

Processors: Structure -> Sequence -> GRN -> Property
Question: "How do ligand types differ in binding mechanisms?"
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Set, Tuple
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
from protos.processing.sequence import SequenceProcessor
from protos.io.ingest.structure_loader import StructureLoader
from protos.analysis.structure_ligand_analysis import (
    extract_all_ligands,
    get_binding_site,
)


# =============================================================================
# Structure Configuration: 8 Adrenergic Receptor Structures
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
        "doi": "10.1038/nature10361",
        "notes": "BI-167107 + Gs protein, fully active (Rasmussen 2011)",
        "color": "#2ca02c",  # Green
    },
    "4LDO": {
        "chain": "A",
        "ligand": "ALE",
        "ligand_name": "Adrenaline",
        "ligand_type": "full_agonist",
        "state": "active",
        "receptor": "ADRB2",
        "doi": "10.1038/nature12572",
        "notes": "Adrenaline + Nb6B9 nanobody, active state (Ring 2013)",
        "color": "#98df8a",  # Light green
    },
    # ADRB1 - Full Agonist (Active-like)
    "2Y02": {
        "chain": "A",
        "ligand": "WHJ",
        "ligand_name": "Isoprenaline",
        "ligand_type": "full_agonist",
        "state": "active_like",
        "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Isoprenaline-bound, agonist-induced changes (Warne 2011)",
        "color": "#17becf",  # Cyan
    },
    # ADRB1 - Partial Agonists (Intermediate)
    "2Y04": {
        "chain": "A",
        "ligand": "68H",
        "ligand_name": "Salbutamol",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Salbutamol-bound, intermediate state (Warne 2011)",
        "color": "#ff7f0e",  # Orange
    },
    "2Y00": {
        "chain": "A",
        "ligand": "Y00",
        "ligand_name": "Dobutamine",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Dobutamine-bound, intermediate state (Warne 2011)",
        "color": "#ffbb78",  # Light orange
    },
    # ADRB2 - Inverse Agonists (Inactive)
    "2RH1": {
        "chain": "A",
        "ligand": "CAU",
        "ligand_name": "Carazolol",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
        "doi": "10.1126/science.1150577",
        "notes": "Carazolol-bound, first high-res GPCR (Cherezov 2007)",
        "color": "#d62728",  # Red
    },
    "3NY9": {
        "chain": "A",
        "ligand": "JSZ",
        "ligand_name": "ICI 118,551",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
        "doi": "10.1021/ja105108q",
        "notes": "ICI 118,551-bound inverse agonist (Wacker 2010)",
        "color": "#ff9896",  # Light red
    },
    # ADRB1 - Antagonist (Inactive)
    "2VT4": {
        "chain": "A",
        "ligand": "P32",
        "ligand_name": "Cyanopindolol",
        "ligand_type": "antagonist",
        "state": "inactive",
        "receptor": "ADRB1",
        "doi": "10.1038/nature07101",
        "notes": "Cyanopindolol-bound, inactive (Warne 2008)",
        "color": "#9467bd",  # Purple
    },
}

# Key GRN positions for binding and hypothesis testing
KEY_GRN_POSITIONS = {
    "3.28": "Asp - salt bridge with amine",
    "3.32": "Asp - key anionic anchor",
    "3.33": "Val",
    "3.36": "Thr/Cys",
    "5.42": "Ser - agonist interaction",
    "5.43": "Ser - H1: agonist binding marker",
    "5.46": "Ser - agonist interaction",
    "6.48": "Trp - H2: toggle switch",
    "6.51": "Phe",
    "6.52": "Phe",
    "6.55": "Asn - H3: water marker",
    "7.35": "Tyr",
    "7.39": "Asn",
}

# Hypothesis-specific GRN positions
H1_GRN = "5.43"  # Agonists bind closer to S5.43
H2_GRN = "6.48"  # Inverse agonists bind closer to W6.48 (toggle switch)
H3_GRN = "6.55"  # Water at N6.55 exclusive to active structures

BINDING_CUTOFF = 5.0  # Angstroms
WATER_CUTOFF = 3.5    # Angstroms for water contacts
DATASET_NAME = "gpcr_mechanism_analysis"
GRN_TABLE_NAME = "gpcr_mechanism_grn"
REFERENCE_TABLE = "gpcrdb_ref"
PROTEIN_FAMILY = "gpcr_a"


# =============================================================================
# Part 1: Structure Loading and Basic Analysis
# =============================================================================

def download_structures() -> Dict[str, pd.DataFrame]:
    """Download all 8 adrenergic receptor structures from PDB."""
    logger.info("=" * 60)
    logger.info("Step 1: Downloading 8 GPCR structures")
    logger.info("=" * 60)

    struct_proc = StructureProcessor()
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

            # Download from PDB
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


def extract_ligands_and_binding_sites(
    structures: Dict[str, pd.DataFrame],
) -> Dict[str, Dict[str, Any]]:
    """Extract ligands and their binding site residues for all structures."""
    logger.info("=" * 60)
    logger.info("Step 2: Extracting ligands and binding sites")
    logger.info("=" * 60)

    struct_proc = StructureProcessor()
    binding_data = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        target_ligand = info["ligand"]

        logger.info(f"  {pdb_id} - {info['ligand_name']} ({info['ligand_type']})...")

        try:
            struct_proc.frames[pdb_id.lower()] = df

            # Find the target ligand
            ligands = extract_all_ligands(
                struct_proc,
                pdb_id.lower(),
                exclude_common=True,
                allowed_res_names={target_ligand},
            )

            if not ligands:
                logger.warning(f"    Ligand {target_ligand} not found")
                continue

            ligand_info = ligands[0]

            # Get binding site
            binding_site = get_binding_site(
                struct_proc,
                pdb_id.lower(),
                ligand_info["atoms"],
                cutoff=BINDING_CUTOFF,
                include_backbone=True,
            )

            binding_residues = binding_site.get("residues", pd.DataFrame())
            n_residues = len(binding_residues) if isinstance(binding_residues, pd.DataFrame) else 0

            logger.info(f"    Ligand: {ligand_info['num_atoms']} atoms, Binding site: {n_residues} residues")

            binding_data[pdb_id] = {
                "ligand": ligand_info,
                "binding_site": binding_site,
                "structure_info": info,
                "n_binding_residues": n_residues,
            }

        except Exception as e:
            logger.warning(f"    Failed: {e}")

    return binding_data


def extract_and_annotate_sequences(
    structures: Dict[str, pd.DataFrame],
    binding_data: Dict[str, Dict[str, Any]],
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Extract sequences and annotate with GRN."""
    logger.info("=" * 60)
    logger.info("Step 3: Extracting sequences and GRN annotation")
    logger.info("=" * 60)

    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()

    sequences = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        chain_id = info["chain"]

        try:
            struct_proc.frames[pdb_id.lower()] = df
            sequence = struct_proc.get_sequence(pdb_id.lower(), chain_id)

            if sequence:
                seq_id = f"{pdb_id}_{chain_id}"
                sequences[seq_id] = sequence

                seq_proc.save_entity(seq_id, sequence, metadata={
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "receptor": info["receptor"],
                    "state": info["state"],
                    "ligand_type": info["ligand_type"],
                })

        except Exception as e:
            logger.warning(f"  {pdb_id}: Failed - {e}")

    # Create dataset and annotate with GRN
    if sequences:
        seq_proc.save_sequences(
            sequences,
            output_file=DATASET_NAME,
            dataset_name=DATASET_NAME,
            metadata={"description": "GPCR mechanism analysis dataset"},
            materialize_entities=True,
        )

        grn_table, summary = seq_proc.annotate_with_grn(
            dataset_name=DATASET_NAME,
            reference_table=REFERENCE_TABLE,
            protein_family=PROTEIN_FAMILY,
            output_table=GRN_TABLE_NAME,
            return_summary=True,
            allow_create=True,
        )

        logger.info(f"  Annotated {summary['global']['annotated']}/{len(sequences)} sequences with GRN")

        # Map GRN to binding residues
        for pdb_id, data in binding_data.items():
            info = STRUCTURES[pdb_id]
            seq_id = f"{pdb_id}_{info['chain']}"

            if seq_id not in grn_table.index:
                continue

            grn_row = grn_table.loc[seq_id]

            # Create residue number to GRN mapping
            res_to_grn = {}
            for grn_pos, residue_info in grn_row.items():
                if residue_info and residue_info != "-" and isinstance(residue_info, str):
                    try:
                        res_num = int("".join(c for c in residue_info if c.isdigit()))
                        res_to_grn[res_num] = grn_pos
                    except ValueError:
                        pass

            data["res_to_grn"] = res_to_grn

            # Annotate binding residues
            binding_site = data.get("binding_site", {})
            binding_residues = binding_site.get("residues", pd.DataFrame())

            if isinstance(binding_residues, pd.DataFrame) and not binding_residues.empty:
                binding_residues = binding_residues.copy()
                res_col = "res_id" if "res_id" in binding_residues.columns else "auth_seq_id"
                if res_col in binding_residues.columns:
                    binding_residues["grn_position"] = binding_residues[res_col].map(res_to_grn)
                data["binding_residues_grn"] = binding_residues

        return sequences, grn_table

    return sequences, pd.DataFrame()


# =============================================================================
# Part 2: Hypothesis Testing
# =============================================================================

def calculate_ligand_to_grn_distance(
    df: pd.DataFrame,
    chain: str,
    ligand_code: str,
    grn_pos: str,
    res_to_grn: Dict[int, str],
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    """Calculate minimum distance from ligand to a GRN-annotated residue."""
    # Get ligand atoms (non-hydrogen)
    ligand_atoms = df[(df["res_name3l"] == ligand_code) & (~df["element"].isin(["H"]))]
    if ligand_atoms.empty:
        return None, None, None

    # Find residue with target GRN
    target_res_id = None
    for res_id, grn in res_to_grn.items():
        if grn == grn_pos:
            target_res_id = res_id
            break

    if target_res_id is None:
        return None, None, None

    # Get residue atoms
    res_atoms = df[
        (df["auth_chain_id"] == chain) &
        (df["auth_seq_id"] == target_res_id) &
        (~df["element"].isin(["H"]))
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


def analyze_hypothesis_1(
    structures: Dict[str, pd.DataFrame],
    binding_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """H1: Agonists bind CLOSER to S5.43 than inverse agonists."""
    logger.info("=" * 60)
    logger.info("H1: Ligand distance to S5.43 (agonist binding marker)")
    logger.info("=" * 60)

    results = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]
        res_to_grn = binding_data.get(pdb_id, {}).get("res_to_grn", {})

        df_reset = df.reset_index() if df.index.name else df

        dist, lig_atom, res_atom = calculate_ligand_to_grn_distance(
            df_reset, chain, ligand, H1_GRN, res_to_grn
        )

        if dist is not None:
            results[pdb_id] = {
                "ligand_type": info["ligand_type"],
                "state": info["state"],
                "receptor": info["receptor"],
                "distance": dist,
                "ligand_atom": lig_atom,
                "residue_atom": res_atom,
            }
            logger.info(f"  {pdb_id} ({info['ligand_type']:15s}): {dist:.2f} A")

    return results


def analyze_hypothesis_2(
    structures: Dict[str, pd.DataFrame],
    binding_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """H2: Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists."""
    logger.info("=" * 60)
    logger.info("H2: Ligand distance to W6.48 (toggle switch)")
    logger.info("=" * 60)

    results = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]
        res_to_grn = binding_data.get(pdb_id, {}).get("res_to_grn", {})

        df_reset = df.reset_index() if df.index.name else df

        dist, lig_atom, res_atom = calculate_ligand_to_grn_distance(
            df_reset, chain, ligand, H2_GRN, res_to_grn
        )

        if dist is not None:
            results[pdb_id] = {
                "ligand_type": info["ligand_type"],
                "state": info["state"],
                "receptor": info["receptor"],
                "distance": dist,
                "ligand_atom": lig_atom,
                "residue_atom": res_atom,
            }
            logger.info(f"  {pdb_id} ({info['ligand_type']:15s}): {dist:.2f} A")

    return results


def analyze_hypothesis_3(
    structures: Dict[str, pd.DataFrame],
    binding_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """H3: Water at N6.55 is EXCLUSIVE to agonist-bound active structures."""
    logger.info("=" * 60)
    logger.info("H3: Water at N6.55 (activation marker)")
    logger.info("=" * 60)
    logger.info("  NOTE: Crystallographic waters depend on resolution")

    results = {}

    for pdb_id, df in structures.items():
        info = STRUCTURES[pdb_id]
        chain = info["chain"]
        res_to_grn = binding_data.get(pdb_id, {}).get("res_to_grn", {})

        df_reset = df.reset_index() if df.index.name else df

        # Find residue with GRN 6.55
        target_res_id = None
        for res_id, grn in res_to_grn.items():
            if grn == H3_GRN:
                target_res_id = res_id
                break

        if target_res_id is None:
            continue

        # Get N6.55 atoms
        n655_atoms = df_reset[
            (df_reset["auth_chain_id"] == chain) &
            (df_reset["auth_seq_id"] == target_res_id) &
            (~df_reset["element"].isin(["H"]))
        ]

        if n655_atoms.empty:
            continue

        # Get water oxygens
        waters = df_reset[df_reset["res_name3l"] == "HOH"]

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
# Part 3: Analysis and Visualization
# =============================================================================

def summarize_by_ligand_type(
    results: Dict[str, Dict[str, Any]],
    metric_key: str = "distance",
) -> Dict[str, Dict[str, Any]]:
    """Summarize results by ligand type."""
    by_type = defaultdict(list)

    for pdb_id, data in results.items():
        lig_type = data["ligand_type"]
        value = data.get(metric_key)
        if value is not None:
            by_type[lig_type].append((pdb_id, value))

    summary = {}
    for lig_type in ["full_agonist", "partial_agonist", "antagonist", "inverse_agonist"]:
        if lig_type in by_type:
            values = [v for _, v in by_type[lig_type]]
            pdbs = [p for p, _ in by_type[lig_type]]
            summary[lig_type] = {
                "mean": np.mean(values),
                "min": min(values),
                "max": max(values),
                "structures": pdbs,
            }

    return summary


def create_hypothesis_visualization(
    h1_results: Dict[str, Dict[str, Any]],
    h2_results: Dict[str, Dict[str, Any]],
    h3_results: Dict[str, Dict[str, Any]],
    output_path: Path,
) -> None:
    """Create visualization of hypothesis testing results."""
    logger.info("=" * 60)
    logger.info("Creating hypothesis visualization")
    logger.info("=" * 60)

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Prepare data
    pdb_ids = list(STRUCTURES.keys())
    ligand_types = [STRUCTURES[p]["ligand_type"] for p in pdb_ids]

    # Colors by ligand type
    type_colors = {
        "full_agonist": "#2ca02c",
        "partial_agonist": "#ff7f0e",
        "antagonist": "#9467bd",
        "inverse_agonist": "#d62728",
    }
    colors = [type_colors.get(t, "#888") for t in ligand_types]

    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "H1: Distance to S5.43 (lower = agonist-like)",
            "H2: Distance to W6.48 (lower = inverse agonist-like)",
            "H3: Water at N6.55 (active state marker)",
            "Summary: Ligand Type Discrimination",
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
            hovertemplate="%{x}<br>S5.43 dist: %{y:.2f}A<extra></extra>",
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
            hovertemplate="%{x}<br>W6.48 dist: %{y:.2f}A<extra></extra>",
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
            hovertemplate="%{x}<br>Water: %{text}<extra></extra>",
        ),
        row=2, col=1,
    )

    # Summary table
    h1_summary = summarize_by_ligand_type(h1_results)
    h2_summary = summarize_by_ligand_type(h2_results)

    summary_data = []
    for lt in ["full_agonist", "partial_agonist", "antagonist", "inverse_agonist"]:
        h1_mean = h1_summary.get(lt, {}).get("mean", "-")
        h2_mean = h2_summary.get(lt, {}).get("mean", "-")
        summary_data.append([
            lt.replace("_", " ").title(),
            f"{h1_mean:.2f}A" if isinstance(h1_mean, float) else "-",
            f"{h2_mean:.2f}A" if isinstance(h2_mean, float) else "-",
        ])

    fig.add_trace(
        go.Table(
            header=dict(
                values=["Ligand Type", "S5.43 Mean", "W6.48 Mean"],
                fill_color="#1f77b4",
                font=dict(color="white", size=11),
                align="center",
            ),
            cells=dict(
                values=list(zip(*summary_data)) if summary_data else [[], [], []],
                fill_color=[["#f0f0f0", "white"] * 2],
                align="center",
            ),
        ),
        row=2, col=2,
    )

    fig.update_layout(
        title=dict(
            text="GPCR Activation Mechanism: Hypothesis Testing (H1-H3)",
            font=dict(size=18),
        ),
        height=800,
        width=1200,
        showlegend=False,
        template="plotly_white",
    )

    fig.update_yaxes(title_text="Distance (A)", row=1, col=1)
    fig.update_yaxes(title_text="Distance (A)", row=1, col=2)
    fig.update_yaxes(title_text="Water Present", tickvals=[0, 1], ticktext=["No", "Yes"], row=2, col=1)

    # Save
    html_path = output_path / "hypothesis_analysis.html"
    fig.write_html(str(html_path))
    logger.info(f"  Saved: {html_path}")

    # Also save as PNG
    try:
        png_path = output_path / "hypothesis_analysis.png"
        fig.write_image(str(png_path), scale=2)
        logger.info(f"  Saved: {png_path}")
    except Exception as e:
        logger.warning(f"  Could not save PNG: {e}")


def generate_pymol_script(
    binding_data: Dict[str, Dict[str, Any]],
    h1_results: Dict[str, Dict[str, Any]],
    h2_results: Dict[str, Dict[str, Any]],
    h3_results: Dict[str, Dict[str, Any]],
    output_path: Path,
) -> None:
    """Generate PyMOL script for visualization."""
    logger.info("=" * 60)
    logger.info("Generating PyMOL visualization script")
    logger.info("=" * 60)

    lines = [
        "# GPCR Mechanism Analysis - PyMOL Visualization",
        "# Generated by ProtOS GPCR Binding Pocket Workflow",
        "#",
        "# Shows 8 adrenergic receptor structures with hypothesis markers",
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
        "set antialias, 2",
        "hide everything",
        "",
        "# Show as cartoon",
    ])

    for pdb_id, info in STRUCTURES.items():
        chain = info["chain"]
        color = info["color"]
        lines.append(f"show cartoon, {pdb_id.lower()} and chain {chain}")
        lines.append(f"color {color}, {pdb_id.lower()}")

    lines.extend([
        "",
        "# Show ligands",
    ])

    for pdb_id, info in STRUCTURES.items():
        ligand = info["ligand"]
        lines.append(f"show sticks, {pdb_id.lower()} and resn {ligand}")
        lines.append(f"color yellow, {pdb_id.lower()} and resn {ligand}")

    lines.extend([
        "",
        "# Highlight key GRN positions",
        "# H1: S5.43 (agonist binding marker)",
        "# H2: W6.48 (toggle switch)",
        "# H3: N6.55 (water marker)",
        "",
    ])

    # Add selections for key residues
    for pdb_id, data in binding_data.items():
        info = STRUCTURES[pdb_id]
        chain = info["chain"]
        res_to_grn = data.get("res_to_grn", {})

        for grn, res_id in [(k, v) for v, k in res_to_grn.items()]:
            if grn in ["5.43", "6.48", "6.55"]:
                lines.append(f"# {pdb_id} {grn}")
                lines.append(f"show sticks, {pdb_id.lower()} and chain {chain} and resi {res_id}")

    lines.extend([
        "",
        "# Create scenes for each hypothesis",
        "",
        "# Scene 1: Overview",
        "orient",
        "zoom all, 5",
        "scene overview, store",
        "",
        "# Scene 2: H1 - S5.43 binding",
        "# Focus on agonist binding site",
        "scene H1_S543, store",
        "",
        "# Scene 3: H2 - W6.48 toggle switch",
        "scene H2_W648, store",
        "",
        "# Scene 4: H3 - N6.55 water",
        "scene H3_N655, store",
        "",
        "# Usage:",
        "# scene overview  - show all structures",
        "# scene H1_S543   - focus on S5.43",
        "# scene H2_W648   - focus on W6.48",
        "# scene H3_N655   - focus on N6.55",
        "",
    ])

    pml_path = output_path / "gpcr_mechanism_analysis.pml"
    with open(pml_path, "w") as f:
        f.write("\n".join(lines))

    logger.info(f"  Saved: {pml_path}")
    logger.info(f"  Run in PyMOL: pymol {pml_path}")


def save_results(
    binding_data: Dict[str, Dict[str, Any]],
    h1_results: Dict[str, Dict[str, Any]],
    h2_results: Dict[str, Dict[str, Any]],
    h3_results: Dict[str, Dict[str, Any]],
    output_path: Path,
) -> None:
    """Save all analysis results."""
    logger.info("=" * 60)
    logger.info("Saving results")
    logger.info("=" * 60)

    # Hypothesis results as JSON
    hypothesis_data = {
        "metadata": {
            "structures": {k: {kk: vv for kk, vv in v.items() if kk != "color"} for k, v in STRUCTURES.items()},
            "hypotheses": {
                "H1": "Agonists bind CLOSER to S5.43 than inverse agonists",
                "H2": "Inverse agonists bind CLOSER to W6.48 than agonists",
                "H3": "Water at N6.55 is EXCLUSIVE to agonist-bound active structures",
            },
        },
        "H1_S543_distance": h1_results,
        "H2_W648_distance": h2_results,
        "H3_N655_water": h3_results,
        "summary": {
            "H1": summarize_by_ligand_type(h1_results),
            "H2": summarize_by_ligand_type(h2_results),
        },
    }

    json_path = output_path / "gpcr_mechanism_results.json"
    with open(json_path, "w") as f:
        json.dump(hypothesis_data, f, indent=2, default=str)
    logger.info(f"  Saved: {json_path}")

    # Summary table as CSV
    rows = []
    for pdb_id, info in STRUCTURES.items():
        h1 = h1_results.get(pdb_id, {})
        h2 = h2_results.get(pdb_id, {})
        h3 = h3_results.get(pdb_id, {})

        rows.append({
            "pdb_id": pdb_id,
            "receptor": info["receptor"],
            "ligand": info["ligand_name"],
            "ligand_type": info["ligand_type"],
            "state": info["state"],
            "H1_S543_dist": round(h1.get("distance", 0), 2) if h1.get("distance") else None,
            "H2_W648_dist": round(h2.get("distance", 0), 2) if h2.get("distance") else None,
            "H3_water": h3.get("has_water", False),
            "H3_water_dist": round(h3.get("min_distance", 0), 2) if h3.get("min_distance") else None,
        })

    summary_df = pd.DataFrame(rows)
    csv_path = output_path / "gpcr_mechanism_summary.csv"
    summary_df.to_csv(csv_path, index=False)
    logger.info(f"  Saved: {csv_path}")


def main() -> int:
    """Run the comprehensive GPCR binding pocket workflow."""
    logger.info("=" * 70)
    logger.info("GPCR BINDING POCKET WORKFLOW: COMPREHENSIVE MECHANISM ANALYSIS")
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
    logger.info("")

    # Initialize protos
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))

    # Part 1: Structure Loading and Basic Analysis
    structures = download_structures()
    if len(structures) < 4:
        logger.error("Not enough structures downloaded, exiting")
        return 1

    binding_data = extract_ligands_and_binding_sites(structures)
    if not binding_data:
        logger.error("No binding data extracted, exiting")
        return 1

    sequences, grn_table = extract_and_annotate_sequences(structures, binding_data)

    # Part 2: Hypothesis Testing
    logger.info("\n" + "=" * 70)
    logger.info("HYPOTHESIS TESTING")
    logger.info("=" * 70)

    h1_results = analyze_hypothesis_1(structures, binding_data)
    h2_results = analyze_hypothesis_2(structures, binding_data)
    h3_results = analyze_hypothesis_3(structures, binding_data)

    # Part 3: Visualization and Output
    logger.info("\n" + "=" * 70)
    logger.info("GENERATING OUTPUTS")
    logger.info("=" * 70)

    create_hypothesis_visualization(h1_results, h2_results, h3_results, OUTPUT_DIR)
    generate_pymol_script(binding_data, h1_results, h2_results, h3_results, FIGURES_DIR)
    save_results(binding_data, h1_results, h2_results, h3_results, OUTPUT_DIR)

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("WORKFLOW COMPLETE")
    logger.info("=" * 70)

    logger.info("\nKey Findings:")

    h1_summary = summarize_by_ligand_type(h1_results)
    if "full_agonist" in h1_summary and "inverse_agonist" in h1_summary:
        fa_dist = h1_summary["full_agonist"]["mean"]
        ia_dist = h1_summary["inverse_agonist"]["mean"]
        logger.info(f"  H1: Full agonists S5.43 mean={fa_dist:.2f}A vs inverse agonists {ia_dist:.2f}A")

    h2_summary = summarize_by_ligand_type(h2_results)
    if "full_agonist" in h2_summary and "inverse_agonist" in h2_summary:
        fa_dist = h2_summary["full_agonist"]["mean"]
        ia_dist = h2_summary["inverse_agonist"]["mean"]
        logger.info(f"  H2: Full agonists W6.48 mean={fa_dist:.2f}A vs inverse agonists {ia_dist:.2f}A")

    active_water = sum(1 for p, r in h3_results.items() if r.get("has_water") and STRUCTURES[p]["state"] in ["active", "active_like"])
    inactive_water = sum(1 for p, r in h3_results.items() if r.get("has_water") and STRUCTURES[p]["state"] in ["inactive", "intermediate"])
    logger.info(f"  H3: Water at N6.55 - Active: {active_water}, Inactive: {inactive_water}")

    logger.info(f"\nOutputs saved to: {OUTPUT_DIR}")
    logger.info(f"PyMOL script: {FIGURES_DIR / 'gpcr_mechanism_analysis.pml'}")
    logger.info(f"Log saved to: {log_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

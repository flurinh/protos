#!/usr/bin/env python3
"""StructureProcessor Example: Type I vs Type II opsin alignment on retinal.

This example demonstrates:
- Loading structures from PDB (Type I microbial and Type II animal opsins)
- Extracting retinal ligand contacts
- Aligning structures based on retinal position
- Exporting aligned structures to CIF format
- Visualizing the alignment

KEY INSIGHT: Type I (microbial) and Type II (animal) opsins both bind retinal
but have DIFFERENT 7TM topologies. This example shows how they compare
structurally around the retinal binding site.

Question: "How do Type I and Type II opsins compare structurally around retinal?"
"""

from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
SRC_DIR = REPO_ROOT / "src"
LOG_DIR = THESIS_DIR / "logs"
OUTPUT_DIR = THESIS_DIR / "outputs" / "structure"

if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

LOG_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Configure logging
log_file = LOG_DIR / f"structure_processor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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
from protos.io.ingest.structure_loader import StructureLoader
from protos.io.formats.cif_utils import write_cif_file

# =============================================================================
# Load Color Scheme
# =============================================================================
with open(THESIS_DIR / "colorscales.yaml") as f:
    COLORS = yaml.safe_load(f)

# =============================================================================
# Structure Dataset
# =============================================================================
# Type I (Microbial) and Type II (Animal) opsins with retinal
# KEY: These have DIFFERENT 7TM topologies but both bind retinal
STRUCTURE_DATASET = {
    # Type I: Microbial opsins (non-GPCR 7TM)
    "1C3W": {
        "name": "Bacteriorhodopsin",
        "type": "Type I",
        "type_full": "Type I (Microbial)",
        "chain": "A",
        "ligand": "RET",
        "organism": "Halobacterium salinarum",
        "function": "proton pump",
        "color_hex": COLORS["structures"]["1C3W"]["hex"],
    },
    "3UG9": {
        "name": "Channelrhodopsin-2 (C1C2)",
        "type": "Type I",
        "type_full": "Type I (Microbial)",
        "chain": "A",
        "ligand": "RET",
        "organism": "Chlamydomonas reinhardtii",
        "function": "cation channel",
        "color_hex": COLORS["structures"]["3UG9"]["hex"],
    },
    # Type II: Animal opsins (GPCR family)
    "1U19": {
        "name": "Bovine Rhodopsin",
        "type": "Type II",
        "type_full": "Type II (Animal/GPCR)",
        "chain": "A",
        "ligand": "RET",
        "organism": "Bos taurus",
        "function": "dim light vision",
        "color_hex": COLORS["structures"]["1U19"]["hex"],
    },
    "2Z73": {
        "name": "Squid Rhodopsin",
        "type": "Type II",
        "type_full": "Type II (Animal/GPCR)",
        "chain": "A",
        "ligand": "RET",
        "organism": "Todarodes pacificus",
        "function": "vision",
        "color_hex": COLORS["structures"]["2Z73"]["hex"],
    },
}

# Reference structure for alignment
REFERENCE_STRUCTURE = "1C3W"
CONTACT_CUTOFF = 4.0  # Angstroms
DATASET_NAME = "opsin_retinal_structures"


def download_structures() -> Dict[str, pd.DataFrame]:
    """Download opsin structures from PDB."""
    logger.info("=" * 60)
    logger.info("Step 1: Downloading opsin structures from PDB")
    logger.info("=" * 60)

    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    structures = {}
    pdb_ids = list(STRUCTURE_DATASET.keys())

    for pdb_id in pdb_ids:
        info = STRUCTURE_DATASET[pdb_id]
        logger.info(f"  Loading {pdb_id}: {info['name']} ({info['type']})...")

        try:
            # Try to load from cache
            df = struct_proc.load_entity(pdb_id)

            if df is None or df.empty:
                # Download from PDB
                success, failed = loader.download_batch(
                    [pdb_id],
                    dataset_name=DATASET_NAME,
                    create_dataset=True,
                    overwrite=False,
                )
                if failed:
                    logger.warning(f"    Failed to download {pdb_id}")
                    continue
                df = struct_proc.load_entity(pdb_id)

            if df is not None and not df.empty:
                structures[pdb_id] = df
                n_atoms = len(df)
                n_residues = df.groupby(["auth_chain_id", "auth_seq_id"]).ngroups
                logger.info(f"    Loaded: {n_atoms} atoms, {n_residues} residues")

        except Exception as e:
            logger.warning(f"    Error loading {pdb_id}: {e}")

    logger.info(f"\nLoaded {len(structures)} structures")
    return structures


def analyze_retinal_contacts(
    struct_proc: StructureProcessor,
    structures: Dict[str, pd.DataFrame],
) -> Dict[str, pd.DataFrame]:
    """Analyze retinal binding pocket for each structure."""
    logger.info("=" * 60)
    logger.info("Step 2: Analyzing retinal binding pockets")
    logger.info("=" * 60)

    contact_data = {}

    for pdb_id, df in structures.items():
        info = STRUCTURE_DATASET[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]

        logger.info(f"\n  {pdb_id} ({info['name']}):")

        # Get ligand interactions using StructureProcessor
        try:
            interactions = struct_proc.get_ligand_interactions(
                pdb_id,
                ligand_id=ligand,
                chain_id=chain,
                cutoff=CONTACT_CUTOFF,
            )

            if interactions:
                binding_residues = interactions.get("binding_residues", [])
                logger.info(f"    Binding residues within {CONTACT_CUTOFF}A: {len(binding_residues)}")

                # Create DataFrame of contacts
                if binding_residues:
                    contacts_df = pd.DataFrame(binding_residues)
                    contact_data[pdb_id] = contacts_df

                    # Show top contacts
                    logger.info("    Top contacts by distance:")
                    for i, res in enumerate(binding_residues[:5]):
                        res_name = res.get("res_name", "UNK")
                        res_id = res.get("res_id", "?")
                        min_dist = res.get("min_distance", 0)
                        logger.info(f"      {res_name}{res_id}: {min_dist:.2f}A")

            else:
                logger.warning(f"    No interactions found for {ligand}")

        except Exception as e:
            logger.warning(f"    Error analyzing contacts: {e}")

    return contact_data


def align_on_retinal(
    structures: Dict[str, pd.DataFrame],
    reference_id: str = REFERENCE_STRUCTURE,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, float]]:
    """Align all structures based on retinal and CA atoms."""
    logger.info("=" * 60)
    logger.info(f"Step 3: Aligning structures on retinal (ref: {reference_id})")
    logger.info("=" * 60)

    from protos.analysis.structure.alignment import kabsch_alignment

    if reference_id not in structures:
        logger.error(f"Reference structure {reference_id} not available")
        return {}, {}

    ref_df = structures[reference_id]
    ref_info = STRUCTURE_DATASET[reference_id]
    ref_chain = ref_info["chain"]
    ref_ligand = ref_info["ligand"]

    # Get reference retinal coordinates (only from specified chain)
    ref_ret = ref_df[
        ((ref_df["res_name3l"] == ref_ligand) | (ref_df["res_name"] == ref_ligand)) &
        (ref_df["auth_chain_id"] == ref_chain)
    ][["x", "y", "z"]].values

    if len(ref_ret) == 0:
        logger.error(f"No retinal found in reference structure")
        return {}, {}

    ref_ret_center = ref_ret.mean(axis=0)
    logger.info(f"  Reference retinal center: [{ref_ret_center[0]:.2f}, {ref_ret_center[1]:.2f}, {ref_ret_center[2]:.2f}]")

    # Get reference CA atoms for alignment
    ref_ca = ref_df[
        (ref_df["atom_name"] == "CA") &
        (ref_df["auth_chain_id"] == ref_chain)
    ][["x", "y", "z"]].values

    logger.info(f"  Reference CA atoms: {len(ref_ca)}")

    aligned_structures = {reference_id: ref_df.copy()}
    rmsd_values = {reference_id: 0.0}

    # Align each structure to reference
    for pdb_id, df in structures.items():
        if pdb_id == reference_id:
            continue

        info = STRUCTURE_DATASET[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]

        logger.info(f"\n  Aligning {pdb_id} ({info['name']})...")

        try:
            # Get target CA atoms
            target_ca = df[
                (df["atom_name"] == "CA") &
                (df["auth_chain_id"] == chain)
            ][["x", "y", "z"]].values

            if len(target_ca) == 0:
                logger.warning(f"    No CA atoms found")
                continue

            # Use the shorter of the two for alignment
            n_common = min(len(ref_ca), len(target_ca))
            logger.info(f"    Using {n_common} CA atoms for alignment")

            # Perform Kabsch alignment
            rotation, translation, rmsd = kabsch_alignment(
                ref_ca[:n_common],
                target_ca[:n_common],
            )

            logger.info(f"    RMSD: {rmsd:.2f}A")

            # Apply transformation to all atoms
            aligned_df = df.copy()
            coords = aligned_df[["x", "y", "z"]].values

            # Apply rotation and translation
            aligned_coords = (coords @ rotation.T) + translation

            aligned_df["x"] = aligned_coords[:, 0]
            aligned_df["y"] = aligned_coords[:, 1]
            aligned_df["z"] = aligned_coords[:, 2]

            aligned_structures[pdb_id] = aligned_df
            rmsd_values[pdb_id] = rmsd

            # Show aligned retinal center (only from specified chain)
            aligned_ret = aligned_df[
                ((aligned_df["res_name3l"] == ligand) | (aligned_df["res_name"] == ligand)) &
                (aligned_df["auth_chain_id"] == chain)
            ][["x", "y", "z"]].values

            if len(aligned_ret) > 0:
                aligned_ret_center = aligned_ret.mean(axis=0)
                logger.info(
                    f"    Aligned retinal center: "
                    f"[{aligned_ret_center[0]:.2f}, {aligned_ret_center[1]:.2f}, {aligned_ret_center[2]:.2f}]"
                )

        except Exception as e:
            logger.warning(f"    Alignment failed: {e}")
            import traceback
            traceback.print_exc()

    return aligned_structures, rmsd_values


def export_aligned_structures(
    aligned_structures: Dict[str, pd.DataFrame],
    output_dir: Path,
) -> List[Path]:
    """Export aligned structures to CIF format."""
    logger.info("=" * 60)
    logger.info("Step 4: Exporting aligned structures to CIF")
    logger.info("=" * 60)

    exported_files = []
    cif_dir = output_dir / "aligned_cif"
    cif_dir.mkdir(parents=True, exist_ok=True)

    for pdb_id, df in aligned_structures.items():
        info = STRUCTURE_DATASET[pdb_id]
        output_file = cif_dir / f"{pdb_id}_aligned.cif"

        try:
            # Add structure_id if not present
            if "structure_id" not in df.columns:
                df = df.copy()
                df["structure_id"] = pdb_id

            # Write CIF file
            write_cif_file(str(output_file), df, force_overwrite=True)
            exported_files.append(output_file)

            logger.info(f"  Exported: {output_file.name}")
            logger.info(f"    Structure: {info['name']}")
            logger.info(f"    Atoms: {len(df)}")

        except Exception as e:
            logger.warning(f"  Failed to export {pdb_id}: {e}")

    logger.info(f"\nExported {len(exported_files)} CIF files to: {cif_dir}")
    return exported_files


def visualize_alignment(
    aligned_structures: Dict[str, pd.DataFrame],
    rmsd_values: Dict[str, float],
    output_dir: Path,
) -> None:
    """Create interactive 3D visualization of aligned structures."""
    logger.info("=" * 60)
    logger.info("Step 5: Creating alignment visualization")
    logger.info("=" * 60)

    import plotly.graph_objects as go

    fig = go.Figure()

    # Color scheme from colorscales.yaml
    type_colors = {
        "Type I": COLORS["types"]["type_i"]["hex"],
        "Type II": COLORS["types"]["type_ii"]["hex"],
    }
    retinal_color = COLORS["ligands"]["retinal"]["hex"]

    # Add each structure
    for pdb_id, df in aligned_structures.items():
        info = STRUCTURE_DATASET[pdb_id]
        chain = info["chain"]
        # Use structure-specific color from colorscales
        color = info.get("color_hex", type_colors.get(info["type"], "#888888"))

        # Get CA atoms for backbone trace
        ca_df = df[
            (df["atom_name"] == "CA") &
            (df["auth_chain_id"] == chain)
        ].sort_values("auth_seq_id")

        if ca_df.empty:
            continue

        rmsd = rmsd_values.get(pdb_id, 0)
        type_label = info.get("type_full", info["type"])
        name_label = f"{pdb_id}: {info['name']} ({type_label}, RMSD: {rmsd:.2f}Å)"

        # Add backbone trace
        fig.add_trace(go.Scatter3d(
            x=ca_df["x"],
            y=ca_df["y"],
            z=ca_df["z"],
            mode="lines+markers",
            name=name_label,
            line=dict(color=color, width=3),
            marker=dict(size=2, color=color),
            hovertemplate=(
                f"{pdb_id}<br>"
                "Res: %{text}<br>"
                "x: %{x:.2f}<br>"
                "y: %{y:.2f}<br>"
                "z: %{z:.2f}<extra></extra>"
            ),
            text=ca_df["res_name3l"].astype(str) + ca_df["auth_seq_id"].astype(str),
        ))

        # Add retinal as larger markers (only from specified chain)
        ligand = info["ligand"]
        ret_df = df[
            ((df["res_name3l"] == ligand) | (df["res_name"] == ligand)) &
            (df["auth_chain_id"] == chain)
        ]

        if not ret_df.empty:
            fig.add_trace(go.Scatter3d(
                x=ret_df["x"],
                y=ret_df["y"],
                z=ret_df["z"],
                mode="markers",
                name=f"{pdb_id}: Retinal",
                marker=dict(
                    size=5,
                    color=retinal_color,
                    symbol="diamond",
                    opacity=0.8,
                ),
                hovertemplate=(
                    f"{pdb_id} Retinal<br>"
                    "Atom: %{text}<br>"
                    "x: %{x:.2f}<br>"
                    "y: %{y:.2f}<br>"
                    "z: %{z:.2f}<extra></extra>"
                ),
                text=ret_df["atom_name"],
            ))

    fig.update_layout(
        title=dict(
            text="Type I vs Type II Opsins: Different Folds, Same Ligand<br>"
                 "<sup>Retinal-based structural alignment reveals distinct 7TM topologies</sup>",
            font=dict(size=18),
        ),
        scene=dict(
            xaxis_title="X (A)",
            yaxis_title="Y (A)",
            zaxis_title="Z (A)",
            aspectmode="data",
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
        ),
        width=1000,
        height=800,
        template="plotly_white",
    )

    # Save
    html_path = output_dir / "opsin_alignment.html"
    fig.write_html(str(html_path))
    logger.info(f"  Saved visualization: {html_path}")


def create_comparison_table(
    contact_data: Dict[str, pd.DataFrame],
    rmsd_values: Dict[str, float],
    output_dir: Path,
) -> pd.DataFrame:
    """Create summary comparison table."""
    logger.info("=" * 60)
    logger.info("Step 6: Creating comparison summary")
    logger.info("=" * 60)

    records = []
    for pdb_id, info in STRUCTURE_DATASET.items():
        contacts = contact_data.get(pdb_id)
        n_contacts = len(contacts) if contacts is not None else 0

        records.append({
            "PDB ID": pdb_id,
            "Name": info["name"],
            "Type": info.get("type_full", info["type"]),
            "Organism": info["organism"],
            "Function": info["function"],
            "Chain": info["chain"],
            "Binding Residues": n_contacts,
            "RMSD (Å)": rmsd_values.get(pdb_id, float("nan")),
        })

    summary_df = pd.DataFrame(records)

    # Save CSV
    csv_path = output_dir / "opsin_comparison.csv"
    summary_df.to_csv(csv_path, index=False)
    logger.info(f"  Saved comparison: {csv_path}")

    # Log table
    logger.info("\n  Comparison Summary:")
    logger.info(summary_df.to_string(index=False))

    return summary_df


def main() -> int:
    """Run the StructureProcessor example."""
    logger.info("=" * 70)
    logger.info("STRUCTURE PROCESSOR EXAMPLE: Type I vs Type II Opsin Alignment")
    logger.info("=" * 70)
    logger.info("")
    logger.info("KEY QUESTION: How do Type I and Type II opsins compare around retinal?")
    logger.info("")
    logger.info(f"Structures ({len(STRUCTURE_DATASET)} total):")
    for pdb_id, info in STRUCTURE_DATASET.items():
        type_label = info.get("type_full", info["type"])
        logger.info(f"  {pdb_id}: {info['name']} - {type_label}")
    logger.info(f"Reference: {REFERENCE_STRUCTURE}")
    logger.info("")

    # Initialize protos
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))

    struct_proc = StructureProcessor()

    # Step 1: Download structures
    structures = download_structures()
    if not structures:
        logger.error("No structures downloaded, exiting")
        return 1

    # Step 2: Analyze retinal contacts
    contact_data = analyze_retinal_contacts(struct_proc, structures)

    # Step 3: Align on retinal
    aligned_structures, rmsd_values = align_on_retinal(structures, REFERENCE_STRUCTURE)

    if not aligned_structures:
        logger.error("No aligned structures, exiting")
        return 1

    # Step 4: Export to CIF
    export_aligned_structures(aligned_structures, OUTPUT_DIR)

    # Step 5: Visualize
    visualize_alignment(aligned_structures, rmsd_values, OUTPUT_DIR)

    # Step 6: Summary table
    create_comparison_table(contact_data, rmsd_values, OUTPUT_DIR)

    logger.info("")
    logger.info("=" * 70)
    logger.info("STRUCTURE PROCESSOR EXAMPLE COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Outputs saved to: {OUTPUT_DIR}")
    logger.info(f"Log saved to: {log_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

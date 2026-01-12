#!/usr/bin/env python3
"""GRNProcessor Example: Structural equivalence across GPCRs via GRN.

Demonstrates ProtOS capabilities:
- StructureLoader: Download GPCR structures from PDB
- SequenceProcessor: GRN annotation via annotate_with_grn()
- GRNProcessor: Load and query GRN tables
- Structure alignment via GRN positions
- CIF export of aligned structures
- PyMOL script generation for visualization

Question: "How does GRN enable structural comparison across different GPCRs?"

KEY CONCEPT: Generic Residue Numbers (GRN) define structurally equivalent
positions across GPCRs. Position 3.50 in rhodopsin is structurally equivalent
to 3.50 in ADRB2, despite different residue numbers. This enables:
- Structure-based alignment without sequence similarity
- Cross-receptor comparison of binding sites
- Transfer of functional annotations between receptors

This example aligns rhodopsin (visual) and ADRB2 (adrenergic) using their
GRN annotations, demonstrating structural equivalence despite low sequence
identity (~25%).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.grn import GRNProcessor
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "grn"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Two structurally different GPCRs to compare
STRUCTURES = {
    "rhodopsin": {
        "pdb": "1U19",
        "chain": "A",
        "ligand": "RET",  # Retinal
        "description": "Bovine rhodopsin (visual GPCR)",
        "color": "marine",
    },
    "adrb2": {
        "pdb": "2RH1",
        "chain": "A",
        "ligand": "CAU",  # Carazolol
        "description": "Human ADRB2 (adrenergic GPCR)",
        "color": "salmon",
    },
}

# Key conserved positions in GPCRs (x.50 numbering) - used for alignment
ALIGNMENT_POSITIONS = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]

# Three-letter to one-letter amino acid mapping
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def extract_sequence_from_structure(structure_df: pd.DataFrame, chain: str) -> str:
    """Extract amino acid sequence from structure dataframe."""
    chain_df = structure_df[structure_df["auth_chain_id"] == chain]
    ca_atoms = chain_df[chain_df["atom_name"] == "CA"]
    residues = ca_atoms.drop_duplicates(subset=["auth_seq_id"]).sort_values("auth_seq_id")
    sequence = "".join([THREE_TO_ONE.get(r, "X") for r in residues["res_name3l"]])
    return sequence


def get_ca_coords(structure_df: pd.DataFrame, chain: str, resnum: int) -> np.ndarray:
    """Get CA coordinates for a specific residue."""
    mask = (
        (structure_df["auth_chain_id"] == chain) &
        (structure_df["auth_seq_id"] == resnum) &
        (structure_df["atom_name"] == "CA")
    )
    row = structure_df[mask]
    if len(row) == 0:
        return None
    return np.array([row.iloc[0]["x"], row.iloc[0]["y"], row.iloc[0]["z"]])


def kabsch_align(P: np.ndarray, Q: np.ndarray) -> tuple:
    """
    Kabsch algorithm to find optimal rotation matrix.
    Returns (rotation_matrix, translation_vector) to align P onto Q.
    """
    # Center both point sets
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Rotation matrix
    R = Vt.T @ U.T

    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Translation
    t = centroid_Q - R @ centroid_P

    return R, t


def apply_transform(structure_df: pd.DataFrame, R: np.ndarray, t: np.ndarray) -> pd.DataFrame:
    """Apply rotation and translation to structure coordinates."""
    df = structure_df.copy()
    coords = df[["x", "y", "z"]].values
    transformed = (R @ coords.T).T + t
    df["x"] = transformed[:, 0]
    df["y"] = transformed[:, 1]
    df["z"] = transformed[:, 2]
    return df


def write_simple_cif(structure_df: pd.DataFrame, filepath: Path, name: str):
    """Write structure to simple CIF format."""
    lines = [
        f"data_{name}",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
    ]

    for atom_idx, (_, row) in enumerate(structure_df.iterrows(), start=1):
        group = row.get("group", "ATOM")
        atom_id = atom_idx
        element = row.get("element", "C")
        atom_name = row.get("atom_name", "CA")
        res_name = row.get("res_name3l", row.get("res_name", "UNK"))
        chain = row.get("auth_chain_id", "A")
        seq_id = row.get("auth_seq_id", 1)
        x, y, z = row["x"], row["y"], row["z"]
        occ = row.get("occupancy", 1.0)
        bfactor = row.get("b_factor", 0.0)

        lines.append(
            f"{group} {atom_id} {element} {atom_name} {res_name} {chain} {seq_id} "
            f"{x:.3f} {y:.3f} {z:.3f} {occ:.2f} {bfactor:.2f}"
        )

    lines.append("#")

    with open(filepath, "w") as f:
        f.write("\n".join(lines))


def main() -> int:
    """Run the GRNProcessor example."""
    print("=" * 70)
    print("GRN PROCESSOR EXAMPLE")
    print("Structural Equivalence Across GPCRs via GRN")
    print("=" * 70)
    print("\nKEY CONCEPT: GRN positions are structurally equivalent across GPCRs.")
    print("Position 3.50 in rhodopsin aligns with 3.50 in ADRB2.")
    print("\nStructures to compare:")
    for name, info in STRUCTURES.items():
        print(f"  {info['pdb']}: {info['description']}")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Download both GPCR structures
    # -------------------------------------------------------------------------
    print("\n[1] Downloading GPCR structures...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    structures = {}
    for name, info in STRUCTURES.items():
        pdb_id = info["pdb"]
        loader.download_and_register(pdb_id, name=f"{name}_{pdb_id}")
        df = struct_proc.load_entity(f"{name}_{pdb_id}")
        if df is not None:
            structures[name] = df
            print(f"  {pdb_id}: {len(df)} atoms")
        else:
            print(f"  {pdb_id}: Failed to load")
            return 1

    # -------------------------------------------------------------------------
    # Step 2: Filter to single chain with ligand
    # -------------------------------------------------------------------------
    print("\n[2] Filtering to single chain with ligand...")
    filtered_structures = {}

    for name, info in STRUCTURES.items():
        df = structures[name]
        chain = info["chain"]
        ligand = info["ligand"]

        # Keep only specified chain and ligand
        chain_mask = df["auth_chain_id"] == chain
        ligand_mask = df["res_name"] == ligand

        # Combine: chain atoms OR ligand atoms
        filtered = df[chain_mask | ligand_mask].copy()
        filtered_structures[name] = filtered

        n_protein = len(filtered[filtered["group"] == "ATOM"])
        n_ligand = len(filtered[ligand_mask])
        print(f"  {info['pdb']} chain {chain}: {n_protein} protein atoms, {n_ligand} ligand atoms ({ligand})")

    # -------------------------------------------------------------------------
    # Step 3: Extract sequences and annotate with GRN
    # -------------------------------------------------------------------------
    print("\n[3] Extracting sequences and annotating with GRN...")
    seq_proc = SequenceProcessor()

    sequences = {}
    for name, info in STRUCTURES.items():
        seq = extract_sequence_from_structure(filtered_structures[name], info["chain"])
        seq_name = f"{name}_{info['pdb']}"
        sequences[seq_name] = seq
        seq_proc.save_entity(seq_name, seq)
        print(f"  {seq_name}: {len(seq)} aa")

    # Save as dataset for GRN annotation
    seq_proc.save_sequences(sequences, output_file="gpcr_comparison", dataset_name="gpcr_comparison")

    # Annotate with GRN
    grn_table, summary = seq_proc.annotate_with_grn(
        dataset_name="gpcr_comparison",
        reference_table="gpcrdb_ref",
        protein_family="gpcr_a",
        output_table="gpcr_comparison_grn",
        return_summary=True,
        allow_create=True,
    )

    print(f"  Annotated {summary['global']['annotated']}/{summary['global']['total']} sequences")

    # -------------------------------------------------------------------------
    # Step 4: Get GRN positions for alignment
    # -------------------------------------------------------------------------
    print("\n[4] Identifying structurally equivalent positions (x.50)...")
    grn_proc = GRNProcessor()
    grn_table = grn_proc.load_table("gpcr_comparison_grn")

    # Build mapping of GRN -> residue number for each structure
    grn_mappings = {}
    for name, info in STRUCTURES.items():
        seq_name = f"{name}_{info['pdb']}"
        grn_mappings[name] = {}

        print(f"\n  {info['pdb']} ({info['description']}):")
        for grn_pos in ALIGNMENT_POSITIONS:
            if grn_pos in grn_table.columns:
                residue_info = grn_table.loc[seq_name, grn_pos]
                if residue_info and residue_info != "-":
                    aa = residue_info[0]
                    resnum = int(residue_info[1:])
                    grn_mappings[name][grn_pos] = {"aa": aa, "resnum": resnum}
                    print(f"    {grn_pos} -> {residue_info}")

    # -------------------------------------------------------------------------
    # Step 5: Align structures using GRN positions
    # -------------------------------------------------------------------------
    print("\n[5] Aligning structures using x.50 positions...")

    # Use rhodopsin as reference
    ref_name = "rhodopsin"
    mobile_name = "adrb2"
    ref_info = STRUCTURES[ref_name]
    mobile_info = STRUCTURES[mobile_name]

    # Collect CA coordinates for alignment positions
    ref_coords = []
    mobile_coords = []
    aligned_positions = []

    for grn_pos in ALIGNMENT_POSITIONS:
        if grn_pos in grn_mappings[ref_name] and grn_pos in grn_mappings[mobile_name]:
            ref_resnum = grn_mappings[ref_name][grn_pos]["resnum"]
            mobile_resnum = grn_mappings[mobile_name][grn_pos]["resnum"]

            ref_ca = get_ca_coords(filtered_structures[ref_name], ref_info["chain"], ref_resnum)
            mobile_ca = get_ca_coords(filtered_structures[mobile_name], mobile_info["chain"], mobile_resnum)

            if ref_ca is not None and mobile_ca is not None:
                ref_coords.append(ref_ca)
                mobile_coords.append(mobile_ca)
                aligned_positions.append(grn_pos)

    ref_coords = np.array(ref_coords)
    mobile_coords = np.array(mobile_coords)

    print(f"  Aligning on {len(aligned_positions)} positions: {', '.join(aligned_positions)}")

    # Compute and apply alignment
    R, t = kabsch_align(mobile_coords, ref_coords)
    aligned_mobile = apply_transform(filtered_structures[mobile_name], R, t)

    # Calculate RMSD after alignment
    aligned_mobile_coords = (R @ mobile_coords.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((aligned_mobile_coords - ref_coords) ** 2, axis=1)))
    print(f"  RMSD after alignment: {rmsd:.2f} A")

    # -------------------------------------------------------------------------
    # Step 6: Export aligned structures as CIF
    # -------------------------------------------------------------------------
    print("\n[6] Exporting aligned structures as CIF...")

    # Save reference structure
    ref_cif = FIGURES_DIR / f"{ref_info['pdb'].lower()}_aligned.cif"
    write_simple_cif(filtered_structures[ref_name], ref_cif, ref_info['pdb'])
    print(f"  Saved: {ref_cif.name}")

    # Save aligned mobile structure
    mobile_cif = FIGURES_DIR / f"{mobile_info['pdb'].lower()}_aligned.cif"
    write_simple_cif(aligned_mobile, mobile_cif, mobile_info['pdb'])
    print(f"  Saved: {mobile_cif.name}")

    # Save GRN table
    grn_table.to_csv(OUTPUT_DIR / "gpcr_grn_comparison.csv")
    print(f"  Saved: {OUTPUT_DIR / 'gpcr_grn_comparison.csv'}")

    # -------------------------------------------------------------------------
    # Step 7: Generate PyMOL visualization script
    # -------------------------------------------------------------------------
    print("\n[7] Generating PyMOL visualization script...")

    pymol_script = f'''# PyMOL script: GRN-based structural alignment of GPCRs
# Generated by ProtOS GRNProcessor example
#
# KEY CONCEPT: GRN positions are structurally equivalent across GPCRs.
# This script loads pre-aligned rhodopsin and ADRB2 structures.
# Despite ~25% sequence identity, the 7TM core superposes well.

# Load aligned structures from CIF files
load {ref_cif.name}, rhodopsin
load {mobile_cif.name}, adrb2

# Basic display settings
bg_color white
set ray_shadows, 0
set antialias, 2
set orthoscopic, 1

# Show both as cartoon
hide everything
show cartoon, all

# Color by receptor
color {ref_info['color']}, rhodopsin
color {mobile_info['color']}, adrb2

# Show ligands as sticks
show sticks, resn {ref_info['ligand']}
show sticks, resn {mobile_info['ligand']}
color yellow, resn {ref_info['ligand']}
color orange, resn {mobile_info['ligand']}
set stick_radius, 0.2

# Highlight x.50 positions used for alignment
'''

    # Add x.50 residue selections
    for grn_pos in ALIGNMENT_POSITIONS:
        if grn_pos in grn_mappings[ref_name] and grn_pos in grn_mappings[mobile_name]:
            ref_res = grn_mappings[ref_name][grn_pos]
            mobile_res = grn_mappings[mobile_name][grn_pos]
            pymol_script += f'''
# {grn_pos}: Rhodopsin {ref_res["aa"]}{ref_res["resnum"]}, ADRB2 {mobile_res["aa"]}{mobile_res["resnum"]}
select x50_{grn_pos.replace(".", "_")}_rho, rhodopsin and resi {ref_res["resnum"]} and name CA
select x50_{grn_pos.replace(".", "_")}_adrb, adrb2 and resi {mobile_res["resnum"]} and name CA
'''

    pymol_script += '''
# Show x.50 positions as spheres
select all_x50, x50_*
show spheres, all_x50
color cyan, all_x50 and rhodopsin
color magenta, all_x50 and adrb2
set sphere_scale, 0.5

# Create nice view
orient
zoom all, 5

# Add labels
set label_size, 14
set label_color, black
'''

    # Add labels for x.50 positions
    for grn_pos in ALIGNMENT_POSITIONS:
        if grn_pos in grn_mappings[ref_name]:
            ref_res = grn_mappings[ref_name][grn_pos]
            pymol_script += f'label rhodopsin and resi {ref_res["resnum"]} and name CA, "{grn_pos}"\n'

    pymol_script += f'''
# Summary info
print("")
print("=" * 60)
print("GRN-BASED GPCR ALIGNMENT")
print("=" * 60)
print("Reference: {ref_info['pdb']} (rhodopsin) - {ref_info['color']}")
print("Mobile: {mobile_info['pdb']} (ADRB2) - {mobile_info['color']}")
print("Aligned on: {', '.join(aligned_positions)}")
print("RMSD: {rmsd:.2f} A")
print("")
print("KEY INSIGHT: Same GRN = Same 3D position")
print("Despite ~25% sequence identity, x.50 positions superpose.")
print("=" * 60)

# Ray trace for publication
# ray 2400, 1800
# png grn_structural_equivalence.png, dpi=300
'''

    # Save PyMOL script
    pymol_path = FIGURES_DIR / "grn_structural_alignment.pml"
    with open(pymol_path, "w") as f:
        f.write(pymol_script)
    print(f"  Saved: {pymol_path}")
    print(f"  Run in PyMOL: cd {FIGURES_DIR} && pymol {pymol_path.name}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)

    print(f"\nKEY RESULT: GRN enables structural alignment across GPCRs")
    print(f"  Rhodopsin ({ref_info['pdb']}) and ADRB2 ({mobile_info['pdb']}) aligned")
    print(f"  Sequence identity: ~25%")
    print(f"  Structural RMSD (x.50 positions): {rmsd:.2f} A")

    print(f"\n  GRN position mapping:")
    for grn_pos in ALIGNMENT_POSITIONS:
        if grn_pos in grn_mappings[ref_name] and grn_pos in grn_mappings[mobile_name]:
            ref_res = grn_mappings[ref_name][grn_pos]
            mobile_res = grn_mappings[mobile_name][grn_pos]
            print(f"    {grn_pos}: Rhodopsin {ref_res['aa']}{ref_res['resnum']} <-> ADRB2 {mobile_res['aa']}{mobile_res['resnum']}")

    print(f"\n→ Same GRN = Same structural position across GPCRs")
    print(f"→ This enables cross-receptor functional annotations")

    print(f"\nOutputs:")
    print(f"  Aligned CIF: {FIGURES_DIR}")
    print(f"  GRN table: {OUTPUT_DIR / 'gpcr_grn_comparison.csv'}")
    print(f"  PyMOL script: {pymol_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

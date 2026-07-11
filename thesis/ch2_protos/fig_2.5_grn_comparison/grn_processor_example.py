#!/usr/bin/env python3
"""GRNProcessor Example: Comparing animal opsins via GRN binding pocket positions.

Demonstrates ProtOS capabilities:
- StructureLoader: Download GPCR structures from PDB
- SequenceProcessor: GRN annotation via annotate_with_grn()
- GRNProcessor: Load and query GRN tables
- Structure alignment via GRN positions
- PyMOL script generation for visualization

Question: "How does GRN enable comparison of spectral tuning sites across animal opsins?"

KEY CONCEPT: Generic Residue Numbers (GRN) define structurally equivalent positions
across GPCRs. This example compares TWO ANIMAL OPSINS (both Type II):
- Bovine rhodopsin (1U19) - vertebrate, dim light vision
- Squid rhodopsin (2Z73) - invertebrate, vision

Focus is on BINDING POCKET positions relevant to spectral tuning:
- 3.28: Counterion position (E113 in bovine)
- 3.32: Spectral tuning site
- 6.48: Rotamer switch (toggle switch)
- 7.42: Schiff base lysine (K296 in bovine)

NOTE: Microbial opsins (Type I) do NOT have GRN - this is the MOGRN gap that
LAMBDA addresses by using graph-based representations instead.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(THESIS_DIR / "shared"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.grn import GRNProcessor
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Load Color Scheme
# =============================================================================
with open(THESIS_DIR / "colorscales.yaml") as f:
    COLORS = yaml.safe_load(f)


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "grn"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Two animal opsins (both Type II GPCRs) to compare
STRUCTURES = {
    "bovine": {
        "pdb": "1U19",
        "chain": "A",
        "ligand": "RET",
        "description": "Bovine rhodopsin (vertebrate)",
        "organism": "Bos taurus",
        "uniprot": "P02699",
        "color_hex": COLORS["structures"]["1U19"]["hex"],
        "color_rgb_norm": COLORS["structures"]["1U19"]["rgb_norm"],
    },
    "squid": {
        "pdb": "2Z73",
        "chain": "A",
        "ligand": "RET",
        "description": "Squid rhodopsin (invertebrate)",
        "organism": "Todarodes pacificus",
        "uniprot": "P31356",
        "color_hex": COLORS["structures"]["2Z73"]["hex"],
        "color_rgb_norm": COLORS["structures"]["2Z73"]["rgb_norm"],
    },
}

# Retinal color
RETINAL_COLOR_HEX = COLORS["molecules"]["retinal"]["hex"]

# Conserved x.50 positions for structural alignment
ALIGNMENT_POSITIONS = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]

# Binding pocket positions for spectral tuning analysis
# Note: Schiff base lysine is at 7.42, not 7.43!
BINDING_POCKET_POSITIONS = {
    "3.28": {"function": "Counterion", "bovine_expected": "E113"},
    "3.32": {"function": "Tuning site", "bovine_expected": "A117"},
    "6.48": {"function": "Rotamer switch", "bovine_expected": "W265"},
    "7.42": {"function": "Schiff base K", "bovine_expected": "K296"},
}

# Three-letter to one-letter amino acid mapping
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def extract_sequence_from_structure(structure_df: pd.DataFrame, chain: str) -> tuple[str, int]:
    """Extract amino acid sequence from structure dataframe.

    Returns:
        tuple: (sequence, first_resnum) where first_resnum is the PDB residue number
               of the first amino acid in the sequence (needed for offset calculation)
    """
    chain_df = structure_df[structure_df["auth_chain_id"] == chain]
    ca_atoms = chain_df[chain_df["atom_name"] == "CA"]
    residues = ca_atoms.drop_duplicates(subset=["auth_seq_id"]).sort_values("auth_seq_id")
    sequence = "".join([THREE_TO_ONE.get(r, "X") for r in residues["res_name3l"]])
    first_resnum = int(residues["auth_seq_id"].iloc[0]) if len(residues) > 0 else 1
    return sequence, first_resnum


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
    """Kabsch algorithm for optimal rotation matrix. Aligns P onto Q."""
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
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
        element = row.get("element", "C")
        atom_name = row.get("atom_name", "CA")
        res_name = row.get("res_name3l", row.get("res_name", "UNK"))
        chain = row.get("auth_chain_id", "A")
        seq_id = row.get("auth_seq_id", 1)
        x, y, z = row["x"], row["y"], row["z"]
        occ = row.get("occupancy", 1.0)
        bfactor = row.get("b_factor", 0.0)

        lines.append(
            f"{group} {atom_idx} {element} {atom_name} {res_name} {chain} {seq_id} "
            f"{x:.3f} {y:.3f} {z:.3f} {occ:.2f} {bfactor:.2f}"
        )

    lines.append("#")
    with open(filepath, "w") as f:
        f.write("\n".join(lines))


def main() -> int:
    """Run the GRNProcessor example."""
    print("=" * 70)
    print("GRN PROCESSOR EXAMPLE")
    print("Comparing Animal Opsins via GRN Binding Pocket Positions")
    print("=" * 70)
    print("\nKEY QUESTION: How does GRN enable comparison of spectral tuning sites?")
    print("\nStructures to compare (both Type II / Animal Opsins):")
    for name, info in STRUCTURES.items():
        print(f"  {info['pdb']}: {info['description']} ({info['organism']})")
        print(f"         UniProt: {info['uniprot']}")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load both opsin structures
    # -------------------------------------------------------------------------
    print("\n[1] Loading opsin structures...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    structures = {}
    for name, info in STRUCTURES.items():
        pdb_id = info["pdb"]
        # Try loading by PDB ID first
        df = struct_proc.load_entity(pdb_id)
        if df is None:
            try:
                loader.download_and_register(f"pdb:{pdb_id}", name=pdb_id)
                df = struct_proc.load_entity(pdb_id)
            except Exception as e:
                print(f"  {pdb_id}: Failed to load - {e}")
                return 1

        if df is not None:
            structures[name] = df
            print(f"  {pdb_id}: {len(df)} atoms")
        else:
            print(f"  {pdb_id}: Failed to load")
            return 1

    # -------------------------------------------------------------------------
    # Step 2: Filter to single chain with ligand
    # -------------------------------------------------------------------------
    print("\n[2] Filtering to chain A with retinal...")
    filtered_structures = {}

    for name, info in STRUCTURES.items():
        df = structures[name]
        chain = info["chain"]
        ligand = info["ligand"]

        chain_mask = df["auth_chain_id"] == chain
        ligand_mask = (df["res_name"] == ligand) & (df["auth_chain_id"] == chain)
        filtered = df[chain_mask | ligand_mask].copy()
        filtered_structures[name] = filtered

        n_protein = len(filtered[filtered["group"] == "ATOM"])
        n_ligand = len(filtered[filtered["res_name"] == ligand])
        print(f"  {info['pdb']} chain {chain}: {n_protein} protein atoms, {n_ligand} ligand atoms ({ligand})")

    # -------------------------------------------------------------------------
    # Step 3: Extract sequences and annotate with GRN
    # -------------------------------------------------------------------------
    print("\n[3] Extracting sequences and annotating with GRN...")
    seq_proc = SequenceProcessor()

    sequences = {}
    first_residue_numbers = {}  # Track PDB numbering offset for each structure
    for name, info in STRUCTURES.items():
        seq, first_resnum = extract_sequence_from_structure(filtered_structures[name], info["chain"])
        seq_name = f"{name}_{info['pdb']}"
        sequences[seq_name] = seq
        first_residue_numbers[name] = first_resnum
        seq_proc.save_entity(seq_name, seq)
        offset = first_resnum - 1
        print(f"  {seq_name}: {len(seq)} aa (PDB starts at res {first_resnum}, offset={offset})")

    # Save as dataset for GRN annotation
    seq_proc.save_sequences(sequences, output_file="opsin_comparison", dataset_name="opsin_comparison")

    # Annotate with GRN
    grn_table, summary = seq_proc.annotate_with_grn(
        dataset_name="opsin_comparison",
        reference_table="gpcrdb_ref",
        protein_family="gpcr_a",
        output_table="opsin_comparison_grn",
        return_summary=True,
        allow_create=True,
    )

    print(f"  Annotated {summary['global']['annotated']}/{summary['global']['total']} sequences")

    # -------------------------------------------------------------------------
    # Step 4: Analyze GRN binding pocket positions
    # -------------------------------------------------------------------------
    print("\n[4] Analyzing binding pocket GRN positions...")
    grn_proc = GRNProcessor()
    grn_table = grn_proc.load_table("opsin_comparison_grn")

    # Build mapping of GRN -> residue number for each structure
    # IMPORTANT: GRN annotation gives sequence positions (1-indexed from extracted sequence)
    # We need to convert to PDB residue numbers using the offset
    grn_mappings = {}
    for name, info in STRUCTURES.items():
        seq_name = f"{name}_{info['pdb']}"
        grn_mappings[name] = {}

        # Calculate offset: PDB_resnum = seq_position + offset
        offset = first_residue_numbers[name] - 1

        # Get all positions (alignment + binding pocket)
        all_positions = ALIGNMENT_POSITIONS + list(BINDING_POCKET_POSITIONS.keys())
        for grn_pos in all_positions:
            if grn_pos in grn_table.columns:
                residue_info = grn_table.loc[seq_name, grn_pos]
                if residue_info and residue_info != "-":
                    aa = residue_info[0]
                    seq_resnum = int(residue_info[1:])  # Position in sequence
                    pdb_resnum = seq_resnum + offset    # Position in PDB structure
                    grn_mappings[name][grn_pos] = {
                        "aa": aa,
                        "resnum": pdb_resnum,  # Use PDB residue number
                        "seq_resnum": seq_resnum,  # Keep original for reference
                        "full": f"{aa}{pdb_resnum}",  # Update to show PDB numbering
                        "original_full": residue_info,  # Keep original GRN output
                    }

    # Print binding pocket comparison
    print("\n  Binding Pocket Position Comparison:")
    print("  " + "-" * 60)
    print(f"  {'GRN':<8} {'Function':<18} {'Bovine (1U19)':<15} {'Squid (2Z73)':<15}")
    print("  " + "-" * 60)

    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        bovine_res = grn_mappings["bovine"].get(grn_pos, {}).get("full", "-")
        squid_res = grn_mappings["squid"].get(grn_pos, {}).get("full", "-")
        print(f"  {grn_pos:<8} {pos_info['function']:<18} {bovine_res:<15} {squid_res:<15}")

    print("  " + "-" * 60)

    # -------------------------------------------------------------------------
    # Step 5: Align structures using x.50 positions
    # -------------------------------------------------------------------------
    print("\n[5] Aligning structures using x.50 positions...")

    ref_name = "bovine"
    mobile_name = "squid"
    ref_info = STRUCTURES[ref_name]
    mobile_info = STRUCTURES[mobile_name]

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

    R, t = kabsch_align(mobile_coords, ref_coords)
    aligned_mobile = apply_transform(filtered_structures[mobile_name], R, t)

    aligned_mobile_coords = (R @ mobile_coords.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((aligned_mobile_coords - ref_coords) ** 2, axis=1)))
    print(f"  RMSD after alignment: {rmsd:.2f} Å")

    # -------------------------------------------------------------------------
    # Step 6: Export aligned structures as CIF
    # -------------------------------------------------------------------------
    print("\n[6] Exporting aligned structures as CIF...")

    ref_cif = FIGURES_DIR / f"{ref_info['pdb'].lower()}_grn_aligned.cif"
    write_simple_cif(filtered_structures[ref_name], ref_cif, ref_info['pdb'])
    print(f"  Saved: {ref_cif.name}")

    mobile_cif = FIGURES_DIR / f"{mobile_info['pdb'].lower()}_grn_aligned.cif"
    write_simple_cif(aligned_mobile, mobile_cif, mobile_info['pdb'])
    print(f"  Saved: {mobile_cif.name}")

    grn_table.to_csv(OUTPUT_DIR / "opsin_grn_comparison.csv")
    print(f"  Saved: opsin_grn_comparison.csv")

    # -------------------------------------------------------------------------
    # Step 7: Generate PyMOL visualization script
    # -------------------------------------------------------------------------
    print("\n[7] Generating PyMOL visualization scripts...")

    # Get RGB values for PyMOL
    bovine_rgb = STRUCTURES["bovine"]["color_rgb_norm"]
    squid_rgb = STRUCTURES["squid"]["color_rgb_norm"]
    retinal_rgb = COLORS["molecules"]["retinal"]["rgb_norm"]

    colorscales_path = THESIS_DIR / "shared" / "colorscales_pymol.pml"
    pymol_script = f'''# PyMOL script: GRN-based comparison of animal opsins
# Generated by ProtOS GRNProcessor example
#
# Compares bovine (1U19) and squid (2Z73) rhodopsin binding pockets
# Using GRN positions to identify structurally equivalent residues

# Load shared colorscales (colors + render settings)
@{colorscales_path}

# Load aligned structures
load {ref_cif.name}, bovine
load {mobile_cif.name}, squid

# Show both as cartoon
hide everything
show cartoon, all

# Color by structure (from colorscales_pymol.pml)
color color_1U19, bovine
color color_2Z73, squid

# Show retinal as sticks
show sticks, resn RET
color retinal_rust, resn RET

'''

    # Add binding pocket residue selections
    pymol_script += "# Binding pocket positions (spectral tuning sites)\n"
    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        grn_safe = grn_pos.replace(".", "_")
        if grn_pos in grn_mappings["bovine"]:
            bovine_res = grn_mappings["bovine"][grn_pos]
            pymol_script += f"select bp_{grn_safe}_bovine, bovine and resi {bovine_res['resnum']}\n"
        if grn_pos in grn_mappings["squid"]:
            squid_res = grn_mappings["squid"][grn_pos]
            pymol_script += f"select bp_{grn_safe}_squid, squid and resi {squid_res['resnum']}\n"

    pymol_script += '''
# Show binding pocket residues as sticks
select all_bp, bp_*
show sticks, all_bp
set stick_radius, 0.12

# Highlight Schiff base lysine (7.42)
show spheres, bp_7_42_*
set sphere_scale, 0.3

# Highlight counterion (3.28)
show spheres, bp_3_28_*
set sphere_scale, 0.3

# Create nice view centered on binding pocket
orient all_bp
zoom all_bp, 8

# Labels (use dark_gray from colorscales)
set label_color, dark_gray
set label_position, [2, 2, 2]
'''

    # Add labels for binding pocket positions
    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        grn_safe = grn_pos.replace(".", "_")
        if grn_pos in grn_mappings["bovine"]:
            bovine_res = grn_mappings["bovine"][grn_pos]
            pymol_script += f'label bovine and resi {bovine_res["resnum"]} and name CA, "{grn_pos}"\n'

    pymol_script += f'''
# Summary
print("")
print("=" * 60)
print("GRN-BASED OPSIN COMPARISON")
print("=" * 60)
print("Bovine rhodopsin (1U19) vs Squid rhodopsin (2Z73)")
print("RMSD on x.50 positions: {rmsd:.2f} A")
print("")
print("BINDING POCKET POSITIONS:")
'''

    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        bovine_res = grn_mappings["bovine"].get(grn_pos, {}).get("full", "-")
        squid_res = grn_mappings["squid"].get(grn_pos, {}).get("full", "-")
        pymol_script += f'print("  {grn_pos} ({pos_info["function"]}): Bovine={bovine_res}, Squid={squid_res}")\n'

    pymol_script += '''print("")
print("KEY INSIGHT: GRN enables systematic comparison of binding pockets")
print("Same GRN position = Same structural role in spectral tuning")
print("")
print("NOTE: Type I opsins (microbial) do NOT have GRN - MOGRN gap!")
print("=" * 60)

# Ray trace for publication
# ray 2400, 1800
# png grn_opsin_comparison.png, dpi=300
'''

    pymol_path = FIGURES_DIR / "grn_structural_alignment.pml"
    with open(pymol_path, "w") as f:
        f.write(pymol_script)
    print(f"  Saved: {pymol_path.name}")

    # -------------------------------------------------------------------------
    # Step 8: Generate spectral tuning focused PyMOL script
    # -------------------------------------------------------------------------
    spectral_script = f'''# PyMOL script: Spectral tuning sites in animal opsins
# Focus on GRN positions involved in wavelength regulation

# Load shared colorscales (colors + render settings)
@{colorscales_path}

# Load structures
load {ref_cif.name}, bovine
load {mobile_cif.name}, squid

# Setup
hide everything
show cartoon, all
set cartoon_transparency, 0.7
color color_1U19, bovine
color color_2Z73, squid

# Retinal - the chromophore
show sticks, resn RET
color retinal_rust, resn RET
set stick_radius, 0.2

'''

    # Add key spectral tuning residues
    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        grn_safe = grn_pos.replace(".", "_")
        spectral_script += f"\n# {grn_pos}: {pos_info['function']}\n"
        if grn_pos in grn_mappings["bovine"]:
            bovine_res = grn_mappings["bovine"][grn_pos]
            spectral_script += f"show sticks, bovine and resi {bovine_res['resnum']}\n"
            spectral_script += f"color color_1U19, bovine and resi {bovine_res['resnum']}\n"
        if grn_pos in grn_mappings["squid"]:
            squid_res = grn_mappings["squid"][grn_pos]
            spectral_script += f"show sticks, squid and resi {squid_res['resnum']}\n"
            spectral_script += f"color color_2Z73, squid and resi {squid_res['resnum']}\n"

    spectral_script += '''
set stick_radius, 0.15

# Center on retinal
center resn RET
zoom resn RET, 12

# Ray trace
# ray 1600, 1200
# png grn_spectral_tuning.png, dpi=300
'''

    spectral_path = FIGURES_DIR / "grn_spectral_tuning.pml"
    with open(spectral_path, "w") as f:
        f.write(spectral_script)
    print(f"  Saved: {spectral_path.name}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("KEY RESULTS")
    print("=" * 70)

    print(f"\n1. STRUCTURAL ALIGNMENT:")
    print(f"   Bovine (1U19) vs Squid (2Z73) rhodopsin")
    print(f"   RMSD on x.50 positions: {rmsd:.2f} Å")

    print(f"\n2. BINDING POCKET COMPARISON (GRN positions):")
    print("   " + "-" * 50)
    print(f"   {'GRN':<8} {'Function':<18} {'Bovine':<12} {'Squid':<12}")
    print("   " + "-" * 50)
    for grn_pos, pos_info in BINDING_POCKET_POSITIONS.items():
        bovine_res = grn_mappings["bovine"].get(grn_pos, {}).get("full", "-")
        squid_res = grn_mappings["squid"].get(grn_pos, {}).get("full", "-")
        print(f"   {grn_pos:<8} {pos_info['function']:<18} {bovine_res:<12} {squid_res:<12}")
    print("   " + "-" * 50)

    print(f"\n3. THE MOGRN GAP:")
    print("   → GRN works for Type II (animal) opsins - they are GPCRs")
    print("   → GRN does NOT work for Type I (microbial) opsins - different fold!")
    print("   → LAMBDA bridges this gap using graph-based representations")

    print(f"\nOutputs:")
    print(f"  CIF files: {FIGURES_DIR}")
    print(f"  GRN table: {OUTPUT_DIR / 'opsin_grn_comparison.csv'}")
    print(f"  PyMOL: {pymol_path.name}, {spectral_path.name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

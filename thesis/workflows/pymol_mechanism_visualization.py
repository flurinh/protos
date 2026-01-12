#!/usr/bin/env python3
"""PyMOL visualization for GPCR Agonist vs Inverse Agonist Mechanism.

Creates views for each hypothesis:
- H1: Agonists bind closer to S5.43
- H2: Inverse agonists bind closer to W6.48 (toggle switch)
- H3: D2.50-W6.48 distance shortest in inverse agonists

Usage:
    pymol -cq pymol_mechanism_visualization.py

Or interactively:
    pymol
    run pymol_mechanism_visualization.py
    setup_views()
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Any

try:
    from pymol import cmd, stored
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    print("PyMOL not available. This script must be run within PyMOL.")


# Configuration
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DATA_ROOT = REPO_ROOT / "data"
STRUCTURE_DIR = DATA_ROOT / "structure" / "entities"
OUTPUT_DIR = SCRIPT_DIR / "figures"

REFERENCE_PDB = "2RH1"

# Structure metadata matching workflow
STRUCTURES: Dict[str, Dict[str, Any]] = {
    "2RH1": {"chain": "A", "ligand": "CAU", "ligand_type": "inverse_agonist", "state": "inactive", "color": "firebrick"},
    "3NY9": {"chain": "A", "ligand": "JSZ", "ligand_type": "inverse_agonist", "state": "inactive", "color": "salmon"},
    "3SN6": {"chain": "R", "ligand": "P0G", "ligand_type": "full_agonist", "state": "active", "color": "forest"},
    "4LDO": {"chain": "A", "ligand": "ALE", "ligand_type": "full_agonist", "state": "active", "color": "lime"},
    "2VT4": {"chain": "A", "ligand": "P32", "ligand_type": "antagonist", "state": "inactive", "color": "slate"},
    "2Y02": {"chain": "A", "ligand": "WHJ", "ligand_type": "full_agonist", "state": "active_like", "color": "chartreuse"},
    "2Y04": {"chain": "A", "ligand": "68H", "ligand_type": "partial_agonist", "state": "intermediate", "color": "splitpea"},
    "2Y00": {"chain": "A", "ligand": "Y00", "ligand_type": "partial_agonist", "state": "intermediate", "color": "olive"},
}

# GRN residue numbers for each structure (from workflow JSON)
# Format: pdb_id -> {grn -> residue_number}
GRN_RESIDUES = {
    "2RH1": {"5.43": 203, "6.48": 286, "2.50": 79},
    "3NY9": {"5.43": 203, "6.48": 286, "2.50": 79},
    "3SN6": {"5.43": 203, "6.48": 286, "2.50": 79},
    "4LDO": {"5.43": 1203, "6.48": 1286, "2.50": 1079},
    "2VT4": {"5.43": 211, "6.48": 294, "2.50": 87},
    "2Y02": {"5.43": 211, "6.48": 294, "2.50": 87},
    "2Y04": {"5.43": 211, "6.48": 294, "2.50": 87},
    "2Y00": {"5.43": 211, "6.48": 294, "2.50": 87},
}

LIGAND_COLORS = {
    "inverse_agonist": "red",
    "antagonist": "orange",
    "partial_agonist": "yellow",
    "full_agonist": "green",
}


def find_structure_file(pdb_id: str) -> Path:
    """Find structure file for a PDB ID."""
    entity_dir = STRUCTURE_DIR / pdb_id.lower()
    if entity_dir.exists():
        for ext in [".cif", ".cif.gz", ".pdb", ".pdb.gz"]:
            for name in [pdb_id, pdb_id.lower()]:
                filepath = entity_dir / f"{name}{ext}"
                if filepath.exists():
                    return filepath
    # Try mmcif directory
    mmcif_dir = DATA_ROOT / "structure" / "mmcif"
    for ext in [".cif", ".cif.gz"]:
        filepath = mmcif_dir / f"{pdb_id.lower()}{ext}"
        if filepath.exists():
            return filepath
    return None


def load_structures() -> List[str]:
    """Load all structures into PyMOL."""
    loaded = []
    for pdb_id in STRUCTURES.keys():
        filepath = find_structure_file(pdb_id)
        if filepath:
            print(f"Loading {pdb_id} from {filepath}")
            cmd.load(str(filepath), pdb_id)
            loaded.append(pdb_id)
        else:
            print(f"Fetching {pdb_id} from PDB...")
            try:
                cmd.fetch(pdb_id, pdb_id, type="cif")
                loaded.append(pdb_id)
            except Exception as e:
                print(f"  Failed: {e}")
    return loaded


def extract_gpcr_chains(loaded: List[str]) -> List[str]:
    """Extract GPCR chains with ligands."""
    extracted = []
    for pdb_id in loaded:
        meta = STRUCTURES[pdb_id]
        chain = meta["chain"]
        obj_name = f"{pdb_id}_gpcr"

        protein_sel = f"{pdb_id} and chain {chain} and polymer.protein"
        ligand_sel = f"{pdb_id} and hetatm and not resn HOH and (within 8 of ({protein_sel}))"

        cmd.select("_temp", f"({protein_sel}) or ({ligand_sel})")
        cmd.create(obj_name, "_temp")
        cmd.delete("_temp")
        cmd.delete(pdb_id)

        extracted.append(pdb_id)
        print(f"  Created {obj_name}")
    return extracted


def align_structures(extracted: List[str]):
    """Align all structures to reference."""
    ref_obj = f"{REFERENCE_PDB}_gpcr"
    print(f"\nAligning to {REFERENCE_PDB}...")

    for pdb_id in extracted:
        if pdb_id == REFERENCE_PDB:
            continue
        obj_name = f"{pdb_id}_gpcr"
        result = cmd.super(f"{obj_name} and name CA", f"{ref_obj} and name CA")
        print(f"  {pdb_id}: RMSD = {result[0]:.2f} A")


def setup_base_representations():
    """Set up base visual representations."""
    cmd.hide("everything")

    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        # Protein cartoon
        cmd.show("cartoon", f"{obj_name} and polymer.protein")
        cmd.color(meta["color"], f"{obj_name} and polymer.protein")
        cmd.set("cartoon_transparency", 0.7, obj_name)

        # Ligand sticks
        ligand = meta["ligand"]
        ligand_sel = f"{obj_name} and resn {ligand}"
        cmd.show("sticks", ligand_sel)
        cmd.color(LIGAND_COLORS.get(meta["ligand_type"], "white"), ligand_sel)
        cmd.set("stick_radius", 0.25, ligand_sel)


def create_grn_selections():
    """Create selections for key GRN positions."""
    for grn in ["5.43", "6.48", "2.50"]:
        grn_safe = grn.replace(".", "_")
        sel_parts = []

        for pdb_id, meta in STRUCTURES.items():
            obj_name = f"{pdb_id}_gpcr"
            if obj_name not in cmd.get_object_list():
                continue

            chain = meta["chain"]
            if pdb_id in GRN_RESIDUES and grn in GRN_RESIDUES[pdb_id]:
                res_id = GRN_RESIDUES[pdb_id][grn]
                sel_parts.append(f"({obj_name} and chain {chain} and resi {res_id})")

        if sel_parts:
            cmd.select(f"grn_{grn_safe}", " or ".join(sel_parts))
            cmd.show("sticks", f"grn_{grn_safe} and sidechain")

    print("Created GRN selections: grn_5_43, grn_6_48, grn_2_50")


def setup_hypothesis_view(hypothesis: str, pdb_ids: List[str], grn_positions: List[str], title: str):
    """Set up a view for a specific hypothesis."""
    # Disable all
    cmd.disable("all")

    # Enable selected structures
    for pdb_id in pdb_ids:
        obj_name = f"{pdb_id}_gpcr"
        if obj_name in cmd.get_object_list():
            cmd.enable(obj_name)

    # Highlight GRN positions
    for grn in grn_positions:
        grn_safe = grn.replace(".", "_")
        cmd.color("magenta", f"grn_{grn_safe} and sidechain")
        cmd.set("stick_radius", 0.35, f"grn_{grn_safe}")

    # Create distance measurements
    dist_count = 0
    for pdb_id in pdb_ids:
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        meta = STRUCTURES[pdb_id]
        chain = meta["chain"]
        ligand = meta["ligand"]

        for grn in grn_positions:
            if pdb_id not in GRN_RESIDUES or grn not in GRN_RESIDUES[pdb_id]:
                continue

            res_id = GRN_RESIDUES[pdb_id][grn]
            dist_name = f"{hypothesis}_{pdb_id}_{grn.replace('.', '_')}"

            # Distance from ligand to GRN residue
            lig_sel = f"{obj_name} and resn {ligand} and not elem H"
            res_sel = f"{obj_name} and chain {chain} and resi {res_id} and sidechain and not elem H"

            cmd.distance(dist_name, lig_sel, res_sel)
            dist_count += 1

    # Style distances
    cmd.set("dash_color", "yellow")
    cmd.set("dash_width", 2.0)
    cmd.set("dash_gap", 0.3)
    cmd.hide("labels", f"{hypothesis}_*")

    # Center on ligands
    ligand_sels = [f"{pdb}_gpcr and resn {STRUCTURES[pdb]['ligand']}" for pdb in pdb_ids if f"{pdb}_gpcr" in cmd.get_object_list()]
    if ligand_sels:
        cmd.center(" or ".join(ligand_sels))
        cmd.zoom(" or ".join(ligand_sels), buffer=8)

    # Store scene
    cmd.scene(hypothesis, "store")
    print(f"Created scene '{hypothesis}': {title}")


def setup_h1_view():
    """H1: Agonists bind closer to S5.43."""
    # Compare agonist (green) vs inverse agonist (red)
    agonists = ["3SN6", "4LDO", "2Y02"]
    inverse = ["2RH1", "3NY9"]

    setup_hypothesis_view("H1_S543", agonists + inverse, ["5.43"], "S5.43 Contact (Agonist vs Inverse)")


def setup_h2_view():
    """H2: Inverse agonists bind closer to W6.48."""
    agonists = ["3SN6", "4LDO", "2Y02"]
    inverse = ["2RH1", "3NY9"]

    setup_hypothesis_view("H2_W648", agonists + inverse, ["6.48"], "W6.48 Toggle Switch (Inverse vs Agonist)")


def setup_h3_view():
    """H3: D2.50-W6.48 distance shortest in inverse agonists."""
    # Show D2.50 to W6.48 distances
    cmd.disable("all")

    agonists = ["3SN6", "4LDO"]
    inverse = ["2RH1", "3NY9"]
    all_pdbs = agonists + inverse

    for pdb_id in all_pdbs:
        obj_name = f"{pdb_id}_gpcr"
        if obj_name in cmd.get_object_list():
            cmd.enable(obj_name)

    # Highlight D2.50 and W6.48
    cmd.color("cyan", "grn_2_50 and sidechain")
    cmd.color("magenta", "grn_6_48 and sidechain")
    cmd.set("stick_radius", 0.35, "grn_2_50 or grn_6_48")

    # Create D2.50 to W6.48 distances
    for pdb_id in all_pdbs:
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        meta = STRUCTURES[pdb_id]
        chain = meta["chain"]

        if pdb_id not in GRN_RESIDUES:
            continue

        d250_res = GRN_RESIDUES[pdb_id].get("2.50")
        w648_res = GRN_RESIDUES[pdb_id].get("6.48")

        if d250_res and w648_res:
            dist_name = f"H3_{pdb_id}_D250_W648"
            d250_sel = f"{obj_name} and chain {chain} and resi {d250_res} and sidechain and not elem H"
            w648_sel = f"{obj_name} and chain {chain} and resi {w648_res} and sidechain and not elem H"
            cmd.distance(dist_name, d250_sel, w648_sel)

    cmd.set("dash_color", "orange")
    cmd.set("dash_width", 2.5)
    cmd.hide("labels", "H3_*")

    # Center
    cmd.center("grn_2_50 or grn_6_48")
    cmd.zoom("grn_2_50 or grn_6_48", buffer=10)

    cmd.scene("H3_D250_W648", "store")
    print("Created scene 'H3_D250_W648': D2.50-W6.48 Distance (Mechanism Marker)")


def setup_overview():
    """Overview showing all structures."""
    cmd.enable("all")
    setup_base_representations()
    cmd.zoom("polymer.protein")
    cmd.scene("overview", "store")
    print("Created scene 'overview': All structures aligned")


def create_ligand_type_groups():
    """Create groups by ligand type."""
    agonist_objs = []
    inverse_objs = []

    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        if meta["ligand_type"] in ["full_agonist", "partial_agonist"]:
            agonist_objs.append(obj_name)
        elif meta["ligand_type"] == "inverse_agonist":
            inverse_objs.append(obj_name)

    if agonist_objs:
        cmd.group("agonists", " ".join(agonist_objs))
    if inverse_objs:
        cmd.group("inverse_agonists", " ".join(inverse_objs))


def render_hypothesis_figures():
    """Render figures for each hypothesis."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.bg_color("white")

    for scene in ["H1_S543", "H2_W648", "H3_D250_W648"]:
        cmd.scene(scene, "recall")
        output_path = OUTPUT_DIR / f"pymol_{scene.lower()}.png"
        cmd.png(str(output_path), width=1800, height=1400, dpi=300, ray=1)
        print(f"Rendered: {output_path}")


def setup_views():
    """Main setup function."""
    print("=" * 60)
    print("GPCR MECHANISM VISUALIZATION")
    print("Agonist vs Inverse Agonist Binding Mode Comparison")
    print("=" * 60)

    # Load and prepare structures
    print("\n[1/6] Loading structures...")
    loaded = load_structures()
    print(f"  Loaded {len(loaded)} structures")

    print("\n[2/6] Extracting GPCR chains...")
    extracted = extract_gpcr_chains(loaded)

    print("\n[3/6] Aligning structures...")
    align_structures(extracted)

    print("\n[4/6] Setting up representations...")
    setup_base_representations()
    create_grn_selections()
    create_ligand_type_groups()

    print("\n[5/6] Creating hypothesis views...")
    setup_overview()
    setup_h1_view()
    setup_h2_view()
    setup_h3_view()

    print("\n[6/6] Final setup...")
    cmd.bg_color("white")
    cmd.scene("overview", "recall")

    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print("\nAvailable scenes:")
    print("  scene overview, recall     - All structures aligned")
    print("  scene H1_S543, recall      - S5.43 contact (agonist closer)")
    print("  scene H2_W648, recall      - W6.48 toggle (inverse closer)")
    print("  scene H3_D250_W648, recall - D2.50-W6.48 distance")
    print("\nGroups:")
    print("  enable agonists            - Show agonist-bound")
    print("  enable inverse_agonists    - Show inverse agonist-bound")
    print("\nTo render figures:")
    print("  render_hypothesis_figures()")


# Convenience functions
def show_agonists():
    cmd.disable("inverse_agonists")
    cmd.enable("agonists")

def show_inverse():
    cmd.disable("agonists")
    cmd.enable("inverse_agonists")

def show_all():
    cmd.enable("all")


if __name__ == "__main__" or "pymol" in sys.modules:
    if PYMOL_AVAILABLE:
        setup_views()
    else:
        print("Run this script within PyMOL:")
        print("  pymol -cq pymol_mechanism_visualization.py")

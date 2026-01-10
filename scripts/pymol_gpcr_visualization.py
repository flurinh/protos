#!/usr/bin/env python3
"""PyMOL visualization script for GPCR binding pocket comparison.

This script creates comprehensive views of aligned GPCR structures showing:
1. All structures aligned to a reference
2. Binding pocket comparison (agonist vs antagonist)
3. Water networks and bridging waters
4. Active vs inactive state comparison

Usage:
    pymol -cq pymol_gpcr_visualization.py

Or interactively:
    pymol
    run pymol_gpcr_visualization.py
    setup_views()
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# Try to import pymol - handle both standalone and embedded modes
try:
    from pymol import cmd, stored
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    print("PyMOL not available. This script must be run within PyMOL.")


# Configuration
DATA_ROOT = Path(__file__).resolve().parents[1] / "data"
STRUCTURE_DIR = DATA_ROOT / "structure" / "entities"
OUTPUT_DIR = DATA_ROOT / "visualizations"

# Reference structure for alignment
REFERENCE_PDB = "2RH1"

# Structure metadata
STRUCTURES: Dict[str, Dict[str, Any]] = {
    "2RH1": {
        "receptor": "beta2_AR",
        "chain": "A",
        "ligand": "CAU",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "color": "marine",  # Blue for inactive
    },
    "3NY9": {
        "receptor": "beta2_AR",
        "chain": "A",
        "ligand": "CLR",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "color": "slate",
    },
    "3SN6": {
        "receptor": "beta2_AR",
        "chain": "R",
        "ligand": "P0G",
        "ligand_type": "full_agonist",
        "state": "active",
        "color": "forest",  # Green for active
    },
    "4LDO": {
        "receptor": "beta2_AR",
        "chain": "A",
        "ligand": "ALE",
        "ligand_type": "full_agonist",
        "state": "active",
        "color": "lime",
    },
    "2VT4": {
        "receptor": "beta1_AR",
        "chain": "A",
        "ligand": "P32",
        "ligand_type": "antagonist",
        "state": "inactive",
        "color": "deepteal",
    },
    "2Y03": {
        "receptor": "beta1_AR",
        "chain": "A",
        "ligand": "Y01",
        "ligand_type": "full_agonist",
        "state": "active_like",
        "color": "chartreuse",
    },
    "2Y04": {
        "receptor": "beta1_AR",
        "chain": "A",
        "ligand": "Y01",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "color": "yellow",
    },
    "2Y00": {
        "receptor": "beta1_AR",
        "chain": "A",
        "ligand": "Y01",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "color": "orange",
    },
}

# Color schemes
STATE_COLORS = {
    "active": "green",
    "active_like": "chartreuse",
    "intermediate": "yellow",
    "inactive": "blue",
}

LIGAND_TYPE_COLORS = {
    "full_agonist": "green",
    "partial_agonist": "yellow",
    "antagonist": "orange",
    "inverse_agonist": "red",
}


def find_structure_file(pdb_id: str) -> Optional[Path]:
    """Find the structure file for a given PDB ID."""
    entity_dir = STRUCTURE_DIR / pdb_id
    if not entity_dir.exists():
        # Try lowercase
        entity_dir = STRUCTURE_DIR / pdb_id.lower()

    if entity_dir.exists():
        # Look for mmCIF or PDB file
        for ext in [".cif", ".cif.gz", ".pdb", ".pdb.gz"]:
            filepath = entity_dir / f"{pdb_id}{ext}"
            if filepath.exists():
                return filepath
            filepath = entity_dir / f"{pdb_id.lower()}{ext}"
            if filepath.exists():
                return filepath

    return None


def load_structures() -> List[str]:
    """Load all GPCR structures into PyMOL."""
    loaded = []

    for pdb_id, meta in STRUCTURES.items():
        filepath = find_structure_file(pdb_id)

        if filepath:
            print(f"Loading {pdb_id} from {filepath}")
            cmd.load(str(filepath), pdb_id)
            loaded.append(pdb_id)
        else:
            # Try fetching from PDB
            print(f"Fetching {pdb_id} from PDB...")
            try:
                cmd.fetch(pdb_id, pdb_id, type="cif")
                loaded.append(pdb_id)
            except Exception as e:
                print(f"  Failed to load {pdb_id}: {e}")

    return loaded


def extract_gpcr_chain(pdb_id: str, chain_id: str, water_cutoff: float = 5.0) -> str:
    """Extract the GPCR chain with nearby waters and ligand.

    Args:
        pdb_id: PDB identifier
        chain_id: Chain ID of the GPCR
        water_cutoff: Distance cutoff for including waters (Angstroms)

    Returns:
        Name of the extracted object
    """
    obj_name = f"{pdb_id}_gpcr"

    # Select the GPCR chain (protein only)
    protein_sel = f"{pdb_id} and chain {chain_id} and polymer.protein"

    # Select waters near the protein chain
    water_sel = f"{pdb_id} and resn HOH and (within {water_cutoff} of ({protein_sel}))"

    # Select all ligands (non-protein, non-water) near the chain
    ligand_sel = f"{pdb_id} and hetatm and not resn HOH and (within {water_cutoff} of ({protein_sel}))"

    # Create the combined selection
    cmd.select("_temp_gpcr", f"({protein_sel}) or ({water_sel}) or ({ligand_sel})")
    cmd.create(obj_name, "_temp_gpcr")
    cmd.delete("_temp_gpcr")

    return obj_name


def align_to_reference(mobile_obj: str, reference_obj: str) -> Tuple[float, int]:
    """Align mobile structure to reference using CA atoms.

    Returns:
        Tuple of (RMSD, number of aligned atoms)
    """
    # Use super for sequence-independent structural alignment
    result = cmd.super(
        f"{mobile_obj} and name CA",
        f"{reference_obj} and name CA",
        object=f"{mobile_obj}_aln"
    )

    # Clean up alignment object
    cmd.delete(f"{mobile_obj}_aln")

    return result[0], result[1]  # RMSD, aligned atoms


def setup_representations():
    """Set up visual representations for all structures."""
    # Hide everything first
    cmd.hide("everything")

    # For each structure
    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"

        if obj_name not in cmd.get_object_list():
            continue

        # Protein: cartoon representation
        cmd.show("cartoon", f"{obj_name} and polymer.protein")
        cmd.color(meta["color"], f"{obj_name} and polymer.protein")

        # Ligand: sticks
        ligand_sel = f"{obj_name} and hetatm and not resn HOH"
        cmd.show("sticks", ligand_sel)
        cmd.color(LIGAND_TYPE_COLORS.get(meta["ligand_type"], "white"), ligand_sel)
        cmd.set("stick_radius", 0.25, ligand_sel)

        # Waters: spheres (small)
        water_sel = f"{obj_name} and resn HOH"
        cmd.show("spheres", water_sel)
        cmd.color("cyan", water_sel)
        cmd.set("sphere_scale", 0.15, water_sel)


def create_binding_pocket_view():
    """Create a view focused on the binding pocket."""
    # Select binding pocket residues (common GRN positions)
    # TM3: 3.28-3.37, TM5: 5.39-5.46, TM6: 6.48-6.55, TM7: 7.35-7.42

    # Center on ligands
    cmd.center("hetatm and not resn HOH")
    cmd.zoom("hetatm and not resn HOH", buffer=8)

    # Set nice orientation
    cmd.turn("x", 20)
    cmd.turn("y", -30)


def create_water_network_view():
    """Create a view highlighting water networks."""
    # Show waters more prominently
    for pdb_id in STRUCTURES.keys():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        water_sel = f"{obj_name} and resn HOH"
        cmd.set("sphere_scale", 0.25, water_sel)
        cmd.set("sphere_transparency", 0.3, water_sel)


def create_state_comparison_groups():
    """Create groups for active vs inactive comparison."""
    active_objs = []
    inactive_objs = []

    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        state = meta["state"]
        if state in ["active", "active_like"]:
            active_objs.append(obj_name)
        else:
            inactive_objs.append(obj_name)

    # Create groups
    if active_objs:
        cmd.group("active_states", " ".join(active_objs))
    if inactive_objs:
        cmd.group("inactive_states", " ".join(inactive_objs))

    return active_objs, inactive_objs


def create_ligand_type_groups():
    """Create groups for agonist vs antagonist comparison."""
    agonist_objs = []
    antagonist_objs = []

    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        ligand_type = meta["ligand_type"]
        if ligand_type in ["full_agonist", "partial_agonist"]:
            agonist_objs.append(obj_name)
        else:
            antagonist_objs.append(obj_name)

    # Create groups
    if agonist_objs:
        cmd.group("agonist_bound", " ".join(agonist_objs))
    if antagonist_objs:
        cmd.group("antagonist_bound", " ".join(antagonist_objs))

    return agonist_objs, antagonist_objs


def show_tm6_movement():
    """Highlight TM6 to show the activation-related movement."""
    # TM6 is approximately residues around position 6.30-6.60
    # We'll color TM6 differently to highlight the movement

    for pdb_id, meta in STRUCTURES.items():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        # Create TM6 selection (approximate - would need GRN mapping for precise)
        # For now, select by secondary structure in the relevant region
        chain = meta["chain"]

        # Color TM6 region distinctly
        state = meta["state"]
        if state in ["active", "active_like"]:
            cmd.color("red", f"{obj_name} and chain {chain} and resi 265-295")
        else:
            cmd.color("blue", f"{obj_name} and chain {chain} and resi 265-295")


def save_views():
    """Save predefined views for easy switching."""
    # View 1: Overall alignment
    cmd.zoom("polymer.protein")
    cmd.view("overall", "store")

    # View 2: Binding pocket
    create_binding_pocket_view()
    cmd.view("binding_pocket", "store")

    # View 3: Extracellular view (looking down into pocket)
    cmd.zoom("hetatm and not resn HOH", buffer=12)
    cmd.turn("x", 90)
    cmd.view("extracellular", "store")

    # View 4: Intracellular view (G-protein interface)
    cmd.zoom("polymer.protein")
    cmd.turn("x", -90)
    cmd.view("intracellular", "store")

    # View 5: Side view (membrane plane)
    cmd.zoom("polymer.protein")
    cmd.turn("x", 0)
    cmd.turn("y", 0)
    cmd.view("side", "store")

    print("\nStored views: overall, binding_pocket, extracellular, intracellular, side")
    print("Use: cmd.view('view_name', 'recall') to switch views")


def create_session_file():
    """Save a PyMOL session file."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    session_path = OUTPUT_DIR / "gpcr_binding_pocket_comparison.pse"
    cmd.save(str(session_path))
    print(f"Session saved to: {session_path}")


def render_images():
    """Render publication-quality images."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Set rendering options
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.bg_color("white")

    views = ["overall", "binding_pocket", "extracellular", "side"]

    for view_name in views:
        cmd.view(view_name, "recall")
        output_path = OUTPUT_DIR / f"gpcr_{view_name}.png"
        cmd.png(str(output_path), width=2400, height=2400, dpi=300, ray=1)
        print(f"Rendered: {output_path}")


def setup_views():
    """Main function to set up all visualizations."""
    print("=" * 60)
    print("GPCR Binding Pocket Visualization Setup")
    print("=" * 60)

    # Step 1: Load structures
    print("\n[1/6] Loading structures...")
    loaded = load_structures()
    print(f"  Loaded {len(loaded)} structures")

    # Step 2: Extract GPCR chains
    print("\n[2/6] Extracting GPCR chains with waters and ligands...")
    extracted = []
    for pdb_id in loaded:
        meta = STRUCTURES[pdb_id]
        obj_name = extract_gpcr_chain(pdb_id, meta["chain"])
        extracted.append(obj_name)
        print(f"  Created {obj_name}")

    # Delete original full structures to reduce clutter
    for pdb_id in loaded:
        cmd.delete(pdb_id)

    # Step 3: Align all structures to reference
    print(f"\n[3/6] Aligning all structures to {REFERENCE_PDB}...")
    ref_obj = f"{REFERENCE_PDB}_gpcr"

    for obj_name in extracted:
        if obj_name == ref_obj:
            continue
        rmsd, n_atoms = align_to_reference(obj_name, ref_obj)
        print(f"  {obj_name}: RMSD = {rmsd:.2f} A ({n_atoms} atoms)")

    # Step 4: Set up representations
    print("\n[4/6] Setting up visual representations...")
    setup_representations()

    # Step 5: Create groups
    print("\n[5/6] Creating comparison groups...")
    active, inactive = create_state_comparison_groups()
    print(f"  Active states: {len(active)} structures")
    print(f"  Inactive states: {len(inactive)} structures")

    agonist, antagonist = create_ligand_type_groups()
    print(f"  Agonist-bound: {len(agonist)} structures")
    print(f"  Antagonist-bound: {len(antagonist)} structures")

    # Step 6: Save views
    print("\n[6/6] Saving predefined views...")
    save_views()

    # Final setup
    cmd.view("binding_pocket", "recall")
    cmd.set("cartoon_transparency", 0.7)
    cmd.bg_color("white")

    print("\n" + "=" * 60)
    print("VISUALIZATION SETUP COMPLETE")
    print("=" * 60)
    print("\nUseful commands:")
    print("  cmd.view('binding_pocket', 'recall')  - Binding pocket view")
    print("  cmd.view('extracellular', 'recall')   - Top-down view")
    print("  cmd.enable('active_states')           - Show active structures")
    print("  cmd.disable('inactive_states')        - Hide inactive structures")
    print("  cmd.enable('agonist_bound')           - Show agonist-bound")
    print("  cmd.disable('antagonist_bound')       - Hide antagonist-bound")
    print("  show_only_active()                    - Show only active state")
    print("  show_only_inactive()                  - Show only inactive state")
    print("  show_water_comparison()               - Compare water networks")
    print("  save_session()                        - Save PyMOL session")
    print("  render_all()                          - Render images")


# Convenience functions for interactive use
def show_only_active():
    """Show only active state structures."""
    cmd.disable("inactive_states")
    cmd.enable("active_states")


def show_only_inactive():
    """Show only inactive state structures."""
    cmd.disable("active_states")
    cmd.enable("inactive_states")


def show_only_agonist():
    """Show only agonist-bound structures."""
    cmd.disable("antagonist_bound")
    cmd.enable("agonist_bound")


def show_only_antagonist():
    """Show only antagonist-bound structures."""
    cmd.disable("agonist_bound")
    cmd.enable("antagonist_bound")


def show_water_comparison():
    """Set up view for water network comparison."""
    # Make protein more transparent
    cmd.set("cartoon_transparency", 0.85)

    # Make waters larger and more visible
    for pdb_id in STRUCTURES.keys():
        obj_name = f"{pdb_id}_gpcr"
        if obj_name not in cmd.get_object_list():
            continue

        water_sel = f"{obj_name} and resn HOH"
        cmd.set("sphere_scale", 0.3, water_sel)
        cmd.set("sphere_transparency", 0.0, water_sel)

    cmd.view("binding_pocket", "recall")
    print("Water comparison view active. Waters are now more prominent.")


def show_ligand_comparison():
    """Set up view for ligand comparison."""
    # Hide protein
    cmd.hide("cartoon")

    # Show ligands as large sticks
    cmd.show("sticks", "hetatm and not resn HOH")
    cmd.set("stick_radius", 0.3)

    # Center on ligands
    cmd.center("hetatm and not resn HOH")
    cmd.zoom("hetatm and not resn HOH", buffer=3)

    print("Ligand comparison view active. Only ligands visible.")


def reset_view():
    """Reset to default view."""
    cmd.show("cartoon", "polymer.protein")
    cmd.set("cartoon_transparency", 0.7)
    cmd.set("stick_radius", 0.25)
    cmd.view("binding_pocket", "recall")


def save_session():
    """Save PyMOL session."""
    create_session_file()


def render_all():
    """Render all standard views."""
    render_images()


def highlight_grn_position(grn: str):
    """Highlight a specific GRN position across all structures.

    Note: This requires GRN annotations to be loaded from the structure data.
    For now, this is a placeholder that would need integration with the
    StructureProcessor GRN data.
    """
    print(f"GRN position highlighting requires loading GRN annotations.")
    print(f"Use StructureProcessor.annotate_with_grn() results for precise GRN selection.")


# Run setup if executed as script
if __name__ == "__main__" or "pymol" in sys.modules:
    if PYMOL_AVAILABLE:
        setup_views()
    else:
        print("This script must be run within PyMOL.")
        print("Usage: pymol -cq pymol_gpcr_visualization.py")

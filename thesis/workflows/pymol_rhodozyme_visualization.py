#!/usr/bin/env python3
"""PyMOL Visualization for Rhodozyme Design.

Visualizes the catalytic triad grafting:
- Trypsin: Original Ser-His-Asp catalytic triad
- Rhodopsin: Target positions for grafting
- Overlay comparison

Run with: pymol -cq pymol_rhodozyme_visualization.py
Or interactively: run pymol_rhodozyme_visualization.py
"""
from __future__ import annotations

import json
from pathlib import Path

# PyMOL commands
from pymol import cmd

# =============================================================================
# Configuration
# =============================================================================
THESIS_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = THESIS_DIR / "outputs" / "rhodozyme"
FIGURES_DIR = THESIS_DIR / "workflows" / "figures"

# Load design data
DESIGN_FILE = OUTPUT_DIR / "rhodozyme_design.json"

# PDB IDs
RHODOPSIN_PDB = "3PQR"
TRYPSIN_PDB = "1S81"

# Trypsin catalytic triad (chymotrypsin numbering)
TRYPSIN_TRIAD = {
    "ser": 195,
    "his": 57,
    "asp": 102,
}


def setup_visualization():
    """Set up PyMOL visualization settings."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.2)


def load_structures():
    """Load rhodopsin and trypsin structures."""
    print("Loading structures...")
    cmd.fetch(RHODOPSIN_PDB, name="rhodopsin", type="pdb")
    cmd.fetch(TRYPSIN_PDB, name="trypsin", type="pdb")

    # Clean up - keep only chain A for both
    cmd.remove("rhodopsin and not chain A")
    cmd.remove("trypsin and not chain A")


def load_design_data() -> dict:
    """Load rhodozyme design data."""
    if DESIGN_FILE.exists():
        with open(DESIGN_FILE) as f:
            return json.load(f)
    else:
        # Fallback to hardcoded values
        return {
            "grafted_positions": {
                "ser_site": {"position": 84},
                "his_site": {"position": 55},
                "asp_site": {"position": 300},
            }
        }


def visualize_trypsin_triad():
    """Visualize the trypsin catalytic triad."""
    print("Visualizing trypsin catalytic triad...")

    # Show trypsin as cartoon
    cmd.show("cartoon", "trypsin")
    cmd.color("palegreen", "trypsin")

    # Highlight catalytic triad
    triad_sel = f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']} and chain A"
    cmd.show("sticks", triad_sel)

    # Color by role
    cmd.color("red", f"trypsin and resi {TRYPSIN_TRIAD['ser']} and chain A")  # Ser - nucleophile
    cmd.color("blue", f"trypsin and resi {TRYPSIN_TRIAD['his']} and chain A")  # His - base
    cmd.color("orange", f"trypsin and resi {TRYPSIN_TRIAD['asp']} and chain A")  # Asp - electrostatic

    # Show CA atoms as spheres
    cmd.show("spheres", f"{triad_sel} and name CA")
    cmd.set("sphere_scale", 0.5, f"{triad_sel} and name CA")

    # Label
    cmd.label(f"trypsin and resi {TRYPSIN_TRIAD['ser']} and name CA and chain A", '"Ser195"')
    cmd.label(f"trypsin and resi {TRYPSIN_TRIAD['his']} and name CA and chain A", '"His57"')
    cmd.label(f"trypsin and resi {TRYPSIN_TRIAD['asp']} and name CA and chain A", '"Asp102"')

    # Draw distance lines between CA atoms
    cmd.distance("ser_his_tryp",
                 f"trypsin and resi {TRYPSIN_TRIAD['ser']} and name CA and chain A",
                 f"trypsin and resi {TRYPSIN_TRIAD['his']} and name CA and chain A")
    cmd.distance("his_asp_tryp",
                 f"trypsin and resi {TRYPSIN_TRIAD['his']} and name CA and chain A",
                 f"trypsin and resi {TRYPSIN_TRIAD['asp']} and name CA and chain A")
    cmd.distance("ser_asp_tryp",
                 f"trypsin and resi {TRYPSIN_TRIAD['ser']} and name CA and chain A",
                 f"trypsin and resi {TRYPSIN_TRIAD['asp']} and name CA and chain A")

    cmd.color("green", "ser_his_tryp")
    cmd.color("green", "his_asp_tryp")
    cmd.color("green", "ser_asp_tryp")


def visualize_rhodopsin_target(design_data: dict):
    """Visualize the target positions in rhodopsin."""
    print("Visualizing rhodopsin target positions...")

    positions = design_data["grafted_positions"]
    ser_pos = positions["ser_site"]["position"]
    his_pos = positions["his_site"]["position"]
    asp_pos = positions["asp_site"]["position"]

    # Show rhodopsin as cartoon
    cmd.show("cartoon", "rhodopsin")
    cmd.color("lightblue", "rhodopsin")

    # Highlight target positions
    target_sel = f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and chain A"
    cmd.show("sticks", target_sel)

    # Color by role (same as trypsin)
    cmd.color("red", f"rhodopsin and resi {ser_pos} and chain A")  # Ser site
    cmd.color("blue", f"rhodopsin and resi {his_pos} and chain A")  # His site
    cmd.color("orange", f"rhodopsin and resi {asp_pos} and chain A")  # Asp site

    # Show CA atoms as spheres
    cmd.show("spheres", f"{target_sel} and name CA")
    cmd.set("sphere_scale", 0.5, f"{target_sel} and name CA")

    # Label with mutation info
    ser_orig = positions["ser_site"].get("original_residue", "X")
    his_orig = positions["his_site"].get("original_residue", "X")
    asp_orig = positions["asp_site"].get("original_residue", "X")

    cmd.label(f"rhodopsin and resi {ser_pos} and name CA and chain A", f'"{ser_orig}{ser_pos}->Ser"')
    cmd.label(f"rhodopsin and resi {his_pos} and name CA and chain A", f'"{his_orig}{his_pos}->His"')
    cmd.label(f"rhodopsin and resi {asp_pos} and name CA and chain A", f'"{asp_orig}{asp_pos}->Asp"')

    # Draw distance lines between target CA atoms
    cmd.distance("ser_his_rhod",
                 f"rhodopsin and resi {ser_pos} and name CA and chain A",
                 f"rhodopsin and resi {his_pos} and name CA and chain A")
    cmd.distance("his_asp_rhod",
                 f"rhodopsin and resi {his_pos} and name CA and chain A",
                 f"rhodopsin and resi {asp_pos} and name CA and chain A")
    cmd.distance("ser_asp_rhod",
                 f"rhodopsin and resi {ser_pos} and name CA and chain A",
                 f"rhodopsin and resi {asp_pos} and name CA and chain A")

    cmd.color("magenta", "ser_his_rhod")
    cmd.color("magenta", "his_asp_rhod")
    cmd.color("magenta", "ser_asp_rhod")


def create_scenes(design_data: dict):
    """Create PyMOL scenes for different views."""
    positions = design_data["grafted_positions"]
    ser_pos = positions["ser_site"]["position"]
    his_pos = positions["his_site"]["position"]
    asp_pos = positions["asp_site"]["position"]

    # Scene 1: Trypsin catalytic triad
    print("Creating Scene 1: Trypsin triad...")
    cmd.hide("everything")
    cmd.show("cartoon", "trypsin")
    cmd.show("sticks", f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']} and chain A")
    cmd.show("spheres", f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']} and name CA and chain A")
    cmd.show("labels")
    # Enable distance objects (they are objects, not selections)
    cmd.enable("ser_his_tryp")
    cmd.enable("his_asp_tryp")
    cmd.enable("ser_asp_tryp")
    cmd.zoom(f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']}", buffer=10)
    cmd.scene("trypsin_triad", "store")

    # Scene 2: Rhodopsin target positions
    print("Creating Scene 2: Rhodopsin target...")
    cmd.hide("everything")
    cmd.show("cartoon", "rhodopsin")
    cmd.show("sticks", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and chain A")
    cmd.show("spheres", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and name CA and chain A")
    cmd.show("labels")
    # Enable distance objects
    cmd.enable("ser_his_rhod")
    cmd.enable("his_asp_rhod")
    cmd.enable("ser_asp_rhod")
    cmd.zoom(f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos}", buffer=10)
    cmd.scene("rhodopsin_target", "store")

    # Scene 3: Full rhodopsin with highlighted ICL region
    print("Creating Scene 3: Full rhodopsin...")
    cmd.hide("everything")
    cmd.show("cartoon", "rhodopsin")
    cmd.color("lightblue", "rhodopsin")
    cmd.color("yellow", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos}")
    cmd.show("spheres", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and name CA and chain A")
    cmd.zoom("rhodopsin")
    cmd.scene("rhodopsin_full", "store")

    # Scene 4: Overlay comparison (CA atoms only)
    print("Creating Scene 4: CA overlay...")
    cmd.hide("everything")

    # Create pseudoatom markers at CA positions for comparison
    # Trypsin triad CAs
    cmd.show("spheres", f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']} and name CA and chain A")
    cmd.set("sphere_scale", 0.6, "trypsin and name CA")
    cmd.color("green", f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']} and name CA")

    # Rhodopsin target CAs
    cmd.show("spheres", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and name CA and chain A")
    cmd.set("sphere_scale", 0.6, "rhodopsin and name CA")
    cmd.color("magenta", f"rhodopsin and resi {ser_pos}+{his_pos}+{asp_pos} and name CA")

    # Enable all distance objects
    cmd.enable("ser_his_tryp")
    cmd.enable("his_asp_tryp")
    cmd.enable("ser_asp_tryp")
    cmd.enable("ser_his_rhod")
    cmd.enable("his_asp_rhod")
    cmd.enable("ser_asp_rhod")
    cmd.show("labels")
    cmd.zoom(f"trypsin and resi {TRYPSIN_TRIAD['ser']}+{TRYPSIN_TRIAD['his']}+{TRYPSIN_TRIAD['asp']}", buffer=15)
    cmd.scene("ca_overlay", "store")


def save_images():
    """Save PNG images of each scene."""
    print("Saving images...")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    cmd.scene("trypsin_triad", "recall")
    cmd.ray(1200, 900)
    cmd.png(str(FIGURES_DIR / "rhodozyme_trypsin_triad.png"))

    cmd.scene("rhodopsin_target", "recall")
    cmd.ray(1200, 900)
    cmd.png(str(FIGURES_DIR / "rhodozyme_rhodopsin_target.png"))

    cmd.scene("rhodopsin_full", "recall")
    cmd.ray(1200, 900)
    cmd.png(str(FIGURES_DIR / "rhodozyme_rhodopsin_full.png"))

    print(f"Images saved to {FIGURES_DIR}")


def save_session():
    """Save PyMOL session."""
    session_path = OUTPUT_DIR / "rhodozyme_visualization.pse"
    cmd.save(str(session_path))
    print(f"Session saved to {session_path}")


def main():
    """Main visualization routine."""
    print("=" * 60)
    print("RHODOZYME DESIGN - PyMOL VISUALIZATION")
    print("=" * 60)

    # Setup
    setup_visualization()

    # Load structures
    load_structures()

    # Load design data
    design_data = load_design_data()
    print(f"Design data loaded from {DESIGN_FILE}")

    # Visualize
    visualize_trypsin_triad()
    visualize_rhodopsin_target(design_data)

    # Create scenes
    create_scenes(design_data)

    # Save outputs
    save_images()
    save_session()

    # Final view - show rhodopsin target
    cmd.scene("rhodopsin_target", "recall")

    print("\n" + "=" * 60)
    print("VISUALIZATION COMPLETE")
    print("=" * 60)
    print("Scenes available:")
    print("  - trypsin_triad: Original catalytic triad")
    print("  - rhodopsin_target: Target positions for grafting")
    print("  - rhodopsin_full: Full rhodopsin structure")
    print("  - ca_overlay: CA atom comparison")
    print("\nUse: cmd.scene('scene_name', 'recall') to switch views")


if __name__ == "__main__":
    main()

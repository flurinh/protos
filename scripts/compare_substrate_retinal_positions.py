#!/usr/bin/env python
"""
Compare substrate and retinal positions between before_refolding and after_refolding structures.
"""

import gemmi
import numpy as np


def get_atom_position(structure, chain_id, res_num, atom_name):
    """Get position of a specific atom."""
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                for residue in chain:
                    if residue.seqid.num == res_num:
                        for atom in residue:
                            if atom.name == atom_name:
                                return np.array([atom.pos.x, atom.pos.y, atom.pos.z])
    return None


def get_chain_center_of_mass(structure, chain_id):
    """Calculate center of mass for an entire chain."""
    positions = []
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                for residue in chain:
                    for atom in residue:
                        positions.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if positions:
        return np.mean(positions, axis=0)
    return None


def get_residue_center_of_mass(structure, chain_id, res_num):
    """Calculate center of mass for a specific residue."""
    positions = []
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                for residue in chain:
                    if residue.seqid.num == res_num:
                        for atom in residue:
                            positions.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if positions:
        return np.mean(positions, axis=0)
    return None


def get_catalytic_triad_center(structure, chain_id):
    """Get center of catalytic triad (His226, Asp229, Ser138) using CA atoms."""
    his226_ca = get_atom_position(structure, chain_id, 226, "CA")
    asp229_ca = get_atom_position(structure, chain_id, 229, "CA")
    ser138_ca = get_atom_position(structure, chain_id, 138, "CA")

    positions = [p for p in [his226_ca, asp229_ca, ser138_ca] if p is not None]
    if len(positions) == 3:
        return np.mean(positions, axis=0)
    return None


def distance(p1, p2):
    """Calculate Euclidean distance between two points."""
    if p1 is None or p2 is None:
        return None
    return np.linalg.norm(p1 - p2)


def analyze_structure(cif_path, label):
    """Analyze a single structure and return distances."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {label}")
    print(f"File: {cif_path}")
    print(f"{'='*60}")

    structure = gemmi.read_structure(cif_path)

    # 1. Distance from Ser138 OG (chain A) to chain B center of mass (substrate)
    ser138_og = get_atom_position(structure, "A", 138, "OG")
    chain_b_com = get_chain_center_of_mass(structure, "B")

    if ser138_og is not None:
        print(f"  Ser138 OG position: ({ser138_og[0]:.3f}, {ser138_og[1]:.3f}, {ser138_og[2]:.3f})")
    else:
        print("  WARNING: Ser138 OG not found!")

    if chain_b_com is not None:
        print(f"  Chain B (substrate) COM: ({chain_b_com[0]:.3f}, {chain_b_com[1]:.3f}, {chain_b_com[2]:.3f})")
    else:
        print("  WARNING: Chain B not found!")

    dist_ser138_substrate = distance(ser138_og, chain_b_com)

    # 2. Distance from Lys296 NZ (chain A) to RET/retinal center of mass (residue 327 in chain A)
    lys296_nz = get_atom_position(structure, "A", 296, "NZ")
    ret_com = get_residue_center_of_mass(structure, "A", 327)

    if lys296_nz is not None:
        print(f"  Lys296 NZ position: ({lys296_nz[0]:.3f}, {lys296_nz[1]:.3f}, {lys296_nz[2]:.3f})")
    else:
        print("  WARNING: Lys296 NZ not found!")

    if ret_com is not None:
        print(f"  RET (res 327) COM: ({ret_com[0]:.3f}, {ret_com[1]:.3f}, {ret_com[2]:.3f})")
    else:
        print("  WARNING: Retinal (residue 327) not found!")

    dist_lys296_retinal = distance(lys296_nz, ret_com)

    # 3. Distance from catalytic triad center to substrate center
    triad_center = get_catalytic_triad_center(structure, "A")

    if triad_center is not None:
        print(f"  Catalytic triad center: ({triad_center[0]:.3f}, {triad_center[1]:.3f}, {triad_center[2]:.3f})")
    else:
        print("  WARNING: Could not compute catalytic triad center!")

    dist_triad_substrate = distance(triad_center, chain_b_com)

    return {
        "ser138_substrate": dist_ser138_substrate,
        "lys296_retinal": dist_lys296_retinal,
        "triad_substrate": dist_triad_substrate,
    }


def main():
    before_path = "/data/fast/projects/protos/data/models/boltzgen/design_job/predictions/final_ranked_designs/final_30_designs/before_refolding/rank1_config_6.cif"
    after_path = "/data/fast/projects/protos/data/models/boltzgen/design_job/predictions/final_ranked_designs/final_30_designs/rank1_config_6.cif"

    before_distances = analyze_structure(before_path, "BEFORE REFOLDING")
    after_distances = analyze_structure(after_path, "AFTER REFOLDING")

    # Print comparison
    print("\n")
    print("=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    print(f"{'Measurement':<45} {'Before':>10} {'After':>10} {'Delta':>10}")
    print("-" * 80)

    measurements = [
        ("Ser138 OG to Substrate (Chain B) COM", "ser138_substrate"),
        ("Lys296 NZ to Retinal (Res 327) COM", "lys296_retinal"),
        ("Catalytic Triad Center to Substrate COM", "triad_substrate"),
    ]

    for label, key in measurements:
        before_val = before_distances[key]
        after_val = after_distances[key]

        if before_val is not None and after_val is not None:
            delta = after_val - before_val
            print(f"{label:<45} {before_val:>10.3f} {after_val:>10.3f} {delta:>+10.3f}")
        else:
            before_str = f"{before_val:.3f}" if before_val is not None else "N/A"
            after_str = f"{after_val:.3f}" if after_val is not None else "N/A"
            print(f"{label:<45} {before_str:>10} {after_str:>10} {'N/A':>10}")

    print("-" * 80)
    print("\nNote: All distances are in Angstroms (Å)")
    print("      Positive delta = distance increased after refolding")
    print("      Negative delta = distance decreased after refolding")


if __name__ == "__main__":
    main()

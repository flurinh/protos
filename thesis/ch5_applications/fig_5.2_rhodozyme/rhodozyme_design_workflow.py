#!/usr/bin/env python3
"""Rhodozyme Design Workflow - Light-Activated Enzyme Engineering.

Biological Question:
Can we graft a catalytic triad onto rhodopsin's conformationally-active
intracellular domain to create a light-switchable enzyme?

Approach:
1. Load active rhodopsin structure (open intracellular domain)
2. Extract trypsin catalytic triad geometry (Ser-His-Asp)
3. Find matching positions in rhodopsin's ICL region
4. Generate chimeric sequence with grafted triad
5. Predict structure with Boltz
6. Validate catalytic geometry is maintained

Demonstrates cross-processor composition:
- StructureProcessor: Load and analyze structures
- GRNProcessor: Identify key positions (7.50)
- SequenceProcessor: Generate chimeric sequences
- ModelManager: Boltz structure prediction
- GraphProcessor: Contact analysis

Question: "Can we engineer a light-controlled protease?"
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "rhodozyme"
FIGURES_DIR = THESIS_DIR / "workflows" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.io.ingest.structure_loader import StructureLoader
from protos.analysis.structure.alignment import kabsch_alignment, calculate_rmsd
from protos.analysis.structure.geometry import calculate_distance, calculate_distance_matrix


# Amino acid 3-letter to 1-letter code mapping
AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


# =============================================================================
# Configuration
# =============================================================================
# Active rhodopsin - "open" intracellular domain
RHODOPSIN_PDB = "3PQR"  # Metarhodopsin II (active state)
RHODOPSIN_CHAIN = "A"

# Geometric parameters
ICL_DISTANCE_FROM_750 = 15.0  # Angstroms from position 7.50 (fallback if no GRN)
DISTANCE_TOLERANCE = 2.0  # Angstrom tolerance for triangle matching
TOP_N_MATCHES = 10  # Number of top designs to generate per enzyme

# =============================================================================
# GRN-based Cytoplasmic Region Definition
# =============================================================================
# TM helices alternate orientation through membrane:
#   Odd TMs (1,3,5,7): N-term extracellular, C-term cytoplasmic → cyto = END (high GRN)
#   Even TMs (2,4,6):  N-term cytoplasmic, C-term extracellular → cyto = START (low GRN)
#
# Cytoplasmic face = TM1_end + ICL1 + TM2_start + TM3_end + ICL2 + TM4_start +
#                    TM5_end + ICL3 + TM6_start + TM7_end + H8

# Standard cytoplasmic boundaries (can be adjusted ±4 for helix bends)
CYTOPLASMIC_BOUNDARIES = {
    # Odd TMs: cytoplasmic END (high GRN positions)
    "TM1_end": {"helix": 1, "boundary": 60, "direction": ">="},   # 1.60 is cytoplasmic end
    "TM3_end": {"helix": 3, "boundary": 55, "direction": ">="},   # 3.55 is cytoplasmic end
    "TM5_end": {"helix": 5, "boundary": 68, "direction": ">="},   # 5.68 is cytoplasmic end
    "TM7_end": {"helix": 7, "boundary": 53, "direction": ">="},   # 7.53+ is cytoplasmic (after NPxxY)

    # Even TMs: cytoplasmic START (low GRN positions)
    "TM2_start": {"helix": 2, "boundary": 42, "direction": "<="},  # 2.42- is cytoplasmic start
    "TM4_start": {"helix": 4, "boundary": 42, "direction": "<="},  # 4.42- is cytoplasmic start
    "TM6_start": {"helix": 6, "boundary": 35, "direction": "<="},  # 6.35- is cytoplasmic start
}

# Intracellular loops (always cytoplasmic)
ICL_REGIONS = {
    "ICL1": 12,  # 12.xx - between TM1 and TM2
    "ICL2": 34,  # 34.xx - between TM3 and TM4
    "ICL3": 56,  # 56.xx - between TM5 and TM6 (main G-protein binding)
    "H8": 8,     # 8.xx - cytoplasmic helix after TM7
}

# For backward compatibility, build the old-style dict
ICL_GRN_REGIONS = {
    # Intracellular loops (always included)
    "ICL1": (12, None),
    "ICL2": (34, None),
    "ICL3": (56, None),
    "H8": (8, None),
    # TM cytoplasmic regions
    "TM1_end": (1, lambda pos: pos >= 56),    # End of TM1 (high positions)
    "TM2_start": (2, lambda pos: pos <= 46),  # Start of TM2 (low positions)
    "TM3_end": (3, lambda pos: pos >= 51),    # End of TM3 (high positions)
    "TM4_start": (4, lambda pos: pos <= 46),  # Start of TM4 (low positions)
    "TM5_end": (5, lambda pos: pos >= 64),    # End of TM5 (high positions)
    "TM6_start": (6, lambda pos: pos <= 39),  # Start of TM6 (low positions)
    "TM7_end": (7, lambda pos: pos >= 53),    # End of TM7 (high positions, after NPxxY)
}

# Activation-dependent helix constraint
# TM6 moves ~10-14A outward during GPCR activation - ideal for light-switchable enzyme
# Require at least one catalytic residue on TM6 to leverage this conformational change
REQUIRE_DYNAMIC_HELIX = True  # Require at least one residue on TM5 or TM6
DYNAMIC_HELICES = {5, 6}  # TM5 and TM6 move significantly during activation

# Retinal cofactor - include in all designs to block the retinal binding pocket
# This prevents substrate from docking in the wrong site
# SMILES without the carbonyl group (Schiff base linkage to Lys)
RETINAL_SMILES = "CC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C"
RETINAL_ID = "RET"

# =============================================================================
# Reference Enzymes with 3-residue active sites
# Each enzyme has: PDB, chain, triad residues, substrate SMILES, short ID
# =============================================================================
REFERENCE_ENZYMES = {
    "trypsin": {
        "id": "TRP",
        "pdb": "1S81",
        "chain": "A",
        "description": "Bovine trypsin (serine protease)",
        "triad": {
            "res1": {"res_num": 195, "res_name": "SER", "role": "nucleophile"},
            "res2": {"res_num": 57, "res_name": "HIS", "role": "base"},
            "res3": {"res_num": 102, "res_name": "ASP", "role": "electrostatic"},
        },
        "substrate_smiles": "NC(=N)c1ccccc1",  # Benzamidine
        "substrate_name": "benzamidine",
    },
    "chymotrypsin": {
        "id": "CHY",
        "pdb": "4CHA",
        "chain": "B",  # Main chain for His57, Asp102; Ser195 is in chain C
        "description": "Bovine chymotrypsin (serine protease)",
        "triad": {
            # Chymotrypsin is cleaved into chains B+C; catalytic residues span both
            "res1": {"res_num": 195, "res_name": "SER", "role": "nucleophile", "chain": "C"},
            "res2": {"res_num": 57, "res_name": "HIS", "role": "base", "chain": "B"},
            "res3": {"res_num": 102, "res_name": "ASP", "role": "electrostatic", "chain": "B"},
        },
        "substrate_smiles": "CC(=O)Nc1ccc(cc1)C(=O)O",  # N-acetyl-p-aminobenzoic acid
        "substrate_name": "NAPABA",
    },
    "papain": {
        "id": "PAP",
        "pdb": "9PAP",
        "chain": "A",
        "description": "Papain (cysteine protease)",
        "triad": {
            "res1": {"res_num": 25, "res_name": "CYS", "role": "nucleophile"},
            "res2": {"res_num": 159, "res_name": "HIS", "role": "base"},
            "res3": {"res_num": 175, "res_name": "ASN", "role": "electrostatic"},
        },
        "substrate_smiles": "NC(Cc1ccccc1)C(=O)O",  # Phenylalanine
        "substrate_name": "phenylalanine",
    },
    "subtilisin": {
        "id": "SUB",
        "pdb": "1SBC",
        "chain": "A",
        "description": "Subtilisin Carlsberg (serine protease, different fold)",
        "triad": {
            "res1": {"res_num": 221, "res_name": "SER", "role": "nucleophile"},
            "res2": {"res_num": 64, "res_name": "HIS", "role": "base"},
            "res3": {"res_num": 32, "res_name": "ASP", "role": "electrostatic"},
        },
        "substrate_smiles": "CC(C)CC(N)C(=O)O",  # Leucine
        "substrate_name": "leucine",
    },
}


def is_cytoplasmic_residue(grn: str, boundary_offset: int = 0) -> bool:
    """Check if a residue is on the cytoplasmic face based on its GRN.

    Args:
        grn: Generic residue number (e.g., "5.68", "56.50")
        boundary_offset: Adjustment to TM helix boundaries (-4 to +4).
            Positive = include more of the helix (expand into membrane)
            Negative = include less of the helix (shrink toward cytoplasm)

    Returns True if the residue is in an ICL loop, H8, or the cytoplasmic
    end/start of a TM helix.
    """
    if pd.isna(grn) or grn == "":
        return False

    try:
        parts = grn.replace("x", ".").split(".")
        helix = int(parts[0])
        pos = int(parts[1]) if len(parts) > 1 else 50

        # Check ICL loops first (always cytoplasmic)
        if helix in ICL_REGIONS.values():
            return True

        # Check TM helix boundaries with offset
        # Odd TMs (1,3,5,7): cytoplasmic = END (high positions, use >=)
        # Even TMs (2,4,6): cytoplasmic = START (low positions, use <=)
        for region_name, config in CYTOPLASMIC_BOUNDARIES.items():
            if helix == config["helix"]:
                adjusted_boundary = config["boundary"]

                # Apply offset: for >=, positive offset lowers threshold (more inclusive)
                #               for <=, positive offset raises threshold (more inclusive)
                if config["direction"] == ">=":
                    adjusted_boundary -= boundary_offset
                    return pos >= adjusted_boundary
                else:  # "<="
                    adjusted_boundary += boundary_offset
                    return pos <= adjusted_boundary

        return False
    except (ValueError, IndexError):
        return False


def get_cytoplasmic_residues_for_offset(
    grn_df: pd.DataFrame,
    boundary_offset: int = 0
) -> pd.DataFrame:
    """Get all cytoplasmic-facing residues for a given boundary offset.

    Args:
        grn_df: DataFrame with GRN annotations
        boundary_offset: Adjustment to helix boundaries (-4 to +4)

    Returns:
        DataFrame of cytoplasmic residues with their GRN annotations
    """
    grn_unique = grn_df[["auth_seq_id", "res_name", "grn"]].drop_duplicates()
    grn_unique["is_cyto"] = grn_unique["grn"].apply(
        lambda g: is_cytoplasmic_residue(g, boundary_offset)
    )
    return grn_unique[grn_unique["is_cyto"]]


def get_helix_number(grn: str) -> Optional[int]:
    """Extract the helix number from a GRN string.

    Returns the helix/region number (1-8 for TM helices, 12/34/56 for loops),
    or None if the GRN is invalid.
    """
    if pd.isna(grn) or grn == "":
        return None

    try:
        parts = grn.replace("x", ".").split(".")
        return int(parts[0])
    except (ValueError, IndexError):
        return None


def is_on_dynamic_helix(grn: str) -> bool:
    """Check if a residue is on a helix that moves during GPCR activation.

    TM5 and TM6 undergo significant movement during activation:
    - TM6: ~10-14 Angstrom outward movement at cytoplasmic end
    - TM5: ~3-5 Angstrom outward movement

    Having catalytic residues on these helices enables light-controlled
    enzyme activity through conformational change.
    """
    helix = get_helix_number(grn)
    return helix in DYNAMIC_HELICES if helix is not None else False


def get_ca_coordinates(df: pd.DataFrame, res_num: int, chain: str = None) -> Optional[np.ndarray]:
    """Extract CA coordinates for a residue."""
    mask = (df["auth_seq_id"] == res_num) & (df["atom_name"] == "CA")
    if chain:
        mask &= df["auth_chain_id"] == chain
    ca = df[mask]
    if len(ca) == 0:
        return None
    return ca[["x", "y", "z"]].values[0]


def get_cb_coordinates(df: pd.DataFrame, res_num: int, chain: str = None) -> Optional[np.ndarray]:
    """Extract CB coordinates for a residue (or CA for glycine)."""
    mask = (df["auth_seq_id"] == res_num) & (df["atom_name"] == "CB")
    if chain:
        mask &= df["auth_chain_id"] == chain
    cb = df[mask]
    if len(cb) == 0:
        # Glycine has no CB, use CA
        return get_ca_coordinates(df, res_num, chain)
    return cb[["x", "y", "z"]].values[0]


def get_backbone_atom(df: pd.DataFrame, res_num: int, atom_name: str, chain: str = None) -> Optional[np.ndarray]:
    """Extract coordinates for a backbone atom (N, CA, C, O)."""
    mask = (df["auth_seq_id"] == res_num) & (df["atom_name"] == atom_name)
    if chain:
        mask &= df["auth_chain_id"] == chain
    atoms = df[mask]
    if len(atoms) == 0:
        return None
    return atoms[["x", "y", "z"]].values[0]


def compute_ideal_cb_direction(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Compute ideal CB direction from backbone N, CA, C atoms.

    Uses standard tetrahedral geometry: CB is placed opposite the N-C bisector,
    tilted by the tetrahedral angle (~109.5°) out of the backbone plane.
    This gives the expected sidechain direction for any residue, including glycine.
    """
    # Normalized vectors from CA to N and C
    ca_n = (n - ca) / np.linalg.norm(n - ca)
    ca_c = (c - ca) / np.linalg.norm(c - ca)

    # Normal to the N-CA-C plane (right-hand rule gives L-amino acid chirality)
    plane_normal = np.cross(ca_n, ca_c)
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Bisector of N-CA-C angle (points "into" the backbone)
    bisector = ca_n + ca_c
    bisector = bisector / np.linalg.norm(bisector)

    # CB direction: opposite to bisector, tilted out of plane
    # For tetrahedral geometry: ~54.75° from plane normal, ~125.25° from bisector
    # Coefficients: sin(54.75°) ≈ 0.816, cos(54.75°) ≈ 0.578
    cb_direction = -bisector * 0.578 + plane_normal * 0.816
    cb_direction = cb_direction / np.linalg.norm(cb_direction)

    return cb_direction


def get_sidechain_vector(df: pd.DataFrame, res_num: int, chain: str = None) -> Optional[np.ndarray]:
    """Get normalized CA→CB vector (side chain direction).

    For glycine (no CB), computes the ideal CB direction from backbone geometry.
    """
    ca = get_ca_coordinates(df, res_num, chain)
    if ca is None:
        return None

    cb = get_cb_coordinates(df, res_num, chain)

    # If CB exists and is different from CA, use actual CB
    if cb is not None:
        vec = cb - ca
        norm = np.linalg.norm(vec)
        if norm > 0.01:
            return vec / norm

    # For glycine or missing CB: compute ideal direction from backbone
    n = get_backbone_atom(df, res_num, "N", chain)
    c = get_backbone_atom(df, res_num, "C", chain)

    if n is not None and c is not None:
        return compute_ideal_cb_direction(n, ca, c)

    return None


def compute_triangle_distances(coords: List[np.ndarray]) -> Tuple[float, float, float]:
    """Compute pairwise distances for 3 points (order: 0-1, 1-2, 0-2).

    Uses protos.analysis.structure.geometry.calculate_distance.
    """
    d01 = float(calculate_distance(coords[0], coords[1]))
    d12 = float(calculate_distance(coords[1], coords[2]))
    d02 = float(calculate_distance(coords[0], coords[2]))
    return d01, d12, d02


def compute_cb_direction_errors(
    source_ca: List[np.ndarray],
    source_cb_vecs: List[np.ndarray],
    target_ca: List[np.ndarray],
    target_cb_vecs: List[np.ndarray],
) -> Tuple[List[float], float]:
    """Compute CB direction errors after aligning CA triangles.

    Uses protos.analysis.structure.alignment.kabsch_alignment to:
    1. Compute rotation R that aligns source CA triangle to target CA triangle
    2. Apply R to source CB vectors to get "expected" directions
    3. Compare expected vs actual target CB vectors

    Returns: (list of angle errors in degrees, RMSD of angles)
    """
    # Stack coordinates (ensure float64 dtype for SVD)
    src_pts = np.array(source_ca, dtype=np.float64)
    tgt_pts = np.array(target_ca, dtype=np.float64)

    # Get rotation from source to target using ProtOS kabsch_alignment
    # kabsch_alignment returns (rotation, translation, rmsd)
    R, _, _ = kabsch_alignment(tgt_pts, src_pts)

    # Rotate source CB vectors to target frame
    angle_errors = []
    for src_cb, tgt_cb in zip(source_cb_vecs, target_cb_vecs):
        if src_cb is None or tgt_cb is None:
            continue
        # Apply rotation to source CB vector
        expected_cb = R @ src_cb
        expected_cb = expected_cb / np.linalg.norm(expected_cb)

        # Angle between expected and actual
        cos_angle = np.clip(np.dot(expected_cb, tgt_cb), -1, 1)
        angle_err = np.degrees(np.arccos(cos_angle))
        angle_errors.append(angle_err)

    if angle_errors:
        rmsd = np.sqrt(np.mean([e**2 for e in angle_errors]))
    else:
        rmsd = 0.0

    return angle_errors, rmsd


def find_top_triplets(
    residue_data: List[Tuple[int, str, np.ndarray, Optional[np.ndarray], str]],  # (res_num, res_name, ca, cb_vec, grn)
    source_ca: List[np.ndarray],  # Trypsin CA coordinates [Ser, His, Asp]
    source_cb_vecs: List[np.ndarray],  # Trypsin CB vectors [Ser, His, Asp]
    target_distances: Tuple[float, float, float],
    distance_tolerance: float = 2.0,
    angle_tolerance: float = 30.0,  # degrees - for CB direction error after alignment
    top_n: int = 10,
    require_dynamic_helix: bool = False,  # Require at least one residue on TM5/TM6
) -> List[Dict]:
    """Find top N triplets of residues matching target geometry via exhaustive search.

    Three-stage filtering:
    1. Dynamic helix filter (optional): At least one residue must be on TM5 or TM6
       to leverage conformational change during activation
    2. Distance filter: CA-CA triangle must match within tolerance
    3. CB direction filter: After Kabsch alignment of CA triangles, compare
       rotated source CB vectors to candidate CB vectors

    This ensures sidechains point in the same relative directions as in the enzyme.
    Returns list of top N matches sorted by score.
    """
    import heapq

    # Use a max-heap (negate scores) to keep top N
    top_matches = []  # List of (-score, match_dict)
    candidates_passing_distance = 0
    candidates_passing_both = 0
    candidates_rejected_no_dynamic = 0

    n = len(residue_data)
    if n < 3:
        return []

    total_triplets = n * (n - 1) * (n - 2) // 6

    # Iterate through all triplets
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                res_i, name_i, ca_i, cb_vec_i, grn_i = residue_data[i]
                res_j, name_j, ca_j, cb_vec_j, grn_j = residue_data[j]
                res_k, name_k, ca_k, cb_vec_k, grn_k = residue_data[k]

                # Stage 0: Dynamic helix filter (if enabled)
                if require_dynamic_helix:
                    has_dynamic = (
                        is_on_dynamic_helix(grn_i) or
                        is_on_dynamic_helix(grn_j) or
                        is_on_dynamic_helix(grn_k)
                    )
                    if not has_dynamic:
                        candidates_rejected_no_dynamic += 1
                        continue

                positions = [
                    (res_i, name_i, ca_i, cb_vec_i, grn_i),
                    (res_j, name_j, ca_j, cb_vec_j, grn_j),
                    (res_k, name_k, ca_k, cb_vec_k, grn_k),
                ]

                # Try all 6 permutations (role assignments: Ser, His, Asp)
                for perm in itertools.permutations(positions):
                    cand_ca = [p[2] for p in perm]
                    cand_cb_vecs = [p[3] for p in perm]

                    # Stage 1: Distance filter
                    cand_dists = compute_triangle_distances(cand_ca)
                    dist_errors = [abs(c - t) for c, t in zip(cand_dists, target_distances)]
                    dist_rmsd = np.sqrt(np.mean([e**2 for e in dist_errors]))

                    if dist_rmsd > distance_tolerance:
                        continue

                    candidates_passing_distance += 1

                    # Stage 2: CB direction filter (after Kabsch alignment)
                    cb_errors, cb_direction_rmsd = compute_cb_direction_errors(
                        source_ca, source_cb_vecs,
                        cand_ca, cand_cb_vecs
                    )

                    if cb_direction_rmsd > angle_tolerance:
                        continue

                    candidates_passing_both += 1

                    # Combined score: distance RMSD + weighted CB direction error
                    score = dist_rmsd + 0.1 * cb_direction_rmsd

                    match = {
                        "ser_pos": perm[0][0],
                        "ser_name": perm[0][1],
                        "ser_cb_vec": perm[0][3],
                        "ser_grn": perm[0][4],
                        "his_pos": perm[1][0],
                        "his_name": perm[1][1],
                        "his_cb_vec": perm[1][3],
                        "his_grn": perm[1][4],
                        "asp_pos": perm[2][0],
                        "asp_name": perm[2][1],
                        "asp_cb_vec": perm[2][3],
                        "asp_grn": perm[2][4],
                        "total_score": score,
                        "distance_rmsd": dist_rmsd,
                        "cb_direction_rmsd": cb_direction_rmsd,
                        "cb_direction_errors": cb_errors,
                        "distances": cand_dists,
                        "has_dynamic_helix": any(is_on_dynamic_helix(p[4]) for p in perm),
                    }

                    # Keep top N using heap
                    if len(top_matches) < top_n:
                        heapq.heappush(top_matches, (-score, id(match), match))
                    elif score < -top_matches[0][0]:
                        heapq.heapreplace(top_matches, (-score, id(match), match))

    # Extract matches sorted by score (best first)
    results = sorted([m[2] for m in top_matches], key=lambda x: x["total_score"])

    # Add metadata to all results
    for match in results:
        match["candidates_passing_distance"] = candidates_passing_distance
        match["candidates_passing_both"] = candidates_passing_both
        match["candidates_rejected_no_dynamic"] = candidates_rejected_no_dynamic
        match["total_triplets"] = total_triplets

    return results


def run_enzyme_design(
    enzyme_key: str,
    enzyme_config: Dict,
    rhodopsin_df: pd.DataFrame,
    chain_df: pd.DataFrame,
    rhodopsin_seq: str,
    icl_residues: List[Tuple],
    coord_750: np.ndarray,
    first_res: int,
    struct_proc: StructureProcessor,
    loader,
    data_root: Path,
) -> Optional[Dict]:
    """Run the Rhodozyme design workflow for a single reference enzyme.

    Returns a dict with all designs for this enzyme, or None on failure.
    """
    enzyme_id = enzyme_config["id"]
    enzyme_pdb = enzyme_config["pdb"]
    enzyme_chain = enzyme_config["chain"]
    triad_config = enzyme_config["triad"]
    substrate_smiles = enzyme_config["substrate_smiles"]
    substrate_name = enzyme_config["substrate_name"]

    print(f"\n{'='*70}")
    print(f"ENZYME: {enzyme_config['description']}")
    print(f"PDB: {enzyme_pdb}, Chain: {enzyme_chain}")
    print(f"Substrate: {substrate_name} (SMILES: {substrate_smiles})")
    print(f"{'='*70}")

    # -------------------------------------------------------------------------
    # Step A: Load enzyme structure and extract triad geometry
    # -------------------------------------------------------------------------
    print(f"\n  [A] Loading enzyme structure {enzyme_pdb}...")
    enzyme_name = f"enzyme_{enzyme_key}"
    loader.download_and_register(enzyme_pdb, name=enzyme_name)
    enzyme_df = struct_proc.load_entity(enzyme_name)

    if enzyme_df is None:
        print(f"  ERROR: Failed to load {enzyme_pdb}")
        return None

    print(f"  Loaded {len(enzyme_df)} atoms")

    # Extract triad coordinates and CB vectors
    print(f"\n  [B] Extracting catalytic triad...")
    triad_coords = {}
    triad_cb_vecs = {}
    role_names = ["res1", "res2", "res3"]  # nucleophile, base, electrostatic

    for role in role_names:
        info = triad_config[role]
        # Allow per-residue chain specification (for multi-chain enzymes like chymotrypsin)
        res_chain = info.get("chain", enzyme_chain)
        coord = get_ca_coordinates(enzyme_df, info["res_num"], res_chain)
        cb_vec = get_sidechain_vector(enzyme_df, info["res_num"], res_chain)
        if coord is None:
            print(f"    ERROR: Could not find {info['res_name']}{info['res_num']} in chain {res_chain}")
            return None
        triad_coords[role] = coord
        triad_cb_vecs[role] = cb_vec
        cb_str = f"CB vec: [{cb_vec[0]:.2f}, {cb_vec[1]:.2f}, {cb_vec[2]:.2f}]" if cb_vec is not None else "no CB"
        chain_info = f" (chain {res_chain})" if res_chain != enzyme_chain else ""
        print(f"    {info['res_name']}{info['res_num']}{chain_info} ({info['role']}): CA={coord.round(2)}, {cb_str}")

    # Compute target triangle distances (CA-CA)
    d_r1_r2 = np.linalg.norm(triad_coords["res1"] - triad_coords["res2"])
    d_r2_r3 = np.linalg.norm(triad_coords["res2"] - triad_coords["res3"])
    d_r1_r3 = np.linalg.norm(triad_coords["res1"] - triad_coords["res3"])
    target_distances = (d_r1_r2, d_r2_r3, d_r1_r3)

    print(f"\n  Triad geometry (CA-CA distances):")
    print(f"    {triad_config['res1']['res_name']}-{triad_config['res2']['res_name']}: {d_r1_r2:.2f} A")
    print(f"    {triad_config['res2']['res_name']}-{triad_config['res3']['res_name']}: {d_r2_r3:.2f} A")
    print(f"    {triad_config['res1']['res_name']}-{triad_config['res3']['res_name']}: {d_r1_r3:.2f} A")

    # -------------------------------------------------------------------------
    # Step C: Find top N matching triplets
    # -------------------------------------------------------------------------
    print(f"\n  [C] Searching for top {TOP_N_MATCHES} matching triplet geometries...")
    print(f"    Distance tolerance: {DISTANCE_TOLERANCE}A, CB direction tolerance: 30°")
    if REQUIRE_DYNAMIC_HELIX:
        print(f"    Dynamic helix requirement: ENABLED (at least one residue on TM5/TM6)")

    source_ca = [triad_coords["res1"], triad_coords["res2"], triad_coords["res3"]]
    source_cb_vecs = [triad_cb_vecs["res1"], triad_cb_vecs["res2"], triad_cb_vecs["res3"]]

    top_matches = find_top_triplets(
        icl_residues,
        source_ca=source_ca,
        source_cb_vecs=source_cb_vecs,
        target_distances=target_distances,
        distance_tolerance=DISTANCE_TOLERANCE,
        angle_tolerance=30.0,
        top_n=TOP_N_MATCHES,
        require_dynamic_helix=REQUIRE_DYNAMIC_HELIX,
    )

    if not top_matches:
        print("    WARNING: No triplets found within tolerance, relaxing...")
        top_matches = find_top_triplets(
            icl_residues,
            source_ca=source_ca,
            source_cb_vecs=source_cb_vecs,
            target_distances=target_distances,
            distance_tolerance=10.0,
            angle_tolerance=90.0,
            top_n=TOP_N_MATCHES,
            require_dynamic_helix=REQUIRE_DYNAMIC_HELIX,
        )

    if not top_matches:
        print("    ERROR: Could not find any matching triplets")
        return None

    print(f"\n    Found {len(top_matches)} matches passing filters")
    if top_matches[0].get('candidates_rejected_no_dynamic', 0) > 0:
        print(f"    Triplets rejected (no TM5/TM6): {top_matches[0]['candidates_rejected_no_dynamic']}")
    print(f"    Candidates passing distance filter: {top_matches[0]['candidates_passing_distance']}")
    print(f"    Candidates passing both filters: {top_matches[0]['candidates_passing_both']}")

    print(f"\n    TOP {len(top_matches)} MATCHES:")
    print("    " + "-" * 80)
    for i, match in enumerate(top_matches):
        dynamic_marker = " [TM5/6]" if match.get('has_dynamic_helix', False) else ""
        grn_info = f"({match.get('ser_grn', '?')}, {match.get('his_grn', '?')}, {match.get('asp_grn', '?')})"
        print(f"    #{i+1}: {match['ser_name']}{match['ser_pos']}->{match['his_name']}{match['his_pos']}->{match['asp_name']}{match['asp_pos']}"
              f"  | score={match['total_score']:.3f} | dist={match['distance_rmsd']:.2f}A | cb={match['cb_direction_rmsd']:.1f}°{dynamic_marker}")
        print(f"         GRN: {grn_info}")
    print("    " + "-" * 80)

    # -------------------------------------------------------------------------
    # Step D: Generate chimeric sequences
    # -------------------------------------------------------------------------
    print(f"\n  [D] Generating {len(top_matches)} chimeric sequences...")

    designs = []
    seq_proc = SequenceProcessor()

    for i, match in enumerate(top_matches):
        design_id = f"Rhodozyme-{enzyme_id}_{i+1:02d}"

        # Create mutations (using generic names ser/his/asp for the 3 sites)
        # Use AA3TO1 to convert 3-letter codes (ASP, SER, etc.) to 1-letter (D, S, etc.)
        mutations = [
            (match["ser_pos"], match["ser_name"], AA3TO1[triad_config["res1"]["res_name"]]),
            (match["his_pos"], match["his_name"], AA3TO1[triad_config["res2"]["res_name"]]),
            (match["asp_pos"], match["asp_name"], AA3TO1[triad_config["res3"]["res_name"]]),
        ]

        # Apply mutations to create chimera
        seq_list = list(rhodopsin_seq)
        mutation_labels = []
        for res_num, old_name, new_aa in mutations:
            seq_idx = res_num - first_res
            if 0 <= seq_idx < len(seq_list):
                old_aa_1letter = seq_list[seq_idx]
                seq_list[seq_idx] = new_aa
                mutation_labels.append(f"{old_aa_1letter}{res_num}{new_aa}")

        chimera_seq = "".join(seq_list)
        seq_proc.save_entity(design_id, chimera_seq)

        design = {
            "id": design_id,
            "sequence": chimera_seq,
            "mutations": mutation_labels,
            "match": match,
            "enzyme": enzyme_key,
            "enzyme_id": enzyme_id,
            "substrate_smiles": substrate_smiles,
            "substrate_name": substrate_name,
        }
        designs.append(design)
        print(f"    {design_id}: {', '.join(mutation_labels)}")

    # -------------------------------------------------------------------------
    # Step E: Prepare Boltz prediction jobs with substrate
    # -------------------------------------------------------------------------
    print(f"\n  [E] Preparing Boltz structure prediction jobs with substrate...")

    # Create output directory for this enzyme
    enzyme_output_dir = OUTPUT_DIR / f"rhodozyme_{enzyme_id}"
    enzyme_output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare sequences dict
    sequences = {d["id"]: d["sequence"] for d in designs}

    # Save to FASTA
    fasta_path = enzyme_output_dir / "sequences.fasta"
    with open(fasta_path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")
    print(f"    Saved: {fasta_path}")

    # Prepare Boltz jobs
    try:
        from protos.models.model_manager import ModelManager

        manager = ModelManager(data_root=data_root)

        # Register dataset for this enzyme's designs
        dataset_name = f"rhodozyme_{enzyme_id}_designs"
        seq_proc.save_sequences(
            sequences,
            output_file=dataset_name,
            dataset_name=dataset_name,
            materialize_entities=True,
        )

        print(f"\n    Preparing Boltz2 jobs for {len(sequences)} sequences...")
        boltz_jobs = []

        for entity_name in sequences.keys():
            invocation = manager.prepare(
                "boltz2",
                inputs={
                    "sequence_dataset": dataset_name,
                    "entity": entity_name,
                },
                config={
                    # Include both retinal (to block its binding pocket) and substrate
                    "ligands": [
                        {
                            "id": RETINAL_ID,
                            "smiles": RETINAL_SMILES,
                            "is_cofactor": True,  # Marks as cofactor, not affinity binder
                        },
                        {
                            "id": "SUB",
                            "smiles": substrate_smiles,
                            "name": substrate_name,
                        },
                    ],
                },
            )

            if invocation.job is not None:
                job = invocation.job
                boltz_jobs.append((entity_name, job))
                print(f"\n    Job for {entity_name}:")
                print(f"      Working dir: {job.working_dir}")

        if boltz_jobs:
            print(f"\n    {len(boltz_jobs)} Boltz jobs prepared")

    except Exception as e:
        print(f"    Boltz job preparation failed: {e}")
        import traceback
        traceback.print_exc()

    # Save enzyme design summary
    best_match = top_matches[0]
    enzyme_summary = {
        "enzyme": {
            "key": enzyme_key,
            "id": enzyme_id,
            "pdb": enzyme_pdb,
            "description": enzyme_config["description"],
            "triad": {
                triad_config["res1"]["role"]: f"{triad_config['res1']['res_name']}{triad_config['res1']['res_num']}",
                triad_config["res2"]["role"]: f"{triad_config['res2']['res_name']}{triad_config['res2']['res_num']}",
                triad_config["res3"]["role"]: f"{triad_config['res3']['res_name']}{triad_config['res3']['res_num']}",
            },
            "substrate": {
                "name": substrate_name,
                "smiles": substrate_smiles,
            },
        },
        "target_geometry": {
            "d_r1_r2": float(d_r1_r2),
            "d_r2_r3": float(d_r2_r3),
            "d_r1_r3": float(d_r1_r3),
        },
        "num_designs": len(designs),
        "best_design": {
            "id": designs[0]["id"],
            "mutations": designs[0]["mutations"],
            "distance_rmsd": float(best_match["distance_rmsd"]),
            "cb_direction_rmsd": float(best_match["cb_direction_rmsd"]),
        },
        "all_designs": [
            {
                "id": d["id"],
                "mutations": d["mutations"],
                "score": float(d["match"]["total_score"]),
            }
            for d in designs
        ],
    }

    summary_path = enzyme_output_dir / "design_summary.json"
    with open(summary_path, "w") as f:
        json.dump(enzyme_summary, f, indent=2)
    print(f"    Saved: {summary_path}")

    return {
        "enzyme_key": enzyme_key,
        "enzyme_id": enzyme_id,
        "designs": designs,
        "top_matches": top_matches,
        "summary": enzyme_summary,
    }


def main() -> int:
    """Run the Rhodozyme design workflow for all reference enzymes."""
    print("=" * 70)
    print("RHODOZYME DESIGN WORKFLOW")
    print("Light-Activated Enzyme Engineering")
    print("=" * 70)
    print()
    print("CONCEPT:")
    print("  Graft catalytic triads from multiple enzymes onto rhodopsin's")
    print("  intracellular domain, which undergoes large conformational changes")
    print("  upon light activation. Goal: light-controlled enzyme activity.")
    print()
    print(f"  Scaffold: {RHODOPSIN_PDB} (active rhodopsin, open ICL)")
    print(f"  Reference enzymes: {', '.join(REFERENCE_ENZYMES.keys())}")
    print("=" * 70)

    # Initialize ProtOS
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    # -------------------------------------------------------------------------
    # Step 1: Load rhodopsin scaffold
    # -------------------------------------------------------------------------
    print("\n[1] Loading rhodopsin scaffold...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    loader.download_and_register(RHODOPSIN_PDB, name="rhodopsin_active")
    rhodopsin_df = struct_proc.load_entity("rhodopsin_active")

    if rhodopsin_df is None:
        print("  ERROR: Failed to load rhodopsin")
        return 1

    print(f"  Rhodopsin ({RHODOPSIN_PDB}): {len(rhodopsin_df)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Identify rhodopsin intracellular domain
    # -------------------------------------------------------------------------
    print("\n[2] Identifying rhodopsin intracellular domain...")

    # Annotate with GRN to find position 7.50
    grn_df = struct_proc.annotate_with_grn("rhodopsin_active", chains=[RHODOPSIN_CHAIN])

    pos_750 = None
    if grn_df is not None and not grn_df.empty:
        grn_unique = grn_df[["auth_seq_id", "grn"]].drop_duplicates()
        grn_750 = grn_unique[grn_unique["grn"] == "7.50"]
        if len(grn_750) == 0:
            grn_750 = grn_unique[grn_unique["grn"].str.contains("7.50|7x50", na=False, regex=True)]
        if len(grn_750) > 0:
            pos_750 = int(grn_750["auth_seq_id"].values[0])
            print(f"  Position 7.50 (NPxxY): residue {pos_750}")

    if pos_750 is None:
        pos_750 = 306  # Fallback for rhodopsin
        print(f"  Position 7.50 (estimated): residue {pos_750}")

    # Get chain data
    chain_df = rhodopsin_df[rhodopsin_df["auth_chain_id"] == RHODOPSIN_CHAIN]
    coord_750 = get_ca_coordinates(rhodopsin_df, pos_750, RHODOPSIN_CHAIN)
    if coord_750 is None:
        print(f"  ERROR: Could not find residue {pos_750}")
        return 1

    print(f"  7.50 coordinates: {coord_750}")

    # -------------------------------------------------------------------------
    # Step 3: Filter to intracellular region using GRN annotation
    # -------------------------------------------------------------------------
    print("\n[3] Filtering to intracellular region using GRN...")

    # Use GRN annotation to identify cytoplasmic-facing residues
    # This is more accurate than distance-based filtering
    if grn_df is not None and not grn_df.empty:
        # Get unique residue-GRN mapping
        grn_unique = grn_df[["auth_seq_id", "res_name", "grn"]].drop_duplicates()

        # Filter to cytoplasmic residues based on GRN
        grn_unique["is_cyto"] = grn_unique["grn"].apply(is_cytoplasmic_residue)
        cyto_residues = grn_unique[grn_unique["is_cyto"]]

        print(f"  GRN-annotated cytoplasmic residues: {len(cyto_residues)}")

        # Count dynamic helix residues
        cyto_residues["is_dynamic"] = cyto_residues["grn"].apply(is_on_dynamic_helix)
        n_dynamic = cyto_residues["is_dynamic"].sum()
        print(f"  Dynamic helix residues (TM5/TM6): {n_dynamic}")

        # Build ICL residue list with coordinates, CB vectors, and GRN
        icl_residues = []
        for _, row in cyto_residues.iterrows():
            res_num = int(row["auth_seq_id"])
            res_name = row["res_name"]
            grn = row["grn"]
            coord = get_ca_coordinates(chain_df, res_num)
            if coord is not None:
                cb_vec = get_sidechain_vector(chain_df, res_num)
                if cb_vec is not None:
                    icl_residues.append((res_num, res_name, coord, cb_vec, grn))

        print(f"  Cytoplasmic residues with valid coordinates: {len(icl_residues)}")

    else:
        # Fallback to distance-based filtering if GRN not available
        print("  WARNING: GRN annotation not available, using distance-based fallback")
        print("  WARNING: Dynamic helix filtering disabled (requires GRN)")
        ca_atoms = chain_df[chain_df["atom_name"] == "CA"].copy()
        ca_coords = ca_atoms[["x", "y", "z"]].values
        distances = calculate_distance(ca_coords, coord_750)
        icl_mask = distances <= ICL_DISTANCE_FROM_750
        icl_ca_atoms = ca_atoms[icl_mask]

        icl_residues = []
        for _, row in icl_ca_atoms.iterrows():
            res_num = int(row["auth_seq_id"])
            res_name = row["res_name"]
            coord = row[["x", "y", "z"]].values
            cb_vec = get_sidechain_vector(chain_df, res_num)
            if cb_vec is not None:
                # No GRN available - use empty string
                icl_residues.append((res_num, res_name, coord, cb_vec, ""))

        print(f"  Residues within {ICL_DISTANCE_FROM_750}A of 7.50: {len(icl_residues)}")

    # Show summary with GRN info
    icl_info = [(r[0], r[1], r[4]) for r in icl_residues]  # (res_num, res_name, grn)
    print(f"  ICL residues: {icl_info[:10]}..." if len(icl_info) > 10 else f"  ICL residues: {icl_info}")

    # Get rhodopsin sequence and first residue number
    rhodopsin_seq = struct_proc.get_sequence("rhodopsin_active", chain_id=RHODOPSIN_CHAIN)
    if not rhodopsin_seq:
        print("  ERROR: Could not extract sequence")
        return 1

    first_res = chain_df["auth_seq_id"].min()
    print(f"  Rhodopsin sequence: {len(rhodopsin_seq)} residues (first res: {first_res})")

    # -------------------------------------------------------------------------
    # Step 4: Run design workflow for each reference enzyme
    # -------------------------------------------------------------------------
    print("\n[4] Running design workflow for each reference enzyme...")

    all_results = {}
    for enzyme_key, enzyme_config in REFERENCE_ENZYMES.items():
        result = run_enzyme_design(
            enzyme_key=enzyme_key,
            enzyme_config=enzyme_config,
            rhodopsin_df=rhodopsin_df,
            chain_df=chain_df,
            rhodopsin_seq=rhodopsin_seq,
            icl_residues=icl_residues,
            coord_750=coord_750,
            first_res=first_res,
            struct_proc=struct_proc,
            loader=loader,
            data_root=data_root,
        )
        if result is not None:
            all_results[enzyme_key] = result

    # -------------------------------------------------------------------------
    # Step 5: Create combined visualization
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("[5] Creating combined visualization...")
    import plotly.graph_objects as go

    # Get ICL coordinates for background
    icl_coords = np.array([r[2] for r in icl_residues])
    icl_labels = [f"{r[1]}{r[0]}" for r in icl_residues]

    fig = go.Figure()

    # Draw ICL region residues (background)
    fig.add_trace(go.Scatter3d(
        x=icl_coords[:, 0],
        y=icl_coords[:, 1],
        z=icl_coords[:, 2],
        mode="markers",
        marker=dict(size=5, color="#cccccc", opacity=0.4),
        text=icl_labels,
        name="ICL region",
        hovertemplate="%{text}<extra></extra>",
    ))

    # Add best design from each enzyme
    colors = {"TRP": "#e74c3c", "CHY": "#3498db", "PAP": "#2ecc71", "SUB": "#f39c12"}

    for enzyme_key, result in all_results.items():
        enzyme_id = result["enzyme_id"]
        best_match = result["top_matches"][0]
        best_design = result["designs"][0]

        # Get triad positions
        ser_coord = get_ca_coordinates(chain_df, best_match["ser_pos"])
        his_coord = get_ca_coordinates(chain_df, best_match["his_pos"])
        asp_coord = get_ca_coordinates(chain_df, best_match["asp_pos"])
        triad_positions = np.array([ser_coord, his_coord, asp_coord])

        color = colors.get(enzyme_id, "#9b59b6")

        # Add triad markers
        fig.add_trace(go.Scatter3d(
            x=triad_positions[:, 0],
            y=triad_positions[:, 1],
            z=triad_positions[:, 2],
            mode="markers",
            marker=dict(size=10, color=color),
            name=f"Rhodozyme-{enzyme_id}",
            hovertext=[
                f"{enzyme_id}: {best_match['ser_name']}{best_match['ser_pos']}",
                f"{enzyme_id}: {best_match['his_name']}{best_match['his_pos']}",
                f"{enzyme_id}: {best_match['asp_name']}{best_match['asp_pos']}",
            ],
            hovertemplate="%{hovertext}<extra></extra>",
        ))

        # Draw triangle edges
        for i, j in [(0, 1), (1, 2), (0, 2)]:
            fig.add_trace(go.Scatter3d(
                x=[triad_positions[i, 0], triad_positions[j, 0]],
                y=[triad_positions[i, 1], triad_positions[j, 1]],
                z=[triad_positions[i, 2], triad_positions[j, 2]],
                mode="lines",
                line=dict(color=color, width=3),
                showlegend=False,
                hoverinfo="skip",
            ))

    # Mark position 7.50
    fig.add_trace(go.Scatter3d(
        x=[coord_750[0]],
        y=[coord_750[1]],
        z=[coord_750[2]],
        mode="markers+text",
        marker=dict(size=10, color="#9b59b6", symbol="diamond"),
        text=["7.50"],
        textposition="bottom center",
        name="Position 7.50 (NPxxY)",
    ))

    fig.update_layout(
        title=dict(
            text="Rhodozyme Designs: Catalytic Triads Grafted onto Rhodopsin ICL<br>"
                 f"<sub>Reference enzymes: {', '.join(all_results.keys())} | "
                 f"Top {TOP_N_MATCHES} designs per enzyme</sub>",
            x=0.5,
        ),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
        ),
        height=700,
        width=900,
    )

    fig.write_html(str(FIGURES_DIR / "rhodozyme_designs_combined.html"))
    print(f"  Saved: {FIGURES_DIR / 'rhodozyme_designs_combined.html'}")

    # -------------------------------------------------------------------------
    # Step 6: Save global summary
    # -------------------------------------------------------------------------
    print("\n[6] Saving global design summary...")

    global_summary = {
        "scaffold": {
            "pdb": RHODOPSIN_PDB,
            "chain": RHODOPSIN_CHAIN,
            "description": "Active rhodopsin (metarhodopsin II)",
            "position_750": pos_750,
            "icl_residue_count": len(icl_residues),
        },
        "reference_enzymes": list(REFERENCE_ENZYMES.keys()),
        "designs_per_enzyme": TOP_N_MATCHES,
        "total_designs": sum(len(r["designs"]) for r in all_results.values()),
        "enzyme_results": {
            enzyme_key: {
                "id": r["enzyme_id"],
                "num_designs": len(r["designs"]),
                "best_design": {
                    "id": r["designs"][0]["id"],
                    "mutations": r["designs"][0]["mutations"],
                    "score": float(r["top_matches"][0]["total_score"]),
                },
            }
            for enzyme_key, r in all_results.items()
        },
    }

    with open(OUTPUT_DIR / "rhodozyme_global_summary.json", "w") as f:
        json.dump(global_summary, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'rhodozyme_global_summary.json'}")

    # -------------------------------------------------------------------------
    # Final Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RHODOZYME DESIGN WORKFLOW COMPLETE")
    print("=" * 70)
    print(f"Scaffold: {RHODOPSIN_PDB} (active rhodopsin)")
    print(f"Enzymes processed: {len(all_results)}/{len(REFERENCE_ENZYMES)}")
    print()
    print("Best designs per enzyme:")
    for enzyme_key, r in all_results.items():
        best = r["designs"][0]
        score = r["top_matches"][0]["total_score"]
        print(f"  {best['id']}: {', '.join(best['mutations'])} (score: {score:.3f})")
    print()
    print(f"Total designs generated: {global_summary['total_designs']}")
    print()
    print("NEXT STEPS:")
    print("  1. Run Boltz2 to predict chimera structures")
    print("  2. Verify catalytic geometry is maintained")
    print("  3. Check solvent accessibility of active sites")
    print("  4. Compare designs across enzyme types")
    print()
    print(f"Outputs: {OUTPUT_DIR}")
    print(f"  Per-enzyme: {OUTPUT_DIR}/rhodozyme_<ID>/")
    print(f"Figures: {FIGURES_DIR}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

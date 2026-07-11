"""
Rhodozyme figure configuration.

Change these variables to switch which design is shown in all figures.
All figure scripts import from here.

Usage:
    from rhodozyme_config import CFG
    print(CFG["design_pdb"])
"""

import os
from pathlib import Path

# ===========================================================================
# CHANGE THESE TO SWITCH CANDIDATE
# ===========================================================================

PLACEMENT_NUM = "00"   # Which placement (00, 02)
DESIGN_NUM = "6"       # Which RFdiffusion design (0-49; best pLDDT + good RMSD)
SEQ_NUM = "6"          # Which LigandMPNN sequence (0-7)

# Theozyme residue positions (rhodopsin numbering, same across all stages)
THEOZYME_SER = 230
THEOZYME_HIS = 245
THEOZYME_ASP = 250

# ===========================================================================
# DERIVED PATHS (do not edit)
# ===========================================================================

REPO_ROOT = Path(__file__).resolve().parents[5]
FIGURE_ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("PROTOS_MODEL_DATA", REPO_ROOT / "data" / "models"))

THEOZYME_RESI = [THEOZYME_HIS, THEOZYME_ASP, THEOZYME_SER]

RUN_DIR = DATA / "rfdiffusion2" / f"20260203_placement_{PLACEMENT_NUM}_production"
BOLTZ_DIR = DATA / "boltz2" / f"20260203_placement_{PLACEMENT_NUM}_production"

CFG = {
    "placement_num": PLACEMENT_NUM,
    "design_num": DESIGN_NUM,
    "seq_num": SEQ_NUM,
    "theozyme_his": THEOZYME_HIS,
    "theozyme_asp": THEOZYME_ASP,
    "theozyme_ser": THEOZYME_SER,
    "theozyme_resi": THEOZYME_RESI,
    # Input structures
    "active_pdb": DATA / "rfdiffusion2" / "input" / "3pqr_chainA_ori.pdb",
    "placement_pdb": DATA / "rfdiffusion2" / "input" / f"placement_{PLACEMENT_NUM}_triad_ori.pdb",
    # RFdiffusion output
    "design_pdb": RUN_DIR / f"rhodozyme_{PLACEMENT_NUM}_{DESIGN_NUM}-atomized-bb-False.pdb",
    # LigandMPNN output
    "design_fa": RUN_DIR / "ligmpnn" / "seqs" / f"rhodozyme_{PLACEMENT_NUM}_{DESIGN_NUM}-atomized-bb-False.fa",
    # Boltz2 output (thesis copy with Ser-230 OG rotated 180 deg around CA-CB)
    "boltz_cif": FIGURE_ROOT / "data" / "boltz2_design6_seq6.cif",
    # Directories
    "run_dir": RUN_DIR,
    "boltz_dir": BOLTZ_DIR,
}

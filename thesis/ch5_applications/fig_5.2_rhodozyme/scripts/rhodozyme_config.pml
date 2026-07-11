# Rhodozyme Figure Configuration
# ===========================================================================
# Change these variables to switch which design is shown in all figures.
# Then re-run the individual figure scripts.
#
# To change candidate: update DESIGN_NUM, SEQ_NUM, and THEOZYME_RESI.
# All figure scripts source this file via: @rhodozyme_config.pml
# ===========================================================================

# All configuration is in the python block below.
# To change candidate: edit PLACEMENT_NUM, DESIGN_NUM, SEQ_NUM, and
# the THEOZYME residue numbers.

python

from pymol import cmd
import os

# Base data directory. Run PyMOL from the repository root, or override this
# with PROTOS_MODEL_DATA when model outputs live elsewhere.
DATA = os.environ.get("PROTOS_MODEL_DATA", os.path.abspath("data/models"))

# Placement
PLACEMENT_NUM = "00"
DESIGN_NUM = "6"
SEQ_NUM = "6"

# Theozyme residues (rhodopsin numbering)
THEOZYME_SER = 230
THEOZYME_HIS = 245
THEOZYME_ASP = 250
THEOZYME_RESI = f"{THEOZYME_HIS}+{THEOZYME_ASP}+{THEOZYME_SER}"

# Derived paths
PLACEMENT_PDB = f"{DATA}/rfdiffusion2/input/placement_{PLACEMENT_NUM}_triad_ori.pdb"
ACTIVE_PDB = f"{DATA}/rfdiffusion2/input/3pqr_chainA_ori.pdb"

RUN_DIR = f"{DATA}/rfdiffusion2/20260203_placement_{PLACEMENT_NUM}_production"
DESIGN_PDB = f"{RUN_DIR}/rhodozyme_{PLACEMENT_NUM}_{DESIGN_NUM}-atomized-bb-False.pdb"
DESIGN_FA = f"{RUN_DIR}/ligmpnn/seqs/rhodozyme_{PLACEMENT_NUM}_{DESIGN_NUM}-atomized-bb-False.fa"

BOLTZ_DIR = f"{DATA}/boltz2/20260203_placement_{PLACEMENT_NUM}_production"
BOLTZ_CIF_ORIG = f"{BOLTZ_DIR}/rhodozyme_{PLACEMENT_NUM}_{DESIGN_NUM}-atomized-bb-False_seq{SEQ_NUM}/outputs/boltz_results_config/predictions/config/config_model_0.cif"
# Thesis copy with Ser-230 OG rotated 180 deg around CA-CB
BOLTZ_CIF = os.path.abspath(
    "thesis/ch5_applications/fig_5.2_rhodozyme/data/boltz2_design6_seq6.cif"
)

# Export to global scope for other python blocks
import builtins
builtins._rz_config = {
    "placement_num": PLACEMENT_NUM,
    "design_num": DESIGN_NUM,
    "seq_num": SEQ_NUM,
    "theozyme_his": THEOZYME_HIS,
    "theozyme_asp": THEOZYME_ASP,
    "theozyme_ser": THEOZYME_SER,
    "theozyme_resi": THEOZYME_RESI,
    "placement_pdb": PLACEMENT_PDB,
    "active_pdb": ACTIVE_PDB,
    "design_pdb": DESIGN_PDB,
    "design_fa": DESIGN_FA,
    "boltz_cif": BOLTZ_CIF,
    "run_dir": RUN_DIR,
    "boltz_dir": BOLTZ_DIR,
}

# Verify files exist
missing = []
for key in ["placement_pdb", "design_pdb", "boltz_cif"]:
    path = builtins._rz_config[key]
    if not os.path.exists(path):
        missing.append(f"  {key}: {path}")

if missing:
    print("WARNING - missing files:")
    for m in missing:
        print(m)
else:
    print("All design files found.")

print(f"Config: placement={PLACEMENT_NUM}, design={DESIGN_NUM}, seq={SEQ_NUM}")
print(f"Theozyme: HIS-{THEOZYME_HIS}, ASP-{THEOZYME_ASP}, SER-{THEOZYME_SER}")

python end

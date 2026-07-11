#!/usr/bin/env python3
"""Minimal GRNProcessor test — matches draft code example."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
protos.set_data_path(str(REPO_ROOT / "data"))

from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor

seq_proc = SequenceProcessor()

# Annotate a dataset with GRN positions (Ballesteros-Weinstein)
grn_table, summary = seq_proc.annotate_with_grn(
    dataset_name="opsin_atlas",
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="atlas_grn",
    return_summary=True,
)
print(f"Annotated: {summary['global']['annotated']}/{summary['global']['total']}")

# Load the GRN table back
grn_proc = GRNProcessor()
table = grn_proc.load_table("atlas_grn")
print(f"GRN table: {table.shape}")

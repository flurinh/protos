#!/usr/bin/env python3
"""Minimal StructureProcessor test — matches draft code example."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
protos.set_data_path(str(REPO_ROOT / "data"))

from protos.processing.structure import StructureProcessor
from protos.io.ingest.structure_loader import StructureLoader

struct_proc = StructureProcessor()
loader = StructureLoader(processor=struct_proc)

# Download structures from PDB
loader.download_batch(["1C3W", "1U19"], dataset_name="opsin_structures")

# Extract retinal binding pocket contacts
contacts = struct_proc.get_ligand_interactions(
    "1C3W", ligand_id="RET", chain_id="A", cutoff=4.0
)
print(f"1C3W binding residues: {len(contacts.get('binding_residues', []))}")

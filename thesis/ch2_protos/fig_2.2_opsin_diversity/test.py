#!/usr/bin/env python3
"""Minimal SequenceProcessor test — matches draft code example."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
protos.set_data_path(str(REPO_ROOT / "data"))

from protos.processing.sequence import SequenceProcessor
from protos.io.ingest.sequence_loader import SequenceLoader
from protos.io.ingest.ncbi_loader import NCBILoader

seq_proc = SequenceProcessor()

# Download and register a sequence from UniProt
loader = SequenceLoader(processor=seq_proc)
loader.download_and_register("uniprot:P02699", name="RHO_BOVIN")

seq = seq_proc.load_entity("RHO_BOVIN")
print(f"RHO_BOVIN: {len(seq)} aa")

# BLAST search for homologs
ncbi = NCBILoader(processor=seq_proc)
hits = ncbi.blast_search(
    sequence=seq,
    database="swissprot",
    hitlist_size=10,
)
print(f"BLAST hits: {hits}")

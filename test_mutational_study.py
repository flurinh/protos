#!/usr/bin/env python3
"""Mutational study helper: generate mutant FASTA datasets from WT + specs.

Zero-config: the script sets ProtosPaths to protos/data and relies on
SequenceProcessor and GRNProcessor to manage paths and registration.

Usage: run as a demo to create a small GPCR mutant set from bovine rhodopsin.
"""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import sys

SRC_DIR = Path(__file__).resolve().parent / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.processing.sequence import SequenceProcessor


def ensure_data_root() -> Path:
    # Always prefer protos/data per repository guidelines
    data_root = Path(__file__).resolve().parent / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    try:
        protos.set_data_path(str(data_root))
    except Exception:
        pass
    return data_root


from protos.processing.grn import GRNProcessor


def demo_mutational_study() -> None:
    ensure_data_root()
    sp = SequenceProcessor()

    # Base WT sequence (bovine rhodopsin)
    seq_id = "OPSD_BOVIN"
    wt = (
        "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGF"
        "TTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERMWVWRTYFTYDENNCECILFMTSALGVPINFLIYFNFSMVFLVCGPNK"
        "NPFEGERQCVMMKEKKLRHVLFGGNKYNKQFRYRMERTACQEEGTHSTSTLQYERLTPQEWYTWQFSSFITIPHSNQSNFDFIVATARHALA"
        "FTSNHICIYPWFGHFTVYVVGALGVGVLFGPVFFDYNKSEHGFPPTNMTIPAFFAKSSAIYNPVIYIMMNKQFRNCMVTTICCGKNPLGDTD"
        "DEAERNETEVDKTVINAFFTLYVTVLAFVFTVPLFVTAIYNVEMKQFRNVEDRMVTHQLVSSMVLNAEDGGNDPVLTPPEQKTSAVLAILY"
        "GNTLVLFMENYIIYVVKEAAAQQQESATTQTAENPSQSEAAPQEPSSQFRASA"
    )

    ds_wt = "rhodopsin_wt"
    sp.save_sequences({seq_id: wt}, output_file=ds_wt, dataset_name=ds_wt, materialize_entities=True)

    # Example mutational spec:
    # - Sequence positions (1-based): 3 -> ['L', 'A']; 10 -> ['R']
    # - GRN labels (for GPCR): '3.50' -> ['A']
    seq_pos = {3: ['L', 'A'], 10: ['R']}
    grn_pos = {'3.50': ['A']}

    # Annotate GRN (GPCR A family)
    grn_table = sp.annotate_with_grn(dataset_name=ds_wt, reference_table="gpcrdb_ref", protein_family="gpcr_a", output_table=f"{ds_wt}_grn", allow_create=True, return_summary=False)

    mutants = sp.generate_mutants_for_sequence(
        seq_id,
        wt,
        seq_positions=seq_pos,
        grn_positions=grn_pos,
        grn_table=grn_table,
        protein_family="gpcr_a",
    )

    ds_mut = f"{ds_wt}__mutational_study"
    payload = dict(mutants)
    payload[seq_id] = wt
    sp.save_sequences(payload, output_file=ds_mut, dataset_name=ds_mut, materialize_entities=True)

    print(f"Created dataset '{ds_mut}' with {len(mutants)} mutants + WT")


def demo_grn_mutational_study_multi() -> None:
    """Demonstrate GRN-level mutational screen across multiple sequences.

    Logic:
    sequences -> GRN annotation -> GRN→seq position map per sequence ->
    generate mutants per sequence using the same GRN positions.
    """
    ensure_data_root()
    sp = SequenceProcessor()

    # Two GPCR sequences (bovine rhodopsin + ADRB2)
    seqs = {
        "OPSD_BOVIN": (
            "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGF"
            "TTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERMWVWRTYFTYDENNCECILFMTSALGVPINFLIYFNFSMVFLVCGPNK"
            "NPFEGERQCVMMKEKKLRHVLFGGNKYNKQFRYRMERTACQEEGTHSTSTLQYERLTPQEWYTWQFSSFITIPHSNQSNFDFIVATARHALA"
            "FTSNHICIYPWFGHFTVYVVGALGVGVLFGPVFFDYNKSEHGFPPTNMTIPAFFAKSSAIYNPVIYIMMNKQFRNCMVTTICCGKNPLGDTD"
            "DEAERNETEVDKTVINAFFTLYVTVLAFVFTVPLFVTAIYNVEMKQFRNVEDRMVTHQLVSSMVLNAEDGGNDPVLTPPEQKTSAVLAILY"
            "GNTLVLFMENYIIYVVKEAAAQQQESATTQTAENPSQSEAAPQEPSSQFRASA"
        ),
        # A compact ADRB2 snippet (for demo); real workflow would use full-length
        "ADRB2_HUMAN": (
            "MDPALINCYENSSSNGNTGEQSGYHVEQEKENKLLCEDLQNMKTWFKVLNSNTVTFHVDTDNQTVLHNLDLQV"
            "HNFSDTFHTLENVQNLSQANETCIPLDGTRSTFGSLSAFSTASFQDI"
        ),
    }

    ds_name = "gpcr_demo_wt"
    sp.save_sequences(seqs, output_file=ds_name, dataset_name=ds_name, materialize_entities=True)

    # Shared GRN mutation spec: mutate the same GRN positions across both sequences
    grn_positions = {
        "3.50": ["A", "G"],
        "6.44": ["V"],
    }

    # Annotate GRN once for the dataset
    grn_table = sp.annotate_with_grn(
        dataset_name=ds_name,
        reference_table="gpcrdb_ref",
        protein_family="gpcr_a",
        output_table=f"{ds_name}_grn",
        allow_create=True,
        return_summary=False,
    )

    all_mutants: Dict[str, str] = {}
    for seq_id, wt in seqs.items():
        mutants = sp.generate_mutants_for_sequence(
            seq_id,
            wt,
            seq_positions=None,
            grn_positions=grn_positions,
            grn_table=grn_table,
            protein_family="gpcr_a",
        )
        all_mutants.update(mutants)

    ds_mut = f"{ds_name}__grn_mutational_study"
    sp.save_sequences(all_mutants, output_file=ds_mut, dataset_name=ds_mut, materialize_entities=True)
    print(
        f"Created GRN-mutant dataset '{ds_mut}' with {len(all_mutants)} mutants (across {len(seqs)} sequences)"
    )


if __name__ == "__main__":
    demo_mutational_study()
    demo_grn_mutational_study_multi()

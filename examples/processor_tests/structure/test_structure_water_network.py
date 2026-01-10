#!/usr/bin/env python3
"""Annotate a GPCR structure with GRNs and report water-mediated networks."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, List

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos

from protos.analysis.structure_water_networks import summarize_water_networks
from protos.processing.structure import StructureProcessor
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.grn import GRNProcessor


STRUCTURE_IDS = ["3sn6", "5d5a"]  # Representative beta2-adrenergic and A2A receptors
REFERENCE_TABLE = "gpcrdb_ref"
PROTEIN_FAMILY = "gpcr_a"
GRN_TABLE = "gpcr_structure_water_networks_grn"


def ensure_data_root() -> Path:
    """Create a data directory for the script and configure Protos paths."""

    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def ensure_structures(processor: StructureProcessor, loader: StructureLoader, structure_ids: Iterable[str]) -> List[str]:
    """Download structures if necessary and return the ones available locally."""

    available: List[str] = []
    for struct_id in structure_ids:
        frame = processor.load_entity(struct_id)
        if frame is not None:
            available.append(struct_id)
            continue

        print(f"Downloading {struct_id} ...", end=" ")
        success, failed = loader.download_batch([struct_id], dataset_name="water_network_structures", create_dataset=True)
        if success:
            available.extend(success)
            print("ok")
        else:
            print("failed")
            if failed:
                print(f"  ↳ Could not fetch {', '.join(failed)}")

    return available


def annotate_with_grn(processor: StructureProcessor, structure_ids: Iterable[str]) -> None:
    """Assign GRN labels to residues using the bundled GPCR reference table."""

    if not structure_ids:
        return

    print("Annotating structures with GRN numbering ...")
    annotations = processor.annotate_structures_with_grn(
        structure_ids,
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        output_table=GRN_TABLE,
        allow_create=True,
        materialize_entries=False,
        return_summary=False,
    )

    if isinstance(annotations, Dict):
        print(f"  Stored GRN annotations for {len(annotations)} sequences")


def apply_grn_to_structures(processor: StructureProcessor, structure_ids: Iterable[str]) -> None:
    """Project the recorded GRN table back onto structure residues."""

    grn_proc = GRNProcessor()
    table = grn_proc.load_table(GRN_TABLE)

    import re
    mapping_per_structure: Dict[str, Dict[tuple[str, int], str]] = {}

    structure_set = set(structure_ids)

    for entity_name, row in table.iterrows():
        if "_chain_" not in entity_name:
            continue
        structure_id, chain_token = entity_name.split("_chain_", 1)
        if structure_id not in structure_set:
            continue
        chain_id = chain_token.replace("__", "_")
        mapping = mapping_per_structure.setdefault(structure_id, {})
        for grn_label, residue_token in row.items():
            if not residue_token or residue_token == "-":
                continue
            match = re.search(r"(-?\d+)$", str(residue_token))
            if not match:
                continue
            try:
                auth_seq_id = int(match.group(1))
            except ValueError:
                continue
            mapping[(chain_id, auth_seq_id)] = str(grn_label)

    for structure_id, mapping in mapping_per_structure.items():
        if not mapping:
            continue
        processor.assign_grns(structure_id, mapping)



def analyze_water_networks(processor: StructureProcessor, structure_ids: Iterable[str]) -> Dict[str, Dict]:
    """Compute water networks for each structure and return the raw analysis."""

    print("Computing water-mediated residue networks ...")
    result = processor.compute_water_networks(
        structure_ids,
        residue_cutoff=3.4,
        water_water_cutoff=3.4,
        hydrogen_bond_cutoff=3.2,
    )

    errors = result.get("errors", {})
    for struct_id, message in errors.items():
        print(f"  ⚠ {struct_id}: {message}")

    return result.get("structures", {})
def report_networks(structure_id: str, analysis: Dict) -> None:
    """Print a GRN-centric report for all water networks in a structure."""

    summary = analysis.get("summary", {})
    print(f"\nStructure {structure_id}")
    print(f"  Networks: {summary.get('network_count', 0)} | Residues: {summary.get('residue_count', 0)} | Waters: {summary.get('water_count', 0)}")

    networks = analysis.get("networks", [])
    if not networks:
        print("  No water-mediated networks detected.")
        return

    for network in networks:
        chains = ", ".join(network.get("chains", [])) or "n/a"
        residue_display = network.get("residues", []) or ["n/a"]
        print(
            f"  Network {network.get('network_id') or '?'}: "
            f"residues={', '.join(residue_display)} | waters={network.get('waters', 0)} "
            f"(bridging={network.get('bridging_waters', 0)}) | max_chain={network.get('max_path_length', 0)} | chains={chains}"
        )

        condensed_paths = network.get("paths", [])
        if not condensed_paths:
            continue

        if network.get('max_path_length', 0) <= 2 or len(condensed_paths) <= 4:
            for path in condensed_paths:
                src = path.get("source")
                tgt = path.get("target")
                seq = path.get("sequence_str") or " → ".join(path.get("sequence", []))
                print(f"    {src} → {tgt} : {seq}")


def main() -> None:
    data_root = ensure_data_root()
    print(f"Protos data root: {data_root}")

    structure_processor = StructureProcessor()
    loader = StructureLoader(processor=structure_processor)

    available = ensure_structures(structure_processor, loader, STRUCTURE_IDS)
    if not available:
        print("No structures available; aborting analysis.")
        return

    annotate_with_grn(structure_processor, available)
    apply_grn_to_structures(structure_processor, available)

    analysis_map = analyze_water_networks(structure_processor, available)
    summary_map = summarize_water_networks(analysis_map)
    for struct_id in available:
        analysis = summary_map.get(struct_id)
        if analysis is None:
            print(f"\nStructure {struct_id}\n  No analysis results (see warnings above).")
            continue
        report_networks(struct_id, analysis)


if __name__ == "__main__":
    main()

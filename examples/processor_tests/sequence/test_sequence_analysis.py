#!/usr/bin/env python3
"""SequenceProcessor analytical capabilities demonstration.

This script demonstrates sequence analysis methods:
- Conservation analysis (per-position entropy and consensus)
- Linkage analysis (co-evolving positions via mutual information)
- Sequence metadata extraction (length, MW, pI, composition)
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def ensure_data_root() -> Path:
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


# Sample aligned GPCR sequences for analysis (same length required for conservation)
ALIGNED_SEQUENCES = {
    "seq_1": "MKTIIALSYIFCLVFADYKDDDDA",
    "seq_2": "MKTVIALSYIFCLVFADYKDDDDA",  # V at position 4
    "seq_3": "MKTIIALSYIFCLVFADYKDDDDA",
    "seq_4": "MKTIIALSYIFCMVFADYKDDDDA",  # M at position 13
    "seq_5": "MKTVIALSYIFCMVFADYKDDDDA",  # V at 4, M at 13 (co-occurring)
    "seq_6": "MKTIIALSYIFCLVFADYKDDDDA",
    "seq_7": "MKTVIALSYIFCMVFADYKDDDDA",  # V at 4, M at 13 (co-occurring)
    "seq_8": "MKTIIALSYIFCLVFADYKDDDDA",
}


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.sequence import SequenceProcessor

    print("=== SequenceProcessor Analytical Demo ===\n")
    processor = SequenceProcessor()

    # 1. Conservation Analysis
    print("1. Conservation Analysis")
    print("-" * 40)
    print("Computing per-position conservation metrics...")

    conservation = processor.compute_conservation(
        sequences=ALIGNED_SEQUENCES,
        normalize_entropy=True,
        pseudocount=0.0,
    )

    print(f"   Analyzed {len(ALIGNED_SEQUENCES)} sequences, {len(conservation)} positions\n")
    print("   First 10 positions:")
    cols = ["position", "consensus", "consensus_frequency", "normalized_entropy"]
    print(conservation.head(10)[cols].to_string(index=False))

    # Find variable positions (low consensus frequency or high entropy)
    variable = conservation[conservation["consensus_frequency"] < 1.0]
    if not variable.empty:
        print(f"\n   Variable positions (consensus_freq < 1.0):")
        print(variable[cols].to_string(index=False))

    print()

    # 2. Linkage Analysis (Mutual Information)
    print("2. Linkage Analysis (Mutual Information)")
    print("-" * 40)
    print("Computing pairwise residue linkage (co-evolution)...")

    linkage = processor.compute_linkage(
        sequences=ALIGNED_SEQUENCES,
        normalize=True,
        top_k=10,
        min_observations=3,
    )

    if not linkage.empty:
        print(f"   Found {len(linkage)} position pairs with significant linkage\n")
        print("   Top linkage pairs:")
        print(linkage.to_string(index=False))

        # Check if positions 4 and 13 show linkage (we designed them to co-vary)
        pos_4_13 = linkage[
            ((linkage["pos_i"] == 4) & (linkage["pos_j"] == 13)) |
            ((linkage["pos_i"] == 13) & (linkage["pos_j"] == 4))
        ]
        if not pos_4_13.empty:
            print("\n   Note: Positions 4 and 13 were designed to co-vary in sample data")
    else:
        print("   No significant linkage pairs found (need more sequence variation)")

    print()

    # 3. Sequence Metadata Extraction
    print("3. Sequence Metadata Extraction")
    print("-" * 40)

    # Register sequences first to use get_sequence_metadata
    dataset_name = "analysis_demo_sequences"
    processor.save_sequences(
        ALIGNED_SEQUENCES,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )

    metadata = processor.get_sequence_metadata(list(ALIGNED_SEQUENCES.keys()))

    if metadata is not None and not metadata.empty:
        print(f"   Extracted metadata for {len(metadata)} sequences\n")
        print("   Sequence properties:")
        # Show available columns
        print(f"   Available columns: {list(metadata.columns)}")
        print(metadata.to_string())
    else:
        print("   Metadata extraction not available or returned empty")

    print()

    # 4. Store conservation results
    print("4. Storing Analysis Results")
    print("-" * 40)

    # Store conservation as a named result
    stored_conservation = processor.compute_conservation(
        sequences=ALIGNED_SEQUENCES,
        store_result=True,
        result_name="demo_conservation_analysis",
    )
    print(f"   Stored conservation analysis as 'demo_conservation_analysis'")

    # Store linkage as a named result
    stored_linkage = processor.compute_linkage(
        sequences=ALIGNED_SEQUENCES,
        store_result=True,
        result_name="demo_linkage_analysis",
        top_k=20,
    )
    print(f"   Stored linkage analysis as 'demo_linkage_analysis'")

    print()
    print("=== SequenceProcessor Analytical Demo Complete ===")


if __name__ == "__main__":
    main()

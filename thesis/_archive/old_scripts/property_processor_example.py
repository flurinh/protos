#!/usr/bin/env python3
"""PropertyProcessor Example: Associating spectral properties with cone opsins.

Demonstrates ProtOS capabilities:
- PropertyProcessor: Store ANY annotation/property for ANY registered entity
- Cross-processor data flow: Link properties to sequences from SequenceProcessor
- Partial coverage: Not all entities need property values (subset of 200 sequences)
- Type-flexible storage: numbers, strings, nested data

Question: "Can we associate experimental lambda_max values with our cone opsin dataset?"

DESIGN PRINCIPLE: The PropertyProcessor enables zero-configuration storage of
arbitrary data linked to any entity in the ProtOS registry. This example shows
how experimental spectral measurements can be associated with the 200-sequence
cone opsin dataset (100 SW + 100 LW) created by SequenceProcessor - even when
we only have measurements for a subset of sequences.

KEY INSIGHT: In real workflows (like the Red-Shift Mutation workflow), predicted
lambda_max values are stored using PropertyProcessor. This means predictions
become *associated data*, not just outputs - enabling queries like "which
sequences have predicted lambda_max > 550 nm?"
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.property import PropertyProcessor


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "property"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
SEQUENCE_OUTPUT_DIR = THESIS_DIR / "outputs" / "sequence"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Source dataset from SequenceProcessor
SOURCE_DATASET = "cone_opsin_diversity"
TABLE_NAME = "opsin_spectral_properties"

# Known spectral data for cone opsins (literature values)
# These UniProt IDs may or may not be in our 200-sequence dataset (SW + LW)
KNOWN_SPECTRAL_DATA = {
    # Short-wave (blue) opsins - λmax ~360-440 nm
    "P03999": {"name": "OPSB_HUMAN", "type": "short_wave", "lambda_max": 420, "species": "Human"},
    "P51491": {"name": "OPSB_MOUSE", "type": "short_wave", "lambda_max": 360, "species": "Mouse"},
    "P51490": {"name": "OPSB_BOVIN", "type": "short_wave", "lambda_max": 438, "species": "Bovine"},
    "O13092": {"name": "OPSB_SAIBB", "type": "short_wave", "lambda_max": 423, "species": "Marmoset"},
    "P60015": {"name": "OPSB_PANTR", "type": "short_wave", "lambda_max": 430, "species": "Chimpanzee"},
    "Q8HY69": {"name": "OPSB_MACFA", "type": "short_wave", "lambda_max": 430, "species": "Macaque"},
    # Long-wave (red) opsins - λmax ~555-565 nm
    "P04000": {"name": "OPSR_HUMAN", "type": "long_wave", "lambda_max": 560, "species": "Human"},
    "P28683": {"name": "OPSR_CEBAP", "type": "long_wave", "lambda_max": 561, "species": "Capuchin"},
    "P34989": {"name": "OPN1LW_PONPY", "type": "long_wave", "lambda_max": 563, "species": "Orangutan"},
}

TYPE_COLORS = {
    "short_wave": "#006d77",   # Teal (matching sequence colors)
    "long_wave": "#e9c46a",    # Gold
}


def wavelength_to_rgb(wavelength: float) -> str:
    """Convert wavelength (nm) to approximate RGB color."""
    if wavelength < 380:
        return "rgb(100, 0, 150)"  # UV -> violet
    elif wavelength < 440:
        r, g, b = -(wavelength - 440) / 60, 0, 1
    elif wavelength < 490:
        r, g, b = 0, (wavelength - 440) / 50, 1
    elif wavelength < 510:
        r, g, b = 0, 1, -(wavelength - 510) / 20
    elif wavelength < 580:
        r, g, b = (wavelength - 510) / 70, 1, 0
    elif wavelength < 645:
        r, g, b = 1, -(wavelength - 645) / 65, 0
    elif wavelength <= 700:
        r, g, b = 1, 0, 0
    else:
        return "rgb(50, 0, 0)"

    factor = 0.3 + 0.7 * min((wavelength - 380) / 40, (700 - wavelength) / 55, 1.0) if wavelength < 420 or wavelength > 645 else 1.0
    return f"rgb({int(255*r*factor)}, {int(255*g*factor)}, {int(255*b*factor)})"


def main() -> int:
    """Run the PropertyProcessor example."""
    print("=" * 70)
    print("PROPERTY PROCESSOR EXAMPLE")
    print("Associating Spectral Properties with Cone Opsin Dataset")
    print("=" * 70)

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load the cone opsin dataset from SequenceProcessor
    # -------------------------------------------------------------------------
    print("\n[1] Loading cone opsin dataset...")
    seq_proc = SequenceProcessor()

    sequences = seq_proc.load_dataset(SOURCE_DATASET)

    if not sequences:
        print(f"  Dataset '{SOURCE_DATASET}' not found!")
        print("  Please run sequence_processor_example.py first.")
        print("  Continuing with spectral data demonstration...")
        sequences = {}

    print(f"  Loaded {len(sequences)} sequences from {SOURCE_DATASET}")

    # Load annotations to get accession IDs
    annotations_file = SEQUENCE_OUTPUT_DIR / "opsin_annotations.json"
    if annotations_file.exists():
        with open(annotations_file) as f:
            annotations = json.load(f)
        print(f"  Loaded annotations for {len(annotations)} sequences")
    else:
        annotations = {}

    # -------------------------------------------------------------------------
    # Step 2: Match spectral data to sequences in our dataset
    # -------------------------------------------------------------------------
    print("\n[2] Matching spectral data to dataset sequences...")
    prop_proc = PropertyProcessor()

    # Find which sequences in our dataset have known spectral values
    matched_properties = []
    unmatched_accessions = set(KNOWN_SPECTRAL_DATA.keys())

    for acc_id in annotations.keys():
        if acc_id in KNOWN_SPECTRAL_DATA:
            spectral_info = KNOWN_SPECTRAL_DATA[acc_id]
            matched_properties.append({
                "accession": acc_id,
                "entity_name": spectral_info["name"],
                "opsin_type": spectral_info["type"],
                "lambda_max": spectral_info["lambda_max"],
                "species": spectral_info["species"],
                "source": "literature",
            })
            unmatched_accessions.discard(acc_id)
            print(f"  MATCHED: {acc_id} -> {spectral_info['name']} (λmax={spectral_info['lambda_max']} nm)")

    print(f"\n  Sequences in dataset: {len(annotations)}")
    print(f"  Sequences with known λmax: {len(matched_properties)}")
    print(f"  Coverage: {len(matched_properties)}/{len(annotations)} ({100*len(matched_properties)/max(len(annotations),1):.1f}%)")

    if unmatched_accessions:
        print(f"\n  Note: {len(unmatched_accessions)} spectral values not in dataset:")
        for acc in list(unmatched_accessions)[:3]:
            info = KNOWN_SPECTRAL_DATA[acc]
            print(f"    - {info['name']} ({acc})")

    # -------------------------------------------------------------------------
    # Step 3: Record properties using ProtOS PropertyProcessor
    # -------------------------------------------------------------------------
    print("\n[3] Recording spectral properties...")

    if matched_properties:
        property_rows = []
        for prop in matched_properties:
            property_rows.append({
                "entity_name": prop["entity_name"],
                "scope": [{"format": "sequence", "name": prop["entity_name"]}],
                "accession": prop["accession"],
                "opsin_type": prop["opsin_type"],
                "lambda_max": prop["lambda_max"],
                "species": prop["species"],
                "source": prop["source"],
            })

        # ProtOS records properties with entity association
        properties_df = prop_proc.record_properties(
            TABLE_NAME,
            property_rows,
            metadata={
                "description": "Cone opsin spectral properties",
                "source_dataset": SOURCE_DATASET,
                "coverage": f"{len(matched_properties)}/{len(annotations)}",
            },
            allow_create=True,
        )
        print(f"  Recorded {len(property_rows)} property entries")
    else:
        # Create properties from known data even if dataset isn't loaded
        property_rows = []
        for acc_id, info in KNOWN_SPECTRAL_DATA.items():
            property_rows.append({
                "entity_name": info["name"],
                "scope": [{"format": "sequence", "name": info["name"]}],
                "accession": acc_id,
                "opsin_type": info["type"],
                "lambda_max": info["lambda_max"],
                "species": info["species"],
                "source": "literature",
            })

        properties_df = prop_proc.record_properties(
            TABLE_NAME, property_rows,
            metadata={"description": "Cone opsin spectral properties"},
            allow_create=True,
        )
        print(f"  Recorded {len(property_rows)} property entries (from reference data)")

    # -------------------------------------------------------------------------
    # Step 4: Demonstrate property queries
    # -------------------------------------------------------------------------
    print("\n[4] Querying properties...")

    # Query by opsin type
    for opsin_type in ["short_wave", "long_wave"]:
        type_props = prop_proc.filter_by_property(TABLE_NAME, "opsin_type", lambda x: x == opsin_type)
        if len(type_props) > 0:
            mean_lambda = type_props["lambda_max"].mean()
            print(f"  {opsin_type.replace('_', ' ').title()}: n={len(type_props)}, mean λmax={mean_lambda:.0f} nm")

    # -------------------------------------------------------------------------
    # Step 5: Create spectral sensitivity visualization
    # -------------------------------------------------------------------------
    print("\n[5] Creating spectral sensitivity visualization...")
    import plotly.graph_objects as go

    wavelengths = np.linspace(300, 700, 400)

    fig = go.Figure()

    # Visible spectrum background
    for wl in range(350, 701, 5):
        fig.add_vrect(x0=wl-2.5, x1=wl+2.5, fillcolor=wavelength_to_rgb(wl), opacity=0.25, line_width=0)

    # Sensitivity curves (Gaussian approximation)
    for _, row in properties_df.iterrows():
        peak = row["lambda_max"]
        # Bandwidth varies by opsin type (UV opsins typically narrower)
        bandwidth = 25 if peak < 400 else 35
        sensitivity = np.exp(-0.5 * ((wavelengths - peak) / bandwidth) ** 2)

        # Shortened name for legend
        short_name = row['entity_name'].split('_')[0] if '_' in row['entity_name'] else row['entity_name']
        fig.add_trace(go.Scatter(
            x=wavelengths,
            y=sensitivity,
            mode="lines",
            name=f"{short_name} ({int(peak)})",
            line=dict(color=TYPE_COLORS.get(row["opsin_type"], "#888888"), width=2.5),
            hovertemplate=f"{row['entity_name']}<br>λmax: {peak} nm<br>Type: {row['opsin_type']}<extra></extra>",
        ))

    fig.update_layout(
        xaxis_title="Wavelength (nm)",
        yaxis_title="Sensitivity",
        height=450,
        width=800,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(x=1.02, y=0.99, font=dict(size=8)),
        margin=dict(t=30, b=50, r=150),
    )
    fig.update_xaxes(showgrid=False, range=[320, 680], title_font_size=10)
    fig.update_yaxes(showgrid=False, range=[0, 1.1], title_font_size=10)
    fig.write_image(str(FIGURES_DIR / "property_sensitivity.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'property_sensitivity.png'}")

    # Save summary
    properties_df.to_csv(OUTPUT_DIR / "spectral_properties.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'spectral_properties.csv'}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)

    print(f"\nKEY CONCEPT: PropertyProcessor associates data with registered entities")
    print(f"\n  Source dataset: {SOURCE_DATASET} ({len(annotations)} sequences)")
    print(f"  Sequences with λmax data: {len(matched_properties)}")
    print(f"  Coverage: {len(matched_properties)}/{len(annotations)} sequences have spectral values")

    print(f"\n  IMPORTANT: Not all entities need property values!")
    print(f"  - Properties link to entities by name")
    print(f"  - Query: 'Which sequences have λmax > 500 nm?'")
    print(f"  - In Red-Shift workflow: predicted λmax stored as properties")
    print(f"  - Predictions become queryable associated data, not just outputs")

    print(f"\nOutputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

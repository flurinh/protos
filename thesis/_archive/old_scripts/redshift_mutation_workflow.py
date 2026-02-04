#!/usr/bin/env python3
"""Prey Animal Vision Enhancement - Red-Shift Mutation Screen.

Biological Question:
Tigers are orange because bovid prey (cattle, gazelles, antelopes) are dichromats -
they have only two cone types (SWS blue ~450nm, LWS yellow-green ~555nm) and cannot
distinguish orange from green. To the gazelle, the tiger's orange fur blends
perfectly with green grass.

Goal: Red-shift the bovine SWS (blue) opsin toward green (~530nm) to create a
"trichromatic-like" sensitivity that could help prey animals distinguish orange
predators from green vegetation.

Demonstrates cross-processor composition:
- SequenceProcessor: Load bovine SWS opsin, generate mutants, annotate with GRN
- EmbeddingProcessor: Generate Ankh-large embeddings
- ModelManager: Run LAMBDA predictions via Docker
- PropertyProcessor: Store predicted lambda_max values

Question: "Which mutations would red-shift bovine blue opsin toward green?"
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "redshift_screen"
FIGURES_DIR = THESIS_DIR / "workflows" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.property import PropertyProcessor
from protos.models.model_manager import ModelManager
from protos.io.ingest.sequence_loader import SequenceLoader


# =============================================================================
# Configuration
# =============================================================================
# Bovine short-wave-sensitive opsin 1 (Blue cone pigment)
# This represents dichromatic prey animal vision (shared by cattle, gazelles, antelopes)
WILD_TYPE_UNIPROT = "P51490"
WILD_TYPE_NAME = "OPSB_BOVIN"
WILD_TYPE_LAMBDA = 451  # nm (experimental blue cone)

# Target: shift toward green (~530nm) for orange/green discrimination
TARGET_LAMBDA = 530  # nm (green sensitivity)

PROTEIN_FAMILY = "opsin"
GRN_REFERENCE_TABLE = "vpod1_2"
GRN_PROTEIN_FAMILY = "gpcr_a"

# Spectral tuning mutation sites for SWS opsins
# Based on comparative studies of SWS1 opsin evolution in vertebrates
# Key sites known to affect SWS1 lambda_max:
MUTATION_SITES = {
    # Site 86: Major tuning site in SWS1 opsins
    # F86Y causes ~7nm red-shift, F86S causes ~40nm red-shift
    86: {"WT": "F", "mutations": ["Y", "S", "C"], "region": "spectral_tuning"},

    # Site 90: Important for UV vs violet sensitivity
    90: {"WT": "S", "mutations": ["C", "A"], "region": "retinal_pocket"},

    # Site 93: Affects chromophore environment
    93: {"WT": "T", "mutations": ["I", "V"], "region": "retinal_pocket"},

    # Site 114: Counterion region
    114: {"WT": "A", "mutations": ["G", "S"], "region": "counterion"},

    # Site 118: Near Schiff base
    118: {"WT": "T", "mutations": ["A", "S"], "region": "schiff_base"},

    # Site 261: TM6 - known spectral shift site
    261: {"WT": "F", "mutations": ["Y", "W"], "region": "TM6"},

    # Site 265: Key site in all opsins
    265: {"WT": "W", "mutations": ["Y"], "region": "retinal_contact"},

    # Site 292: TM7 position
    292: {"WT": "A", "mutations": ["S"], "region": "TM7"},
}

EMBEDDING_MODEL = "ankh_large"
LAMBDA_RUN_ID = "007061"


def generate_mutant_sequences(wt_sequence: str, mutation_sites: dict) -> Dict[str, Tuple[str, str]]:
    """Generate mutant sequences from wild-type."""
    mutants = {}

    for pos, site_info in mutation_sites.items():
        wt_aa = site_info["WT"]

        if pos > len(wt_sequence):
            print(f"  Warning: Position {pos} out of range")
            continue

        actual_aa = wt_sequence[pos - 1]
        if actual_aa != wt_aa:
            print(f"  Note: Position {pos} is {actual_aa} (expected {wt_aa})")
            wt_aa = actual_aa

        for mut_aa in site_info["mutations"]:
            if mut_aa == actual_aa:
                continue

            mut_seq = wt_sequence[:pos-1] + mut_aa + wt_sequence[pos:]
            mut_id = f"{actual_aa}{pos}{mut_aa}"
            mut_desc = f"{actual_aa}{pos}{mut_aa} ({site_info['region']})"
            mutants[mut_id] = (mut_seq, mut_desc)

    return mutants


def main() -> int:
    """Run the prey vision enhancement mutation screen."""
    print("=" * 70)
    print("PREY ANIMAL VISION ENHANCEMENT - RED-SHIFT MUTATION SCREEN")
    print("=" * 70)
    print()
    print("BIOLOGICAL CONTEXT:")
    print("  Tigers are orange because their prey (gazelles, cattle) are dichromats.")
    print("  These animals cannot distinguish orange from green - the tiger hides")
    print("  in plain sight against green vegetation.")
    print()
    print(f"  Wild-type: Bovine SWS opsin ({WILD_TYPE_NAME})")
    print(f"  Current sensitivity: {WILD_TYPE_LAMBDA} nm (blue)")
    print(f"  Target sensitivity: {TARGET_LAMBDA} nm (green)")
    print(f"  Required shift: +{TARGET_LAMBDA - WILD_TYPE_LAMBDA} nm")
    print()
    print("  GOAL: Find mutations that red-shift blue opsin toward green,")
    print("        enabling orange/green discrimination.")
    print("=" * 70)

    # Initialize ProtOS
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    # -------------------------------------------------------------------------
    # Step 1: Load bovine SWS opsin sequence
    # -------------------------------------------------------------------------
    print("\n[1] Loading bovine blue cone opsin (SWS1)...")
    seq_proc = SequenceProcessor()
    loader = SequenceLoader(processor=seq_proc)

    loader.download_and_register(
        f"uniprot:{WILD_TYPE_UNIPROT}",
        name=WILD_TYPE_NAME,
        materialize_entities=True,
    )

    wt_sequence = seq_proc.load_entity(WILD_TYPE_NAME)
    if not wt_sequence:
        print("  ERROR: Failed to load sequence")
        return 1
    print(f"  Loaded: {len(wt_sequence)} residues")

    # -------------------------------------------------------------------------
    # Step 2: Generate spectral tuning mutants
    # -------------------------------------------------------------------------
    print("\n[2] Generating spectral tuning mutants...")
    mutants = generate_mutant_sequences(wt_sequence, MUTATION_SITES)
    print(f"  Generated {len(mutants)} mutants at {len(MUTATION_SITES)} sites")

    # Build sequences dict
    all_sequences = {WILD_TYPE_NAME: wt_sequence}
    mutation_info = {WILD_TYPE_NAME: ("WT", "Wild-type (blue, 451nm)")}

    for mut_id, (mut_seq, mut_desc) in mutants.items():
        seq_name = f"{WILD_TYPE_NAME}_{mut_id}"
        all_sequences[seq_name] = mut_seq
        mutation_info[seq_name] = (mut_id, mut_desc)

    print(f"  Total sequences: {len(all_sequences)}")

    # -------------------------------------------------------------------------
    # Step 3: Save sequences as dataset
    # -------------------------------------------------------------------------
    print("\n[3] Saving sequence dataset...")
    dataset_name = "prey_vision_screen"

    seq_proc.save_sequences(
        all_sequences,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )
    print(f"  Dataset: {dataset_name}")

    # -------------------------------------------------------------------------
    # Step 4: Annotate with GRN
    # -------------------------------------------------------------------------
    print("\n[4] Annotating with GRN (GPCRdb numbering)...")
    grn_table_name = f"{dataset_name}_grn"

    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table=GRN_REFERENCE_TABLE,
        protein_family=GRN_PROTEIN_FAMILY,
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )
    print(f"  GRN table: {grn_table_name}")

    # -------------------------------------------------------------------------
    # Step 5: Generate Ankh-large embeddings
    # -------------------------------------------------------------------------
    print("\n[5] Generating Ankh-large embeddings...")
    embedding_dataset_name = f"{dataset_name}__ankh_large__per_residue"

    emb_proc = EmbeddingProcessor(model_name=EMBEDDING_MODEL, device="cpu")
    emb_proc.embed_sequences(
        all_sequences,
        embedding_type="per_residue",
        save_dataset=embedding_dataset_name,
    )
    print(f"  Embeddings: {embedding_dataset_name}")

    # -------------------------------------------------------------------------
    # Step 6: Run LAMBDA predictions via Docker
    # -------------------------------------------------------------------------
    print("\n[6] Running LAMBDA spectral predictions via Docker...")

    manager = ModelManager(data_root=data_root)

    invocation = manager.prepare(
        "lambda",
        inputs={
            "sequence_dataset": dataset_name,
            "grn_table": grn_table_name,
            "embedding_dataset": embedding_dataset_name,
            "protein_family": GRN_PROTEIN_FAMILY,
        },
        config={
            "run_id": LAMBDA_RUN_ID,
            "batch_size": 4,
            "embedding_model": EMBEDDING_MODEL,
            "embedding_type": "per_residue",
            "use_docker": True,
            "use_gpu": False,
            "job_name": f"lambda_{dataset_name}",
        },
    )

    job = invocation.job
    if job is None:
        print("  ERROR: LAMBDA adapter did not produce a Docker job")
        return 1

    print(f"  Running Docker container...")
    proc = subprocess.run(
        job.command,
        cwd=job.working_dir,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        print(f"  ERROR: Docker failed (exit {proc.returncode})")
        print(f"  stderr: {proc.stderr[:500]}")
        return 1

    print("  LAMBDA predictions complete")

    # -------------------------------------------------------------------------
    # Step 7: Load and analyze predictions
    # -------------------------------------------------------------------------
    print("\n[7] Analyzing predictions...")
    output_dir = Path(job.metadata.get("output_dir"))
    predictions_path = output_dir / "predictions.csv"

    if not predictions_path.exists():
        print(f"  ERROR: Predictions not found: {predictions_path}")
        return 1

    predictions_df = pd.read_csv(predictions_path)
    print(f"  Loaded {len(predictions_df)} predictions")

    # Find lambda column
    lambda_col = None
    for col in predictions_df.columns:
        if "11cis" in col.lower() or "lambda" in col.lower():
            lambda_col = col
            break
    if lambda_col is None:
        numeric_cols = predictions_df.select_dtypes(include=[np.number]).columns
        lambda_col = numeric_cols[0] if len(numeric_cols) > 0 else None
    if lambda_col is None:
        print("  ERROR: No prediction column found")
        return 1

    # Get wild-type prediction
    wt_pred_row = predictions_df[predictions_df["protein_id"] == WILD_TYPE_NAME]
    wt_pred_lambda = wt_pred_row[lambda_col].values[0] if len(wt_pred_row) > 0 else WILD_TYPE_LAMBDA

    # Build results
    results = []
    for _, row in predictions_df.iterrows():
        protein_id = row.get("protein_id", "")
        pred_lambda = row[lambda_col]
        mut_id, mut_desc = mutation_info.get(protein_id, ("unknown", "unknown"))
        shift = pred_lambda - wt_pred_lambda
        toward_green = TARGET_LAMBDA - WILD_TYPE_LAMBDA
        progress = (shift / toward_green * 100) if toward_green > 0 else 0

        results.append({
            "sequence_id": protein_id,
            "mutation": mut_id,
            "description": mut_desc,
            "lambda_pred": pred_lambda,
            "shift_nm": shift,
            "progress_%": progress,
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("shift_nm", ascending=False)
    results_df.to_csv(OUTPUT_DIR / "prey_vision_screen_results.csv", index=False)

    # -------------------------------------------------------------------------
    # Step 8: Report results
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("MUTATION SCREEN RESULTS")
    print("=" * 70)
    print(f"Wild-type predicted: {wt_pred_lambda:.1f} nm")
    print(f"Target (green): {TARGET_LAMBDA} nm")
    print(f"Required shift: +{TARGET_LAMBDA - wt_pred_lambda:.1f} nm")
    print()
    print(f"{'Mutation':<12} {'λmax (nm)':<12} {'Shift':<10} {'Progress':<12} {'Site'}")
    print("-" * 70)

    mutant_results = results_df[results_df["mutation"] != "WT"]
    for _, row in mutant_results.head(12).iterrows():
        region = row["description"].split("(")[-1].rstrip(")") if "(" in row["description"] else ""
        progress_str = f"{row['progress_%']:.0f}%" if row['shift_nm'] > 0 else "-"
        print(f"{row['mutation']:<12} {row['lambda_pred']:<12.1f} {row['shift_nm']:+.1f} nm    {progress_str:<12} {region}")

    # -------------------------------------------------------------------------
    # Step 9: Visualize results
    # -------------------------------------------------------------------------
    print("\n[8] Creating visualization...")
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "Spectral Shift Distribution",
            "Progress Toward Green Sensitivity"
        ),
    )

    # Histogram
    shifts = mutant_results["shift_nm"]
    fig.add_trace(
        go.Histogram(x=shifts, nbinsx=12, marker_color="#2ecc71", name="Mutations"),
        row=1, col=1,
    )
    fig.add_vline(x=0, line_dash="dash", line_color="black", row=1, col=1)

    # Bar chart - progress toward green
    top_mutants = mutant_results.head(10)
    colors = ["#27ae60" if s > 0 else "#e74c3c" for s in top_mutants["shift_nm"]]

    fig.add_trace(
        go.Bar(
            x=top_mutants["mutation"],
            y=top_mutants["shift_nm"],
            marker_color=colors,
            name="Shift",
        ),
        row=1, col=2,
    )

    # Add target line
    target_shift = TARGET_LAMBDA - wt_pred_lambda
    fig.add_hline(y=target_shift, line_dash="dot", line_color="green",
                  annotation_text=f"Target (+{target_shift:.0f}nm)", row=1, col=2)

    fig.update_layout(
        title=dict(
            text="Prey Vision Enhancement: Red-Shifting Blue Opsin Toward Green<br>"
                 "<sub>Goal: Help dichromatic prey animals distinguish orange predators from green vegetation</sub>",
            x=0.5,
        ),
        height=500,
        width=1100,
        showlegend=False,
    )
    fig.update_xaxes(title_text="Shift (nm)", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_xaxes(title_text="Mutation", row=1, col=2)
    fig.update_yaxes(title_text="Shift (nm)", row=1, col=2)

    fig.write_html(str(FIGURES_DIR / "prey_vision_enhancement.html"))
    print(f"  Saved: {FIGURES_DIR / 'prey_vision_enhancement.html'}")

    # Compute summary stats needed for Step 9
    n_redshift = len(shifts[shifts > 0])
    n_blueshift = len(shifts[shifts < 0])
    best_shift = shifts.max()
    best_mut = mutant_results.iloc[0]["mutation"] if len(mutant_results) > 0 else "N/A"

    # -------------------------------------------------------------------------
    # Step 9: Simulate dichromatic vs enhanced vision
    # -------------------------------------------------------------------------
    print("\n[9] Simulating prey animal vision...")
    from PIL import Image

    tiger_path = FIGURES_DIR / "tiger_green.jpg"
    if tiger_path.exists():
        img = Image.open(tiger_path)
        img_array = np.array(img, dtype=np.float32) / 255.0

        # Cone sensitivity peaks (nm)
        sws_wt = WILD_TYPE_LAMBDA  # 451 nm (blue)
        sws_mut = wt_pred_lambda + best_shift  # ~460 nm (shifted toward green)
        lws = 555  # Long-wave (yellow-green) - unchanged

        def simulate_dichromat_vision(rgb_img: np.ndarray, sws_peak: float, lws_peak: float) -> np.ndarray:
            """Simulate dichromatic vision based on cone sensitivity peaks.

            Bovids have SWS (blue) and LWS (yellow-green) cones only.
            They cannot distinguish colors that differ only in the M-cone (green) region.
            Orange and green appear similar because both primarily stimulate LWS.
            """
            # Approximate cone response weights based on peak wavelengths
            # RGB roughly maps to: R~600nm, G~550nm, B~450nm

            # SWS cone response (primarily blue channel, some green)
            sws_weight_b = 1.0
            sws_weight_g = 0.1 * max(0, 1 - abs(sws_peak - 500) / 100)  # Small green contribution if shifted
            sws_weight_r = 0.0

            # LWS cone response (yellow-green: R and G channels)
            lws_weight_r = 0.7  # Strong red contribution
            lws_weight_g = 1.0  # Strong green contribution
            lws_weight_b = 0.0

            # Calculate cone responses
            r, g, b = rgb_img[:,:,0], rgb_img[:,:,1], rgb_img[:,:,2]

            sws_response = sws_weight_r * r + sws_weight_g * g + sws_weight_b * b
            lws_response = lws_weight_r * r + lws_weight_g * g + lws_weight_b * b

            # Normalize responses
            sws_response = np.clip(sws_response, 0, 1)
            lws_response = np.clip(lws_response, 0, 1)

            # Map back to RGB for visualization
            # SWS -> Blue channel, LWS -> Yellow (R+G)
            # This shows what colors are distinguishable
            out_r = lws_response * 0.8
            out_g = lws_response * 0.7 + sws_response * 0.1
            out_b = sws_response * 0.6

            return np.stack([out_r, out_g, out_b], axis=-1)

        def simulate_enhanced_dichromat(rgb_img: np.ndarray, sws_shift: float, full_shift: float = 80) -> np.ndarray:
            """Simulate dichromatic vision with red-shifted SWS opsin.

            The shifted SWS cone provides slightly better green sensitivity,
            improving orange/green discrimination.
            """
            r, g, b = rgb_img[:,:,0], rgb_img[:,:,1], rgb_img[:,:,2]

            # Enhanced SWS with better green sensitivity
            green_boost = min(sws_shift / full_shift, 1.0)  # Proportion of shift achieved
            sws_response = b + green_boost * g * 0.4  # Gains green sensitivity
            lws_response = 0.7 * r + g

            sws_response = np.clip(sws_response, 0, 1)
            lws_response = np.clip(lws_response, 0, 1)

            # Better separation between orange (high R) and green (high G)
            orange_signal = r - g * 0.5  # Orange has more R than G
            orange_signal = np.clip(orange_signal, 0, 1)

            # With shifted SWS, orange appears slightly more distinct
            out_r = lws_response * 0.8 + orange_signal * green_boost * 0.4
            out_g = lws_response * 0.7 + sws_response * 0.15
            out_b = sws_response * 0.5

            return np.stack([out_r, out_g, out_b], axis=-1)

        # Generate simulations
        dichromat_view = simulate_dichromat_vision(img_array, sws_wt, lws)
        enhanced_view = simulate_enhanced_dichromat(img_array, best_shift)
        target_view = simulate_enhanced_dichromat(img_array, 80, 80)  # Full green shift

        # Create comparison figure - 2x2 layout
        fig_vision = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "Human Vision (Trichromat)",
                f"Current Prey Vision (Dichromat)<br><sub>SWS={sws_wt:.0f}nm + LWS={lws}nm - Tiger is camouflaged!</sub>",
                f"With {best_mut} Mutation ({best_shift:+.0f}nm)<br><sub>SWS={sws_mut:.0f}nm - {best_shift/80*100:.0f}% progress</sub>",
                f"Target: Full Green Sensitivity<br><sub>SWS=530nm - Tiger becomes visible!</sub>",
            ),
            horizontal_spacing=0.02,
            vertical_spacing=0.08,
        )

        # Original (human trichromatic view)
        fig_vision.add_trace(
            go.Image(z=(img_array * 255).astype(np.uint8)),
            row=1, col=1,
        )

        # Dichromat view (current prey vision)
        fig_vision.add_trace(
            go.Image(z=(np.clip(dichromat_view, 0, 1) * 255).astype(np.uint8)),
            row=1, col=2,
        )

        # Enhanced dichromat view (with mutant opsin)
        fig_vision.add_trace(
            go.Image(z=(np.clip(enhanced_view, 0, 1) * 255).astype(np.uint8)),
            row=2, col=1,
        )

        # Target view (full green sensitivity)
        fig_vision.add_trace(
            go.Image(z=(np.clip(target_view, 0, 1) * 255).astype(np.uint8)),
            row=2, col=2,
        )

        fig_vision.update_layout(
            title=dict(
                text="Tiger Camouflage: Engineering Prey Vision to See the Predator<br>"
                     f"<sub>Dichromats confuse orange with green | {best_mut} provides {best_shift/80*100:.0f}% "
                     f"progress toward green sensitivity</sub>",
                x=0.5,
            ),
            height=800,
            width=1100,
            showlegend=False,
        )

        # Remove axes
        for i in range(1, 3):
            for j in range(1, 3):
                fig_vision.update_xaxes(visible=False, row=i, col=j)
                fig_vision.update_yaxes(visible=False, row=i, col=j)

        fig_vision.write_html(str(FIGURES_DIR / "prey_vision_simulation.html"))
        print(f"  Saved: {FIGURES_DIR / 'prey_vision_simulation.html'}")

        # Also save as static images for comparison
        Image.fromarray((np.clip(dichromat_view, 0, 1) * 255).astype(np.uint8)).save(
            FIGURES_DIR / "tiger_dichromat_view.jpg"
        )
        Image.fromarray((np.clip(enhanced_view, 0, 1) * 255).astype(np.uint8)).save(
            FIGURES_DIR / "tiger_enhanced_view.jpg"
        )
        Image.fromarray((np.clip(target_view, 0, 1) * 255).astype(np.uint8)).save(
            FIGURES_DIR / "tiger_target_view.jpg"
        )
        print(f"  Saved: tiger_dichromat_view.jpg, tiger_enhanced_view.jpg, tiger_target_view.jpg")
    else:
        print(f"  Tiger image not found: {tiger_path}")

    # Summary (stats computed above before Step 9)
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Mutants screened: {len(mutants)}")
    print(f"Red-shifted (toward green): {n_redshift}")
    print(f"Blue-shifted: {n_blueshift}")
    print(f"Best mutation: {best_mut} ({best_shift:+.1f} nm)")
    print(f"Progress to green: {best_shift / (TARGET_LAMBDA - wt_pred_lambda) * 100:.1f}%")
    print()
    print("BIOLOGICAL INTERPRETATION:")
    if best_shift > 20:
        print(f"  {best_mut} achieves significant red-shift toward green sensitivity!")
        print("  This could substantially improve orange/green discrimination.")
    elif best_shift > 0:
        print(f"  {best_mut} shows modest red-shift. Multiple mutations may be needed")
        print("  to achieve full trichromatic-like vision.")
    else:
        print("  No significant red-shift achieved with single mutations.")
        print("  Combinatorial mutations may be required.")
    print()
    print(f"Outputs: {OUTPUT_DIR}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

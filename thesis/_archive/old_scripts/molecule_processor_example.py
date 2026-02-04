#!/usr/bin/env python3
"""MoleculeProcessor Example: Ligand similarity search and ranking.

Demonstrates ProtOS capabilities:
- MoleculeProcessor: Entity storage and property calculation
- Database querying for compound similarity search
- RDKit fingerprint-based similarity ranking

Question: "Which compounds are similar to carazolol (the inverse agonist in 2RH1)?"

This example ties directly into the GPCR binding workflow by finding similar
compounds to carazolol - the inverse agonist from the beta-2 adrenergic receptor
structure 2RH1. Carazolol is a beta-blocker that stabilizes the inactive state.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.molecule import MoleculeProcessor

# RDKit for similarity calculations
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs, Descriptors
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("Warning: RDKit not available.")


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "molecule"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Query compound: Carazolol - the inverse agonist in 2RH1 (from GPCR binding workflow)
# Carazolol is a non-selective beta-blocker with carbazole scaffold
QUERY_COMPOUND = {
    "name": "carazolol",
    "smiles": "CC(C)NCC(O)COc1cccc2[nH]c3ccccc3c12",
    "source": "2RH1 crystallographic ligand",
    "target": "ADRB2",
    "type": "inverse_agonist",
    "pdb_ligand_id": "CAU",
}

# Compound database: Beta-blockers, adrenergic ligands, and related drugs
# This simulates querying a chemical database like ChEMBL for adrenergic compounds
COMPOUND_DATABASE = {
    # ==========================================================================
    # NON-SELECTIVE BETA-BLOCKERS (same class as carazolol)
    # ==========================================================================
    "propranolol": {"smiles": "CC(C)NCC(O)COc1cccc2ccccc12", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "nadolol": {"smiles": "CC(C)NCC(O)c1ccc2c(c1)C(O)C(O)CC2", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "timolol": {"smiles": "CC(C)(C)NCC(O)COc1nsnc1N1CCOCC1", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "pindolol": {"smiles": "CC(C)NCC(O)COc1cccc2[nH]ccc12", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "sotalol": {"smiles": "CC(C)NCC(O)c1ccc(NS(C)(=O)=O)cc1", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "labetalol": {"smiles": "CC(CCc1ccccc1)NCC(O)c1ccc(O)c(C(N)=O)c1", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2/ADRA1"},
    "carvedilol": {"smiles": "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2/ADRA1"},
    "oxprenolol": {"smiles": "CC(C)NCC(O)COc1ccccc1OCC=C", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "penbutolol": {"smiles": "CC(C)(C)NCC(O)COc1ccccc1C1CCCC1", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "levobunolol": {"smiles": "CC(C)(C)NCC(O)COc1ccc2CC(=O)CCc2c1", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "metipranolol": {"smiles": "CC(C)NCC(O)COc1ccc(C(C)=O)c(OC(C)=O)c1OC(C)=O", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},
    "bupranolol": {"smiles": "CC(C)NCC(O)COc1ccc(Cl)cc1Cl", "class": "nonselective_beta_blocker", "target": "ADRB1/ADRB2"},

    # ==========================================================================
    # BETA-1 SELECTIVE BLOCKERS (cardioselective)
    # ==========================================================================
    "metoprolol": {"smiles": "COCCc1ccc(OCC(O)CNC(C)C)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "atenolol": {"smiles": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "bisoprolol": {"smiles": "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "esmolol": {"smiles": "COC(=O)CCc1ccc(OCC(O)CNC(C)C)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "acebutolol": {"smiles": "CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(C(C)=O)c1", "class": "beta1_selective", "target": "ADRB1"},
    "betaxolol": {"smiles": "CC(C)NCC(O)COc1ccc(CCOCC2CC2)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "nebivolol": {"smiles": "OC(CNCC(O)c1ccc2OCCc2c1)c1ccc2OCCc2c1", "class": "beta1_selective", "target": "ADRB1"},
    "celiprolol": {"smiles": "CCN(CC)C(=O)Nc1ccc(OCC(O)CNC(C)(C)C)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "bevantolol": {"smiles": "COc1ccccc1OCCNCC(O)COc1ccc(C)cc1", "class": "beta1_selective", "target": "ADRB1"},
    "practolol": {"smiles": "CC(C)NCC(O)COc1ccc(NC(C)=O)cc1", "class": "beta1_selective", "target": "ADRB1"},

    # ==========================================================================
    # BETA-2 AGONISTS (opposing pharmacology - for comparison)
    # ==========================================================================
    "isoproterenol": {"smiles": "CC(C)NCC(O)c1ccc(O)c(O)c1", "class": "beta_agonist", "target": "ADRB1/ADRB2"},
    "salbutamol": {"smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "terbutaline": {"smiles": "CC(C)(C)NCC(O)c1cc(O)cc(O)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "formoterol": {"smiles": "COc1ccc(CC(C)NCC(O)c2ccc(O)c(NC=O)c2)cc1", "class": "beta2_agonist", "target": "ADRB2"},
    "salmeterol": {"smiles": "OCc1ccc(O)c(C(O)CNCCCCCCOCCCCc2ccccc2)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "fenoterol": {"smiles": "CC(Cc1ccc(O)cc1)NCC(O)c1ccc(O)c(O)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "clenbuterol": {"smiles": "CC(C)(C)NCC(O)c1cc(Cl)c(N)c(Cl)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "pirbuterol": {"smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)n1", "class": "beta2_agonist", "target": "ADRB2"},
    "procaterol": {"smiles": "CC(C)c1cc(C(O)CNC(C)C)c(O)c2ccccc12", "class": "beta2_agonist", "target": "ADRB2"},
    "indacaterol": {"smiles": "CCc1cc2c(cc1CC)C(O)C(O)CNCCCOc1cccc3[nH]ccc13", "class": "beta2_agonist", "target": "ADRB2"},
    "vilanterol": {"smiles": "OCCCCCCNCCc1ccc(O)c(CO)c1OCC(O)c1cccc(Cl)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "olodaterol": {"smiles": "COc1ccc(CC(C)NCC(O)c2ccc(O)c3OCCc23)cc1OC", "class": "beta2_agonist", "target": "ADRB2"},

    # ==========================================================================
    # NATURAL CATECHOLAMINES (endogenous ligands)
    # ==========================================================================
    "adrenaline": {"smiles": "CNCC(O)c1ccc(O)c(O)c1", "class": "catecholamine", "target": "ADRB1/ADRB2/ADRA"},
    "noradrenaline": {"smiles": "NCC(O)c1ccc(O)c(O)c1", "class": "catecholamine", "target": "ADRB1/ADRA"},
    "dopamine": {"smiles": "NCCc1ccc(O)c(O)c1", "class": "catecholamine", "target": "D1-D5"},
    "dobutamine": {"smiles": "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "class": "catecholamine", "target": "ADRB1"},

    # ==========================================================================
    # ALPHA-1 BLOCKERS (related GPCR pharmacology)
    # ==========================================================================
    "prazosin": {"smiles": "COc1cc2nc(N3CCN(C(=O)c4ccco4)CC3)nc(N)c2cc1OC", "class": "alpha1_blocker", "target": "ADRA1"},
    "doxazosin": {"smiles": "COc1cc2nc(N3CCN(C(=O)C4COc5ccccc5O4)CC3)nc(N)c2cc1OC", "class": "alpha1_blocker", "target": "ADRA1"},
    "terazosin": {"smiles": "COc1cc2nc(N3CCN(C(=O)C4CCCO4)CC3)nc(N)c2cc1OC", "class": "alpha1_blocker", "target": "ADRA1"},
    "tamsulosin": {"smiles": "COc1ccc(CCNC(C)COCCCC2ccc(OC)c(S(N)(=O)=O)c2)cc1", "class": "alpha1_blocker", "target": "ADRA1A"},
    "alfuzosin": {"smiles": "COc1cc2nc(N3CCN(C(=O)N4CCCC4)CC3)nc(N)c2cc1OC", "class": "alpha1_blocker", "target": "ADRA1"},
    "silodosin": {"smiles": "CCN(CC)CCN(Cc1ccccc1F)c1ccc(C(N)=O)c2ccc(OC)cc12", "class": "alpha1_blocker", "target": "ADRA1A"},

    # ==========================================================================
    # ALPHA-2 AGONISTS (related GPCR pharmacology)
    # ==========================================================================
    "clonidine": {"smiles": "Clc1cccc(Cl)c1NC1=NCCN1", "class": "alpha2_agonist", "target": "ADRA2"},
    "guanfacine": {"smiles": "NC(=N)NC(=O)Cc1cccc(Cl)c1Cl", "class": "alpha2_agonist", "target": "ADRA2A"},
    "brimonidine": {"smiles": "Brc1ccc2ncc(NC3=NCCN3)nc2c1", "class": "alpha2_agonist", "target": "ADRA2"},
    "dexmedetomidine": {"smiles": "Cc1cccc(C)c1C(C)c1cnc[nH]1", "class": "alpha2_agonist", "target": "ADRA2"},
    "methyldopa": {"smiles": "CC(N)(Cc1ccc(O)c(O)c1)C(=O)O", "class": "alpha2_agonist", "target": "ADRA2"},

    # ==========================================================================
    # CARBAZOLE DERIVATIVES (similar scaffold to carazolol)
    # ==========================================================================
    "carprofen": {"smiles": "CC(C(=O)O)c1ccc2[nH]c3ccc(Cl)cc3c2c1", "class": "carbazole", "target": "COX"},
    "rimcazole": {"smiles": "c1ccc2c(c1)c1ccccc1n2CCCN1CCN(c2ccccc2)CC1", "class": "carbazole", "target": "sigma"},
    "carvedilol_metabolite": {"smiles": "Oc1ccc2[nH]c3ccccc3c2c1OCCNCC(O)COc1ccccc1OC", "class": "carbazole", "target": "ADRB"},

    # ==========================================================================
    # INDOLE/INDOLINE DERIVATIVES (structural similarity)
    # ==========================================================================
    "pindolol": {"smiles": "CC(C)NCC(O)COc1cccc2[nH]ccc12", "class": "indole_beta_blocker", "target": "ADRB1/ADRB2"},
    "bopindolol": {"smiles": "CC(=O)Oc1cccc2[nH]ccc12.CC(C)NCC(O)CO", "class": "indole_beta_blocker", "target": "ADRB1/ADRB2"},
    "mepindolol": {"smiles": "Cc1ccc2[nH]ccc2c1OCC(O)CNC(C)C", "class": "indole_beta_blocker", "target": "ADRB1/ADRB2"},

    # ==========================================================================
    # OTHER GPCR LIGANDS (for comparison)
    # ==========================================================================
    "isoprenaline": {"smiles": "CC(C)NCC(O)c1ccc(O)c(O)c1", "class": "beta_agonist", "target": "ADRB"},
    "ritodrine": {"smiles": "CC(NCCc1ccc(O)cc1)C(O)c1ccc(O)cc1", "class": "beta2_agonist", "target": "ADRB2"},
    "bambuterol": {"smiles": "CN(C)C(=O)Oc1cc(C(O)CNC(C)(C)C)cc(OC(=O)N(C)C)c1", "class": "beta2_agonist", "target": "ADRB2"},
    "orciprenaline": {"smiles": "CC(C)NCC(O)c1cc(O)cc(O)c1", "class": "beta_agonist", "target": "ADRB"},
    "reproterol": {"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2CCNCC(O)c1cc(O)cc(O)c1", "class": "beta2_agonist", "target": "ADRB2"},
}


def calculate_fingerprint(smiles: str) -> Optional[object]:
    """Calculate Morgan fingerprint for a SMILES string."""
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def calculate_similarity(fp1, fp2) -> float:
    """Calculate Tanimoto similarity between two fingerprints."""
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def calculate_properties(smiles: str) -> Optional[Dict]:
    """Calculate molecular properties."""
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "rings": Descriptors.RingCount(mol),
    }


def main() -> int:
    """Run the MoleculeProcessor example."""
    print("=" * 70)
    print("MOLECULE PROCESSOR EXAMPLE")
    print("Ligand Similarity Search and Ranking")
    print("=" * 70)
    print(f"\nQuery: {QUERY_COMPOUND['name']} (from {QUERY_COMPOUND['source']})")
    print(f"Database: {len(COMPOUND_DATABASE)} compounds")

    if not HAS_RDKIT:
        print("\nError: RDKit required for this example")
        return 1

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Register query compound using ProtOS MoleculeProcessor
    # -------------------------------------------------------------------------
    print("\n[1] Registering query compound...")
    mol_proc = MoleculeProcessor()

    mol_proc.save_entity(
        name=QUERY_COMPOUND["name"],
        data={
            "smiles": QUERY_COMPOUND["smiles"],
            "kind": "query_ligand",
            "source": QUERY_COMPOUND["source"],
        },
        metadata={
            "target": QUERY_COMPOUND["target"],
            "type": QUERY_COMPOUND["type"],
        },
    )

    query_props = calculate_properties(QUERY_COMPOUND["smiles"])
    print(f"  {QUERY_COMPOUND['name']}: MW={query_props['mw']:.1f}, LogP={query_props['logp']:.2f}")
    print(f"  SMILES: {QUERY_COMPOUND['smiles']}")

    # Calculate query fingerprint
    query_fp = calculate_fingerprint(QUERY_COMPOUND["smiles"])

    # -------------------------------------------------------------------------
    # Step 2: Register database compounds with ProtOS
    # -------------------------------------------------------------------------
    print("\n[2] Registering compound database...")

    for name, info in COMPOUND_DATABASE.items():
        mol_proc.save_entity(
            name=name,
            data={
                "smiles": info["smiles"],
                "kind": "database_compound",
                "drug_class": info["class"],
            },
            metadata={
                "target": info["target"],
                "class": info["class"],
            },
        )

    print(f"  Registered {len(COMPOUND_DATABASE)} compounds")

    # -------------------------------------------------------------------------
    # Step 3: Similarity search and ranking
    # -------------------------------------------------------------------------
    print("\n[3] Calculating similarity to all database compounds...")

    import pandas as pd

    results = []
    for name, info in COMPOUND_DATABASE.items():
        fp = calculate_fingerprint(info["smiles"])
        similarity = calculate_similarity(query_fp, fp)
        props = calculate_properties(info["smiles"])

        results.append({
            "name": name,
            "smiles": info["smiles"],
            "class": info["class"],
            "target": info["target"],
            "similarity": similarity,
            "mw": props["mw"] if props else None,
            "logp": props["logp"] if props else None,
            "tpsa": props["tpsa"] if props else None,
            "rings": props["rings"] if props else None,
        })

    # Create DataFrame and rank by similarity
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("similarity", ascending=False)

    # Save full results
    results_df.to_csv(OUTPUT_DIR / "similarity_search_results.csv", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'similarity_search_results.csv'}")

    # -------------------------------------------------------------------------
    # Step 4: Display top hits
    # -------------------------------------------------------------------------
    print("\n[4] Top 20 most similar compounds:")
    print("-" * 70)

    top_20 = results_df.head(20)
    for i, (_, row) in enumerate(top_20.iterrows(), 1):
        print(f"  {i:2d}. {row['name']:25s} | Sim: {row['similarity']:.3f} | {row['class']}")

    # -------------------------------------------------------------------------
    # Step 5: Analyze similarity by drug class
    # -------------------------------------------------------------------------
    print("\n[5] Similarity by drug class:")
    print("-" * 70)

    class_stats = results_df.groupby("class").agg({
        "similarity": ["mean", "max", "count"],
    }).round(3)
    class_stats.columns = ["mean_sim", "max_sim", "count"]
    class_stats = class_stats.sort_values("mean_sim", ascending=False)

    for class_name, row in class_stats.iterrows():
        print(f"  {class_name:25s} | Mean: {row['mean_sim']:.3f} | Max: {row['max_sim']:.3f} | n={int(row['count'])}")

    # -------------------------------------------------------------------------
    # Step 6: Create visualizations
    # -------------------------------------------------------------------------
    print("\n[6] Creating visualizations...")

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Figure 1: Similarity distribution histogram
    fig = go.Figure()

    fig.add_trace(go.Histogram(
        x=results_df["similarity"],
        nbinsx=30,
        marker_color="#9467bd",  # Purple (molecule color)
        opacity=0.8,
    ))

    # Mark similarity thresholds
    fig.add_vline(x=0.7, line_dash="dash", line_color="#c5b0d5",
                  annotation_text=">0.7", annotation_position="top", annotation_font_size=9)
    fig.add_vline(x=0.5, line_dash="dash", line_color="#7f7f7f",
                  annotation_text=">0.5", annotation_position="top", annotation_font_size=9)

    fig.update_layout(
        xaxis_title="Tanimoto Similarity",
        yaxis_title="Count",
        height=380,
        width=600,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(t=40, b=50),
    )
    fig.update_xaxes(showgrid=False, range=[0, 1], title_font_size=10)
    fig.update_yaxes(showgrid=False, title_font_size=10)
    fig.write_image(str(FIGURES_DIR / "molecule_similarity_distribution.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'molecule_similarity_distribution.png'}")

    # Figure 2: Top hits bar chart - purple gradient by similarity
    top_30 = results_df.head(30)

    # Purple gradient: high similarity = dark purple, low = light
    def sim_to_purple(sim):
        # Map 0-1 to light purple (#c5b0d5) to dark purple (#9467bd)
        if sim >= 0.5:
            return "#9467bd"  # High similarity - dark purple
        elif sim >= 0.3:
            return "#c5b0d5"  # Moderate - light purple
        else:
            return "#7f7f7f"  # Low - gray

    colors = [sim_to_purple(s) for s in top_30["similarity"]]

    fig = go.Figure(go.Bar(
        y=top_30["name"][::-1],
        x=top_30["similarity"][::-1],
        orientation="h",
        marker_color=colors[::-1],
        text=[f"{s:.2f}" for s in top_30["similarity"][::-1]],
        textposition="outside",
        textfont_size=9,
    ))

    fig.update_layout(
        xaxis_title="Tanimoto Similarity",
        height=700,
        width=600,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(l=120, t=30, b=50),
    )
    fig.update_xaxes(showgrid=False, range=[0, 1.0], title_font_size=10)
    fig.update_yaxes(showgrid=False, tickfont_size=8)
    fig.write_image(str(FIGURES_DIR / "molecule_top_hits.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'molecule_top_hits.png'}")

    # Figure 3: Similarity vs MW scatter - single color with query marked
    fig = go.Figure()

    # All hits in purple
    fig.add_trace(go.Scatter(
        x=results_df["similarity"],
        y=results_df["mw"],
        mode="markers",
        name="Hits",
        marker=dict(size=7, color="#9467bd", opacity=0.6),
        text=results_df["name"],
        hovertemplate="%{text}<br>Sim: %{x:.3f}<br>MW: %{y:.1f}<extra></extra>",
    ))

    # Mark query compound
    fig.add_trace(go.Scatter(
        x=[1.0],
        y=[query_props["mw"]],
        mode="markers",
        name="Query",
        marker=dict(size=12, color="black", symbol="star"),
    ))

    fig.update_layout(
        xaxis_title="Tanimoto Similarity",
        yaxis_title="MW (Da)",
        height=450,
        width=550,
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=dict(x=0.02, y=0.98, font_size=9),
        margin=dict(t=30, b=50),
    )
    fig.update_xaxes(showgrid=False, range=[0, 1.05], title_font_size=10)
    fig.update_yaxes(showgrid=False, title_font_size=10)
    fig.write_image(str(FIGURES_DIR / "molecule_similarity_mw.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'molecule_similarity_mw.png'}")

    # Figure 4: 2D molecule structures (query + top 3 hits)
    print("\n  Creating 2D molecule visualization...")
    from rdkit.Chem import Draw

    # Prepare molecules for visualization
    top_3 = results_df.head(3)
    mols_to_draw = []
    legends = []

    # Query compound first
    query_mol = Chem.MolFromSmiles(QUERY_COMPOUND["smiles"])
    if query_mol:
        mols_to_draw.append(query_mol)
        legends.append(f"QUERY: {QUERY_COMPOUND['name']}")

    # Top 3 hits
    for i, (_, row) in enumerate(top_3.iterrows(), 1):
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol:
            mols_to_draw.append(mol)
            legends.append(f"#{i}: {row['name']}\nSim={row['similarity']:.3f}")

    # Draw grid of molecules
    if mols_to_draw:
        img = Draw.MolsToGridImage(
            mols_to_draw,
            molsPerRow=2,
            subImgSize=(400, 400),
            legends=legends,
            legendFontSize=14,
        )
        img.save(str(FIGURES_DIR / "molecule_2d_structures.png"))
        print(f"  Saved: {FIGURES_DIR / 'molecule_2d_structures.png'}")

    # -------------------------------------------------------------------------
    # Step 7: Save summary
    # -------------------------------------------------------------------------
    high_sim = results_df[results_df["similarity"] >= 0.7]
    moderate_sim = results_df[(results_df["similarity"] >= 0.5) & (results_df["similarity"] < 0.7)]

    summary = {
        "query": QUERY_COMPOUND,
        "database_size": len(COMPOUND_DATABASE),
        "high_similarity_count": len(high_sim),
        "moderate_similarity_count": len(moderate_sim),
        "top_hit": {
            "name": results_df.iloc[0]["name"],
            "similarity": float(results_df.iloc[0]["similarity"]),
            "class": results_df.iloc[0]["class"],
        },
        "class_summary": class_stats.to_dict(),
    }

    with open(OUTPUT_DIR / "similarity_search_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  Query: {QUERY_COMPOUND['name']} (ADRB2 inverse agonist from 2RH1)")
    print(f"  Database: {len(COMPOUND_DATABASE)} compounds searched")
    print(f"  High similarity (>=0.7): {len(high_sim)} compounds")
    print(f"  Moderate similarity (0.5-0.7): {len(moderate_sim)} compounds")
    print(f"  Top hit: {results_df.iloc[0]['name']} (Sim={results_df.iloc[0]['similarity']:.3f})")
    print(f"\nOutputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

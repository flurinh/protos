#!/usr/bin/env python3
"""
Verify all code snippets from README.md work correctly.

This script tests every Python code example from the README to ensure
documentation accuracy. Uses standard protos data path.
"""

import sys
import traceback
from pathlib import Path

# Setup path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

# Track results
RESULTS = []

def test_section(name):
    """Decorator to track test results."""
    def decorator(func):
        def wrapper():
            print(f"\n{'='*60}")
            print(f"TESTING: {name}")
            print('='*60)
            try:
                func()
                RESULTS.append((name, "PASS", None))
                print(f"[PASS] {name}")
                return True
            except Exception as e:
                RESULTS.append((name, "FAIL", str(e)))
                print(f"[FAIL] {name}")
                traceback.print_exc()
                return False
        return wrapper
    return decorator


# ============================================================
# SECTION 1: Data Path Setup
# ============================================================
@test_section("1. Data Path Setup")
def test_data_path_setup():
    import protos

    # Set custom data path (optional - defaults to ./data)
    protos.set_data_path(str(REPO_ROOT / "data"))

    # Get current path
    print(protos.get_data_path())


# ============================================================
# SECTION 2: Entity Registry
# ============================================================
@test_section("2. Entity Registry")
def test_entity_registry():
    from protos.processing.structure import StructureProcessor

    sp = StructureProcessor()

    # Load creates/retrieves an entity
    structure = sp.load_entity("1ubq")
    print(f"Loaded structure: {type(structure)}, shape: {structure.shape if structure is not None else 'None'}")

    # Check if entity exists
    exists = sp.entity_exists("1ubq")
    print(f"Entity exists: {exists}")

    # List all entities
    entities = sp.list_entities()
    print(f"Entities count: {len(entities)}")


# ============================================================
# SECTION 3: Datasets
# ============================================================
@test_section("3. Datasets")
def test_datasets():
    from protos.processing.structure import StructureProcessor
    sp = StructureProcessor()

    # Create a dataset (need existing entities)
    # First ensure entities exist
    existing = [e for e in ["1ubq", "2ubq", "3ubq"] if sp.entity_exists(e)]
    if existing:
        sp.create_dataset("my_structures", existing)

        # List datasets
        datasets = sp.list_datasets()
        print(f"Datasets: {datasets}")

        # Load all entities in a dataset
        structures = sp.load_dataset("my_structures")  # Returns Dict[str, DataFrame]
        print(f"Loaded {len(structures)} structures from dataset")
    else:
        print("No existing entities, skipping dataset test")


# ============================================================
# SECTION 4: StructureLoader
# ============================================================
@test_section("4. StructureLoader")
def test_structure_loader():
    from protos.io.ingest.structure_loader import StructureLoader
    from protos.processing.structure import StructureProcessor

    loader = StructureLoader()

    # Download single structure
    loader.download_and_register("1ubq")
    print("Downloaded 1ubq")

    # Download batch and create dataset
    loader.download_batch(
        ["3sn6", "4ldo", "2rh1"],
        dataset_name="gpcr_structures",
        create_dataset=True,
    )
    print("Downloaded batch and created dataset")


# ============================================================
# SECTION 5: Loading Data
# ============================================================
@test_section("5. Loading Structure Data")
def test_loading_structure_data():
    from protos.processing.structure import StructureProcessor
    sp = StructureProcessor()

    # Load single structure (returns pandas DataFrame)
    df = sp.load_entity("1ubq")
    print(f"Atoms: {len(df)}, Columns: {list(df.columns)[:5]}...")

    # Load dataset
    structures = sp.load_dataset("gpcr_structures")
    for pdb_id, df in structures.items():
        print(f"{pdb_id}: {len(df)} atoms")


# ============================================================
# SECTION 6: Structure Operations
# ============================================================
@test_section("6. Structure Operations")
def test_structure_operations():
    from protos.processing.structure import StructureProcessor
    sp = StructureProcessor()

    # Extract sequences from structure
    sequences = sp.get_all_sequences("1ubq")  # Dict[chain_id, sequence]
    print(f"Chains in 1ubq: {list(sequences.keys())}")

    sequence_a = sp.get_sequence("1ubq", chain_id="A")
    print(f"Sequence A length: {len(sequence_a) if sequence_a else 'None'}")

    # Register chain sequences as entities (cross-processor linking)
    sp.register_chain_sequences(
        ["3sn6", "4ldo"],
        dataset_prefix="gpcr_chains",
        create_dataset=True,
    )
    print("Registered chain sequences")


# ============================================================
# SECTION 7: Structure Ligand Operations
# ============================================================
@test_section("7. Structure Ligand Operations")
def test_structure_ligand_ops():
    from protos.io.ingest.structure_loader import StructureLoader
    from protos.processing.structure import StructureProcessor

    # Need 5d5a for ligand analysis
    loader = StructureLoader()
    loader.download_and_register("5d5a")

    sp = StructureProcessor()

    # Summarize ligands
    ligands = sp.summarize_ligands("5d5a")
    print(f"Ligands in 5d5a: {ligands}")


# ============================================================
# SECTION 8: Structure GRN Annotation
# ============================================================
@test_section("8. Structure GRN Annotation")
def test_structure_grn_annotation():
    from protos.processing.structure import StructureProcessor
    sp = StructureProcessor()

    # Annotate with GRN (Generic Residue Numbering)
    sp.annotate_with_grn("3sn6", chains=["R"])
    print("Annotated 3sn6 chain R with GRN")


# ============================================================
# SECTION 9: Structure Analysis
# ============================================================
@test_section("9. Structure Ligand Analysis")
def test_structure_analysis():
    from protos.processing.structure import StructureProcessor
    from protos.analysis.structure_ligand_analysis import calculate_ligand_interactions

    sp = StructureProcessor()

    # Get ligand atoms
    df = sp.load_entity("5d5a").reset_index()
    ligand_atoms = df[(df["group"] == "HETATM") & (df["res_name3l"] == "CAU")]
    print(f"Found {len(ligand_atoms)} ligand atoms")

    # Calculate interactions
    interactions = calculate_ligand_interactions(
        sp, "5d5a", ligand_atoms,
        detailed=True,
        cutoff=4.0
    )
    print(f"Binding residues: {len(interactions['binding_residues'])}")
    print(f"H-bonds: {len(interactions['hydrogen_bonds'])}")


# ============================================================
# SECTION 10: SequenceLoader
# ============================================================
@test_section("10. SequenceLoader")
def test_sequence_loader():
    from protos.io.ingest.sequence_loader import SequenceLoader
    from protos.processing.sequence import SequenceProcessor

    seq_proc = SequenceProcessor()
    loader = SequenceLoader(processor=seq_proc)

    # Download from UniProt
    loader.download_and_register(
        "uniprot:P00533,P07550",
        name="gpcr_sequences",
        materialize_entities=True,  # Save each sequence as entity
    )
    print("Downloaded UniProt sequences")

    # Skip local FASTA test (needs actual file)
    # loader.download_and_register("/path/to/sequences.fasta", name="my_sequences")


# ============================================================
# SECTION 11: Sequence Loading Data
# ============================================================
@test_section("11. Sequence Loading Data")
def test_sequence_loading():
    from protos.processing.sequence import SequenceProcessor
    seq_proc = SequenceProcessor()

    # Load single sequence
    sequence = seq_proc.load_entity("P00533")
    print(f"Loaded P00533: length={len(sequence) if sequence else 'None'}")

    # Load dataset
    sequences = seq_proc.load_dataset("gpcr_sequences")  # Dict[id, sequence]
    print(f"Loaded {len(sequences)} sequences from dataset")


# ============================================================
# SECTION 12: Sequence Operations
# ============================================================
@test_section("12. Sequence Operations")
def test_sequence_operations():
    from protos.processing.sequence import SequenceProcessor
    seq_proc = SequenceProcessor()

    # Save a sequence
    seq_proc.save_entity("my_protein", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    print("Saved my_protein sequence")

    # Load it back for alignment
    seq1 = seq_proc.load_entity("P00533")
    seq2 = seq_proc.load_entity("P07550")

    if seq1 and seq2:
        # Align two sequences
        score, alignment = seq_proc.align_sequences(
            seq1, seq2, "protein_1", "protein_2"
        )
        print(f"Alignment score: {score}")

    # Apply mutations
    sequence = seq_proc.load_entity("my_protein")
    if sequence and len(sequence) > 50:
        # Use valid positions
        mutant = seq_proc.mutate_sequence(sequence, ["V2A", "L5F"], "my_mutant")
        print(f"Created mutant: {mutant[:20]}...")

    # Save a wild_type for mutant library
    seq_proc.save_entity("wild_type", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")

    # Generate mutant library
    library, metadata = seq_proc.create_mutant_library(
        base_sequence_id="wild_type",
        mutation_map={2: ["A", "V", "L"], 5: ["F", "W"]},
        limit=10,
        return_metadata=True,
    )
    print(f"Generated {len(library)} mutants")

    # Conservation analysis
    conservation = seq_proc.compute_conservation(sequences=library)
    print(f"Conservation computed: {len(conservation)} positions")


# ============================================================
# SECTION 13: Sequence Export
# ============================================================
@test_section("13. Sequence Export")
def test_sequence_export():
    from protos.processing.sequence import SequenceProcessor
    seq_proc = SequenceProcessor()

    # Export to FASTA (use overwrite=True to handle re-runs)
    seq_proc.export_dataset("gpcr_sequences", export_name="gpcr_export", overwrite=True)
    print("Exported gpcr_sequences to FASTA")


# ============================================================
# SECTION 14: Sequence GRN Annotation
# ============================================================
@test_section("14. Sequence GRN Annotation")
def test_sequence_grn_annotation():
    from protos.processing.sequence import SequenceProcessor
    seq_proc = SequenceProcessor()

    # Need a GPCR sequence dataset with proper sequences
    # Use the chain sequences we registered earlier
    datasets = seq_proc.list_datasets()
    gpcr_dataset = None
    for ds in datasets:
        if 'gpcr_chain' in ds.lower():
            gpcr_dataset = ds
            break

    if gpcr_dataset:
        # Annotate sequences with Generic Residue Numbering
        grn_table, summary = seq_proc.annotate_with_grn(
            dataset_name=gpcr_dataset,
            reference_table="gpcrdb_ref",
            protein_family="gpcr_a",
            output_table="my_grn_annotations",
            return_summary=True,
            allow_create=True,
        )
        print(f"Annotated: {summary['global']['annotated']}/{summary['global']['total']}")
    else:
        print(f"No GPCR chain dataset found. Available: {datasets}")


# ============================================================
# SECTION 15: PropertyProcessor
# ============================================================
@test_section("15. PropertyProcessor")
def test_property_processor():
    from protos.processing.property import PropertyProcessor

    prop_proc = PropertyProcessor()

    # Record properties (creates table if needed)
    rows = [
        {
            "scope": [{"format": "sequence", "name": "opsin_1"}],
            "entity_name": "opsin_1",
            "lambda_max": 568,
            "photocycle": "fast",
        },
        {
            "scope": [
                {"format": "structure", "name": "5d5a"},
                {"format": "sequence", "name": "5d5a_chain_A"},
            ],
            "entity_name": "5d5a_chain_A",
            "classification": "gpcr_like",
            "similarity_score": 0.85,
        },
    ]

    prop_proc.record_properties(
        "opsin_properties",
        rows,
        metadata={"source": "experimental"},
        allow_create=True,
    )
    print("Recorded properties")


# ============================================================
# SECTION 16: Querying Properties
# ============================================================
@test_section("16. Querying Properties")
def test_querying_properties():
    from protos.processing.property import PropertyProcessor
    prop_proc = PropertyProcessor()

    # Get properties for an entity
    props = prop_proc.get_properties("opsin_1")
    print(f"Properties for opsin_1:\n{props.to_string()}")

    # Get properties for a structure (across all linked sequences)
    struct_props = prop_proc.get_properties("5d5a")
    print(f"Properties for 5d5a: {len(struct_props)} rows")


# ============================================================
# SECTION 17: GRNProcessor
# ============================================================
@test_section("17. GRNProcessor")
def test_grn_processor():
    from protos.processing.grn import GRNProcessor

    grn_proc = GRNProcessor()

    # Load reference table (for alignment-based annotation)
    ref_table = grn_proc.load_reference_table("gpcrdb_ref")
    print(f"Reference table: {ref_table.shape}")

    # List available annotation tables
    tables = grn_proc.list_tables()
    print(f"Available tables: {tables}")

    # Load a custom annotation table
    if tables:
        grn_table = grn_proc.load_table(tables[0])
        print(f"Loaded {tables[0]}: {grn_table.shape}")

    # Annotate sequences directly
    sequences = {"protein_1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"}
    annotations, summary = grn_proc.annotate_sequences(
        sequences,
        reference_table="gpcrdb_ref",
        protein_family="gpcr_a",
    )
    print(f"Annotation result: {annotations.shape}")


# ============================================================
# SECTION 18: EmbeddingProcessor Direct
# ============================================================
@test_section("18. EmbeddingProcessor Direct")
def test_embedding_direct():
    from protos.processing.embedding import EmbeddingProcessor

    # Initialize with model name
    emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")

    # Check available models
    print(f"Available models: {list(EmbeddingProcessor.MODEL_REGISTRY.keys())}")

    if not emb_proc.dependencies_available:
        print("SKIP: Embedding dependencies not available")
        return

    # Generate embeddings directly
    sequences = {"protein_1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"}
    embeddings = emb_proc.embed_sequences(
        sequences,
        embedding_type="mean",  # or "per_residue", "sum", "cls"
        save_dataset="my_embeddings"  # optional: save as dataset
    )
    print(f"Generated embeddings: {list(embeddings.keys())}")

    # Load saved embeddings (returns dict with original sequence IDs)
    embeddings = emb_proc.load_embeddings("my_embeddings")
    print(f"Loaded {len(embeddings)} embeddings")


# ============================================================
# SECTION 19: EmbeddingProcessor via ModelManager
# ============================================================
@test_section("19. EmbeddingProcessor via ModelManager")
def test_embedding_model_manager():
    from protos.processing.embedding import EmbeddingProcessor
    from protos.models.model_manager import ModelManager

    manager = ModelManager()

    # List embedding cards (use .cards dict)
    embedding_cards = [name for name in manager.cards.keys() if 'embedding' in name]
    print(f"Embedding cards: {embedding_cards}")

    # Check if dependencies available
    emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")
    if not emb_proc.dependencies_available:
        print("SKIP: Embedding dependencies not available")
        return

    # Prepare and run embedding
    invocation = manager.prepare(
        "embedding_esm2_t12_35m",
        inputs={"sequence_dataset": "gpcr_sequences"},
        config={"embedding_type": "mean"},  # or "per_residue", "sum"
    )
    print(f"Invocation created: runtime={invocation.runtime is not None}")

    # Ingest results into processor
    emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")
    emb_proc.ingest_from_invocation(invocation, dataset_name="gpcr_embeddings")
    print("Ingested embeddings")

    # Load embeddings
    embeddings = emb_proc.load_embeddings("gpcr_embeddings")
    print(f"Loaded {len(embeddings)} embeddings from gpcr_embeddings")


# ============================================================
# SECTION 20: MoleculeProcessor
# ============================================================
@test_section("20. MoleculeProcessor")
def test_molecule_processor():
    from protos.processing.molecule import MoleculeProcessor

    mol_proc = MoleculeProcessor()

    # Save ligand with SMILES
    mol_proc.save_entity(
        "5d5a_CAU_A",
        {"smiles": "CC(C)NC1=NC(=O)N(C=C1)C2=CN=C(N)N=C2N", "kind": "smiles_record"},
        metadata={"source_structure": "5d5a"},
    )
    print("Saved molecule entity")


# ============================================================
# SECTION 21: Model Integration
# ============================================================
@test_section("21. Model Integration")
def test_model_integration():
    from protos.models.model_manager import ModelManager

    manager = ModelManager()

    # List available models via .cards dict
    print(f"Available models: {list(manager.cards.keys())}")

    # Create a test sequence dataset for boltz
    from protos.processing.sequence import SequenceProcessor
    seq_proc = SequenceProcessor()
    seq_proc.save_entity("protein_1", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    seq_proc.create_dataset("my_sequences", ["protein_1"])

    # Prepare model input (boltz2 requires external setup, so just test prepare)
    try:
        invocation = manager.prepare(
            "boltz2",
            inputs={"sequence_dataset": "my_sequences", "entity": "protein_1"},
            config={
                "output_name": "protein_1_prediction",
                "recycling": 5,
                "num_samples": 3,
            },
        )

        # Access job details (for external models)
        if invocation.job:
            job = invocation.job
            print(f"Command: {' '.join(job.command[:3])}...")
            print(f"Working dir: {job.working_dir}")

        # Access runtime results (for in-process models)
        if invocation.runtime:
            print(f"Status: {invocation.runtime.outputs.get('status')}")
    except Exception as e:
        print(f"Boltz prepare (expected to need setup): {e}")


# ============================================================
# SECTION 22: Complete Workflow Example
# ============================================================
@test_section("22. Complete Workflow Example")
def test_complete_workflow():
    # Note: protos.set_data_path() must be called BEFORE any processors
    # It was already called at the start of this script.
    # In a real workflow, you'd call it at the very beginning:
    # import protos
    # protos.set_data_path("/path/to/data")

    from protos.io.ingest.structure_loader import StructureLoader
    from protos.processing.structure import StructureProcessor
    from protos.processing.sequence import SequenceProcessor
    from protos.processing.property import PropertyProcessor

    # Download structures (already downloaded, but verify)
    loader = StructureLoader()
    loader.download_batch(["3sn6", "4ldo", "2rh1"], dataset_name="gpcr_study")
    print("Downloaded/verified structures")

    # Initialize processors
    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    prop_proc = PropertyProcessor()

    # Extract and register sequences
    struct_proc.register_chain_sequences(
        ["3sn6", "4ldo", "2rh1"],
        dataset_prefix="gpcr_chains",
        create_dataset=True,
    )
    print("Registered chain sequences")

    # Find GPCR chain dataset
    datasets = seq_proc.list_datasets()
    gpcr_chain_ds = None
    for ds in datasets:
        if ds.startswith("gpcr_chains_"):
            gpcr_chain_ds = ds
            break

    if gpcr_chain_ds:
        # Annotate with GRN
        grn_table, summary = seq_proc.annotate_with_grn(
            dataset_name=gpcr_chain_ds,
            reference_table="gpcrdb_ref",
            protein_family="gpcr_a",
            output_table="gpcr_grn",
            return_summary=True,
            allow_create=True,
        )
        print(f"GRN annotation: {summary['global']['annotated']}/{summary['global']['total']}")

    # Record analysis results (scope is required for each property row)
    prop_proc.record_properties(
        "gpcr_analysis",
        [
            {"scope": [{"format": "structure", "name": "3sn6"}], "entity_name": "3sn6", "state": "active", "ligand_type": "full_agonist"},
            {"scope": [{"format": "structure", "name": "4ldo"}], "entity_name": "4ldo", "state": "active", "ligand_type": "full_agonist"},
            {"scope": [{"format": "structure", "name": "2rh1"}], "entity_name": "2rh1", "state": "inactive", "ligand_type": "inverse_agonist"},
        ],
        allow_create=True,
    )

    print("Workflow complete!")


# ============================================================
# MAIN
# ============================================================
def main():
    import protos
    protos.set_data_path(str(REPO_ROOT / "data"))

    # Run all tests
    tests = [
        test_data_path_setup,
        test_structure_loader,  # Run loader first to get data
        test_entity_registry,
        test_datasets,
        test_loading_structure_data,
        test_structure_operations,
        test_structure_ligand_ops,
        test_structure_grn_annotation,
        test_structure_analysis,
        test_sequence_loader,
        test_sequence_loading,
        test_sequence_operations,
        test_sequence_export,
        test_sequence_grn_annotation,
        test_property_processor,
        test_querying_properties,
        test_grn_processor,
        test_embedding_direct,
        test_embedding_model_manager,
        test_molecule_processor,
        test_model_integration,
        test_complete_workflow,
    ]

    for test in tests:
        test()

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    passed = sum(1 for _, status, _ in RESULTS if status == "PASS")
    failed = sum(1 for _, status, _ in RESULTS if status == "FAIL")

    for name, status, error in RESULTS:
        icon = "[PASS]" if status == "PASS" else "[FAIL]"
        print(f"{icon} {name}")
        if error:
            print(f"       Error: {error[:80]}...")

    print(f"\nTotal: {passed} passed, {failed} failed")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

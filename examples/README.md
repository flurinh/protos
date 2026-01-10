# Protos Examples

This directory contains example scripts and test workflows organized into three categories.

## Directory Structure

```
examples/
  processor_tests/     # Single-processor focused scripts (12 working)
    embedding/         # EmbeddingProcessor examples
    grn/              # GRNProcessor examples
    sequence/         # SequenceProcessor examples
    structure/        # StructureProcessor examples

  workflow_tests/      # Multi-processor and model-based workflows (25 working)
    cross_processor/  # Cross-processor integration
    embedding/        # Embedding feature workflows
    graph/            # Graph-based workflows
    lambda/           # Lambda model Docker workflows
    ligand/           # Ligand dataset workflows
    models/           # ModelManager and model integration
    property/         # Property integration workflows
    sequence/         # Sequence analysis workflows
    structure/        # Structure analysis workflows
    visualization/    # Visualization scripts

  deprecated/          # Old/obsolete scripts (71 files)
```

## processor_tests/ (12 scripts)

Working scripts focused on individual processor functionality:

- **embedding/**: Sequence embedding with ESM2/Ankh models
- **grn/**: Generic Residue Numbering annotation and config
- **sequence/**: Sequence alignment, MMseqs2 integration
- **structure/**: Structure loading, alignment, water networks, format validation

## workflow_tests/ (25 scripts)

Working multi-step workflows combining processors or using models:

- **cross_processor/**: Integration tests across data types
- **embedding/**: End-to-end embedding workflows
- **graph/**: Graph construction workflows
- **lambda/**: Lambda spectral prediction via Docker
- **ligand/**: Ligand dataset processing
- **models/**: Boltz, BoltzGen, model manager workflows
- **property/**: Property dataset creation
- **sequence/**: Mutational studies, sequence workflows
- **structure/**: Structure annotation and analysis
- **visualization/**: GRN visualization for GPCRs

## deprecated/ (71 scripts)

Scripts using old APIs, hardcoded paths, or removed modules. Kept for reference.

## Running Examples

Most scripts can be run directly:

```bash
# Activate the protos environment first
conda activate protos

# Run a processor test
python examples/processor_tests/embedding/test_sequence_embedding.py

# Run a workflow
python examples/workflow_tests/lambda/test_docker_opsin.py --run --gpu
```

Many scripts support pytest for automated testing:

```bash
pytest examples/processor_tests/grn/test_grn_workflow.py -v
pytest examples/workflow_tests/lambda/ -m "not docker"  # Skip Docker tests
```

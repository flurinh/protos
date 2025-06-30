# Conversions

Protos provides seamless conversion capabilities between different biological data formats, enabling cross-format workflows and data integration. The framework handles format conversions automatically while preserving data integrity and relationships.

## Overview

Conversion types supported:
- **Structure → Sequence**: Extract sequences from 3D structures
- **Sequence → Structure**: Generate structures via AlphaFold
- **Sequence → Embeddings**: Generate ML representations
- **GRN ↔ Sequence**: Convert between numbering schemes
- **Format conversions**: PDB ↔ mmCIF, FASTA variants
- **Property mappings**: Cross-format property transfer

## Structure to Sequence Conversions

### Basic Sequence Extraction

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
cp = CifBaseProcessor()

# Extract sequence from single chain
structure = cp.load_structure("1ubq")
sequence = cp.extract_sequence("1ubq", chain="A")
# Returns: "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

# Extract all chains
sequences = cp.get_seq_dict()
# Returns: {"1ubq_A": "MQIFVKTLTG...", "2gb1_A": "MTYKLILNGK...", ...}

# Save extracted sequences
from protos.processing.sequence.seq_processor import SeqProcessor
sp = SeqProcessor()
sp.save_sequences(sequences, "structure_derived_sequences.fasta")
```

### Advanced Extraction Options

```python
# Extract with missing residues handled
sequence_complete = cp.extract_sequence(
    "1abc",
    chain="A",
    include_missing=True,  # Include missing residues as 'X'
    missing_char='X'
)

# Extract specific region
partial_sequence = cp.extract_sequence_region(
    "1abc",
    chain="A",
    start_residue=50,
    end_residue=150
)

# Extract with modifications
sequence_with_mods = cp.extract_sequence(
    "1abc",
    include_modifications=True,
    modification_notation="parentheses"  # e.g., "M(ox)KTAYIAK"
)
```

### Multi-Model Structures

```python
# Handle NMR structures with multiple models
sequences_by_model = {}
structure = cp.load_structure("2koc")  # NMR structure

for model in cp.get_models("2koc"):
    seq = cp.extract_sequence("2koc", model=model)
    sequences_by_model[f"2koc_model{model}"] = seq

# Check sequence consistency across models
unique_sequences = set(sequences_by_model.values())
if len(unique_sequences) > 1:
    print("Warning: Sequence variations across models")
```

## Sequence to Structure Conversions

### AlphaFold Structure Prediction

```python
from protos.loaders.alphafold_utils import download_alphafold_structures

# Predict structure from sequence
sequences = {
    "MY_PROTEIN": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
}

# Option 1: Download from AlphaFold DB (if available)
af_structures = download_alphafold_structures(
    uniprot_ids=["P12345"],  # If sequence has UniProt ID
    output_dir=cp.path_structure_dir
)

# Option 2: Use ColabFold/LocalColabFold
predicted_structures = cp.predict_structures(
    sequences,
    method="colabfold",
    num_models=5,
    use_templates=True
)

# Register predicted structures
for name, structure_path in predicted_structures.items():
    cp.register_structure(
        name=f"{name}_predicted",
        file_path=structure_path,
        metadata={
            "method": "AlphaFold2",
            "confidence": "pLDDT",
            "date": datetime.now().isoformat()
        }
    )
```

### Homology-Based Structure Generation

```python
# Find template structures
def find_structure_templates(sequence, identity_threshold=0.3):
    """Find potential templates for homology modeling."""
    
    # Search PDB for similar sequences
    blast_results = sp.blast_sequence(
        sequence,
        database="pdb",
        e_value=0.001
    )
    
    # Filter by identity
    templates = []
    for hit in blast_results:
        if hit['identity'] >= identity_threshold:
            templates.append({
                'pdb_id': hit['accession'].split('_')[0],
                'chain': hit['accession'].split('_')[1],
                'identity': hit['identity'],
                'coverage': hit['coverage']
            })
    
    return templates

# Download template structures
templates = find_structure_templates(sequence)
for template in templates[:5]:  # Top 5 templates
    cp.download_structure(template['pdb_id'])
```

## Sequence to Embedding Conversions

### Generate Embeddings

```python
from protos.processing.embedding.embedding_processor import EmbeddingProcessor
ep = EmbeddingProcessor()

# Single sequence to embedding
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
embedding = ep.embed_sequence(sequence, model="esm2_t33_650M_UR50D")

# Batch conversion
sequences = sp.load_sequences(["PROT1", "PROT2", "PROT3"])
embeddings = ep.embed_sequences(
    sequences,
    model="esm2_t33_650M_UR50D",
    embedding_type="mean",
    register_entities=True  # Auto-register in entity system
)

# Save embeddings
for name, embedding in embeddings.items():
    ep.save_embedding(name, embedding)
```

### Structure-Aware Embeddings

```python
# Combine sequence and structure information
def generate_structure_aware_embeddings(protein_id):
    """Generate embeddings that incorporate structural context."""
    
    # Get sequence
    sequence = sp.load_sequence(protein_id)
    
    # Get structure if available
    if cp.has_entity(protein_id):
        structure = cp.load_structure(protein_id)
        
        # Extract structural features
        features = {
            'secondary_structure': cp.get_secondary_structure(protein_id),
            'solvent_accessibility': cp.calculate_sasa(protein_id),
            'contact_map': cp.calculate_contact_map(protein_id)
        }
        
        # Generate structure-conditioned embedding
        embedding = ep.embed_with_structure(
            sequence=sequence,
            structural_features=features,
            model="structure_aware_bert"  # Hypothetical model
        )
    else:
        # Fallback to sequence-only
        embedding = ep.embed_sequence(sequence)
    
    return embedding
```

## GRN Conversions

### Sequence to GRN

```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
gp = GRNBaseProcessor()

# Convert sequence positions to GRN
sequence = "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRT"
grn_annotation = gp.sequence_to_grn(
    sequence=sequence,
    sequence_id="NEW_GPCR",
    family="gpcr_class_a",
    reference="RHO_HUMAN"
)

# Convert specific positions
positions_of_interest = [45, 83, 135, 296]  # Sequence positions
grn_positions = gp.convert_positions_to_grn(
    sequence_id="RHO_HUMAN",
    positions=positions_of_interest
)
# Returns: {"45": "1.45", "83": "2.50", "135": "3.50", "296": "7.43"}
```

### GRN to Sequence

```python
# Extract sequence from GRN table
grn_table = gp.load_grn_table("gpcr_alignment")
sequences = gp.get_seq_dict()

# Convert GRN position to sequence positions
grn_position = "3.50"  # DRY motif position
sequence_positions = gp.grn_to_sequence_positions(
    grn_position,
    sequences=sequences
)
# Returns: {"RHO_HUMAN": 135, "ADRB2_HUMAN": 134, "BACR_HALSA": 129}

# Extract residues at GRN position
residues_at_350 = gp.get_position_column("3.50")
# Returns: ["R135", "R134", "L129", ...]
```

## File Format Conversions

### Structure Format Conversions

```python
# PDB to mmCIF
cp.convert_structure_format(
    input_file="structure.pdb",
    output_file="structure.cif",
    input_format="pdb",
    output_format="mmcif"
)

# mmCIF to PDB (may lose some annotations)
cp.convert_structure_format(
    input_file="structure.cif",
    output_file="structure.pdb",
    input_format="mmcif",
    output_format="pdb",
    preserve_annotations=False  # PDB format limitations
)

# Batch conversion
pdb_files = Path("old_structures").glob("*.pdb")
for pdb_file in pdb_files:
    cif_file = pdb_file.with_suffix(".cif")
    cp.convert_structure_format(pdb_file, cif_file)
```

### Sequence Format Conversions

```python
# FASTA to different formats
alignment = sp.read_alignment("sequences.fasta", format="fasta")

# Convert to Clustal
sp.write_alignment(alignment, "sequences.aln", format="clustal")

# Convert to Stockholm (rich annotations)
sp.write_alignment(alignment, "sequences.sto", format="stockholm")

# Convert to PHYLIP (for phylogenetic analysis)
sp.write_alignment(alignment, "sequences.phy", format="phylip")

# GenBank to FASTA
sp.convert_sequence_format(
    input_file="sequences.gb",
    output_file="sequences.fasta",
    input_format="genbank",
    output_format="fasta",
    id_source="accession"  # or "locus", "gi"
)
```

### Property Format Conversions

```python
from protos.processing.property.property_processor import PropertyProcessor
pp = PropertyProcessor()

# CSV to Excel with multiple sheets
property_data = pp.load_property_dataset("screening_results")
pp.export_dataset(
    dataset_id="screening_results",
    output_file="results.xlsx",
    format="excel",
    sheets={
        "raw_data": property_data,
        "statistics": pp.get_dataset_statistics("screening_results"),
        "correlations": pp.calculate_property_correlations("screening_results")
    }
)

# JSON to CSV (flatten nested properties)
pp.convert_property_format(
    input_file="properties.json",
    output_file="properties.csv",
    flatten_nested=True,
    separator="_"  # For nested keys: "binding.affinity" → "binding_affinity"
)
```

## Cross-Format Property Transfer

### Structure to Sequence Properties

```python
def transfer_structure_properties_to_sequence(protein_id):
    """Transfer structural properties to sequence processor."""
    
    # Calculate structural properties
    structure = cp.load_structure(protein_id)
    
    struct_properties = {
        "structure_length": len(cp.extract_sequence(protein_id)),
        "num_chains": len(cp.get_chains(protein_id)),
        "resolution": cp.get_metadata(protein_id).get("resolution"),
        "secondary_structure_content": cp.calculate_ss_content(protein_id),
        "disorder_regions": cp.predict_disorder(protein_id)
    }
    
    # Transfer to sequence processor
    sp.add_properties(protein_id, struct_properties)
    
    # Also save to property processor
    pp.assign_properties(
        entity_name=protein_id,
        properties=struct_properties,
        dataset_id="structural_properties"
    )
    
    return struct_properties
```

### Embedding to Property Conversions

```python
def embedding_to_properties(embeddings, property_names=None):
    """Convert embeddings to interpretable properties."""
    
    # Use dimensionality reduction
    from sklearn.decomposition import PCA
    
    # Reduce to interpretable dimensions
    pca = PCA(n_components=10)
    reduced_embeddings = pca.fit_transform(list(embeddings.values()))
    
    # Create property dataset
    property_data = []
    for (name, _), reduced_emb in zip(embeddings.items(), reduced_embeddings):
        properties = {
            "protein_id": name,
            **{f"emb_PC{i+1}": val for i, val in enumerate(reduced_emb)}
        }
        
        # Add interpretable properties if model available
        if property_names:
            predicted_props = predict_properties_from_embedding(reduced_emb)
            properties.update(predicted_props)
        
        property_data.append(properties)
    
    # Save as property dataset
    pp.create_property_dataset(
        dataset_id="embedding_derived_properties",
        properties_df=pd.DataFrame(property_data),
        entity_column="protein_id"
    )
```

## Conversion Pipelines

### Complete Analysis Pipeline

```python
def complete_conversion_pipeline(uniprot_id):
    """Demonstrate all conversion types in analysis pipeline."""
    
    results = {"uniprot_id": uniprot_id}
    
    # 1. Download sequence
    sequences = sp.download_from_uniprot([uniprot_id])
    sequence = sequences[uniprot_id]
    results["sequence_length"] = len(sequence)
    
    # 2. Predict/download structure  
    try:
        # Try AlphaFold first
        af_structure = download_alphafold_structures([uniprot_id])
        structure_id = f"AF-{uniprot_id}-F1"
        results["structure_source"] = "AlphaFold"
    except:
        # Fallback to homology search
        templates = find_structure_templates(sequence)
        if templates:
            structure_id = templates[0]['pdb_id']
            cp.download_structure(structure_id)
            results["structure_source"] = "PDB template"
    
    # 3. Extract sequence from structure (validation)
    struct_sequence = cp.extract_sequence(structure_id)
    sequence_match = struct_sequence == sequence
    results["sequence_structure_match"] = sequence_match
    
    # 4. Generate embeddings
    embedding = ep.embed_sequence(sequence)
    results["embedding_dimension"] = embedding.shape
    
    # 5. GRN annotation (if applicable)
    if is_gpcr(sequence):  # Hypothetical function
        grn_annotation = gp.annotate_sequence(sequence, family="gpcr")
        results["grn_positions"] = len(grn_annotation)
    
    # 6. Calculate properties
    properties = {
        "mw": sp.calculate_molecular_weight(sequence),
        "pi": sp.calculate_isoelectric_point(sequence),
        "gravy": sp.calculate_gravy(sequence),
        "instability": sp.calculate_instability_index(sequence)
    }
    
    pp.create_property_dataset(
        dataset_id=f"{uniprot_id}_properties",
        properties_df=pd.DataFrame([properties]),
        entity_column=None
    )
    
    results["properties"] = properties
    
    return results
```

### Batch Conversion Workflow

```python
def batch_convert_formats(input_dir, output_dir, conversions):
    """Batch convert files between formats."""
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    results = []
    
    for conversion in conversions:
        input_pattern = conversion["input_pattern"]
        output_format = conversion["output_format"]
        processor = conversion["processor"]
        
        # Find matching files
        input_files = input_path.glob(input_pattern)
        
        for input_file in input_files:
            try:
                # Determine output file
                output_file = output_path / input_file.stem + f".{output_format}"
                
                # Perform conversion
                if processor == "structure":
                    cp.convert_structure_format(
                        input_file,
                        output_file,
                        output_format=output_format
                    )
                elif processor == "sequence":
                    sp.convert_sequence_format(
                        input_file,
                        output_file,
                        output_format=output_format
                    )
                
                results.append({
                    "input": str(input_file),
                    "output": str(output_file),
                    "status": "success"
                })
                
            except Exception as e:
                results.append({
                    "input": str(input_file),
                    "error": str(e),
                    "status": "failed"
                })
    
    return results

# Example usage
conversions = [
    {
        "input_pattern": "*.pdb",
        "output_format": "cif",
        "processor": "structure"
    },
    {
        "input_pattern": "*.gb",
        "output_format": "fasta",
        "processor": "sequence"
    }
]

results = batch_convert_formats("raw_data", "converted_data", conversions)
```

## Best Practices

### 1. Preserve Metadata

```python
# Always preserve metadata during conversions
def convert_with_metadata(input_file, output_file, input_format, output_format):
    """Convert file while preserving metadata."""
    
    # Read original metadata
    original_metadata = read_metadata(input_file, input_format)
    
    # Perform conversion
    convert_file(input_file, output_file, input_format, output_format)
    
    # Reattach metadata
    attach_metadata(output_file, original_metadata, output_format)
    
    # Verify conversion
    verify_conversion(input_file, output_file, check_metadata=True)
```

### 2. Validate Conversions

```python
# Always validate after conversion
def validate_structure_conversion(original_file, converted_file):
    """Validate structure file conversion."""
    
    # Load both structures
    original = cp.load_structure_raw(original_file)
    converted = cp.load_structure_raw(converted_file)
    
    # Check atom counts
    assert len(original) == len(converted), "Atom count mismatch"
    
    # Check sequences
    orig_seq = cp.extract_sequence_from_file(original_file)
    conv_seq = cp.extract_sequence_from_file(converted_file)
    assert orig_seq == conv_seq, "Sequence mismatch"
    
    # Check coordinates (allowing small numerical differences)
    coord_diff = np.abs(original[['x','y','z']] - converted[['x','y','z']])
    assert coord_diff.max().max() < 0.001, "Coordinate mismatch"
```

### 3. Handle Format Limitations

```python
# Be aware of format limitations
def convert_with_warnings(input_file, output_format):
    """Convert with format limitation warnings."""
    
    warnings = []
    
    if output_format == "pdb":
        # PDB format limitations
        if has_long_chain_ids(input_file):
            warnings.append("Chain IDs will be truncated to 1 character")
        
        if has_many_atoms(input_file, threshold=99999):
            warnings.append("Atom numbers will wrap after 99999")
        
        if has_long_residue_names(input_file):
            warnings.append("Residue names will be truncated to 3 characters")
    
    if warnings:
        print("Format conversion warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    
    # Proceed with conversion
    return convert_file(input_file, output_format)
```

### 4. Batch Processing

```python
# Efficient batch conversions
def parallel_batch_conversion(files, conversion_func, n_workers=4):
    """Perform batch conversions in parallel."""
    
    from concurrent.futures import ProcessPoolExecutor
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        
        for file in files:
            future = executor.submit(conversion_func, file)
            futures.append((file, future))
        
        results = []
        for file, future in futures:
            try:
                result = future.result()
                results.append({"file": file, "status": "success", "result": result})
            except Exception as e:
                results.append({"file": file, "status": "failed", "error": str(e)})
    
    return results
```

## Troubleshooting

### Common Issues

1. **Sequence mismatch after conversion**
```python
# Check for missing residues
missing = cp.find_missing_residues(structure_id)
if missing:
    print(f"Missing residues: {missing}")
    # Use include_missing=True in extraction
```

2. **Loss of information in format conversion**
```python
# Use format-specific fields
if converting_to_limited_format:
    # Save additional info separately
    save_supplementary_info(original_data, "supplementary.json")
```

3. **Memory issues with large conversions**
```python
# Use streaming for large files
def stream_convert_large_file(input_file, output_file):
    with open(input_file) as inp, open(output_file, 'w') as out:
        for record in parse_records_lazy(inp):
            converted = convert_record(record)
            write_record(out, converted)
```

## Summary

Protos conversion capabilities include:
- Seamless structure ↔ sequence conversions
- Sequence → embedding generation
- GRN ↔ sequence position mapping
- File format conversions
- Cross-format property transfer
- Validation and metadata preservation

The conversion system enables fluid workflows across different data types while maintaining data integrity and relationships.
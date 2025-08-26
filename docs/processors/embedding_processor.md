# EmbeddingProcessor

The EmbeddingProcessor generates and manages protein sequence embeddings using state-of-the-art language models, enabling machine learning applications and sequence analysis based on learned representations.

## Overview

EmbeddingProcessor provides:
- Integration with multiple embedding models (ESM, ProtTrans, Ankh)
- Efficient batch processing and caching
- Multiple embedding types (per-residue, mean, cls)
- Model comparison and ensemble methods
- GPU acceleration support
- Cross-format integration

## Basic Usage

### Initialization

```python
from protos.processing.embedding.embedding_processor import EmbeddingProcessor

# Create processor instance
ep = EmbeddingProcessor(name="embeddings")

# With specific model
ep = EmbeddingProcessor(
    name="esm_embeddings",
    default_model="esm2_t33_650M_UR50D",
    device="cuda"  # Use GPU if available
)
```

### Checking Dependencies

```python
# Check if embedding models are available
deps = ep.check_dependencies()
print(deps)
# {
#     'torch': True,
#     'transformers': True,
#     'fair-esm': True,
#     'bio-embeddings': False,
#     'ready': True
# }

# List available models
models = ep.list_available_models()
for model, info in models.items():
    print(f"{model}: {info['embedding_dim']}D, max_length={info['max_length']}")
```

### Generating Embeddings

```python
# Single sequence
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
embedding = ep.embed_sequence(sequence, model="esm2_t6_8M_UR50D")

# Multiple sequences
sequences = {
    "PROT1": "MKTAYIAKQRQISFVKSHFSRQLE",
    "PROT2": "MVLSPADKTNVKAAWGKVGAHAGE",
    "PROT3": "MGSSHHHHHHSSGLVPRGSHMASMT"
}
embeddings = ep.embed_sequences(sequences)

# With entity registration
embeddings = ep.embed_sequences(
    sequences,
    register_entities=True  # Automatically register in entity system
)
```

## Embedding Types

### Mean Embeddings

Average representation across all residues:

```python
# Generate mean embeddings (default)
mean_embeddings = ep.embed_sequences(
    sequences,
    embedding_type="mean"
)
# Shape: (n_sequences, embedding_dim)
```

### CLS Token Embeddings

Special classification token representation:

```python
# CLS token embeddings
cls_embeddings = ep.embed_sequences(
    sequences,
    embedding_type="cls"
)
# Shape: (n_sequences, embedding_dim)
```

### Per-Residue Embeddings

Full residue-level representations:

```python
# Per-residue embeddings
residue_embeddings = ep.embed_sequences(
    sequences,
    embedding_type="per_residue"
)
# Shape: (n_sequences, max_length, embedding_dim)
```

### Contact Predictions

Some models provide contact predictions:

```python
# Get contacts from ESM models
embeddings, contacts = ep.embed_with_contacts(
    sequence,
    model="esm2_t33_650M_UR50D"
)
# contacts shape: (length, length)
```

## Model Management

### Available Models

```python
# ESM models (Facebook AI)
esm_models = [
    "esm2_t6_8M_UR50D",      # 6 layers, 8M parameters
    "esm2_t12_35M_UR50D",    # 12 layers, 35M parameters
    "esm2_t30_150M_UR50D",   # 30 layers, 150M parameters
    "esm2_t33_650M_UR50D",   # 33 layers, 650M parameters
    "esm2_t36_3B_UR50D"      # 36 layers, 3B parameters
]

# ProtTrans models
prottrans_models = [
    "prot_t5_xl_uniref50",
    "prot_t5_xl_bfd",
    "prot_bert_bfd"
]

# Ankh models
ankh_models = [
    "ankh_base",
    "ankh_large"
]
```

### Model Selection

```python
# Choose model based on requirements
def select_model(sequence_length, quality_priority=True):
    """Select appropriate model based on constraints."""
    
    if sequence_length > 1000:
        # Use model with longer context
        return "esm2_t12_35M_UR50D"
    elif quality_priority:
        # Use larger model for quality
        return "esm2_t33_650M_UR50D"
    else:
        # Use smaller model for speed
        return "esm2_t6_8M_UR50D"

# Apply selection
model = select_model(len(sequence))
embedding = ep.embed_sequence(sequence, model=model)
```

### Model Comparison

```python
# Compare embeddings from different models
models_to_compare = [
    "esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D",
    "esm2_t30_150M_UR50D"
]

comparison_results = ep.compare_models(
    sequences,
    models=models_to_compare,
    metric="cosine_similarity"
)

# Ensemble embeddings
ensemble_embedding = ep.ensemble_embed(
    sequences,
    models=models_to_compare,
    method="mean"  # or "weighted", "concat"
)
```

## Batch Processing

### Efficient Batch Generation

```python
# Process large sequence sets
large_sequence_set = sp.load_sequences(sequence_list)  # 1000s of sequences

# Batch processing with progress
embeddings = ep.batch_embed_sequences(
    large_sequence_set,
    batch_size=32,
    show_progress=True,
    save_intermediate=True  # Save after each batch
)

# Resume interrupted batch
embeddings = ep.batch_embed_sequences(
    large_sequence_set,
    resume=True,  # Continue from last checkpoint
    checkpoint_dir="embeddings_checkpoint"
)
```

### Memory Management

```python
# Stream processing for very large datasets
def stream_embed(sequence_file, output_file):
    """Process sequences without loading all into memory."""
    
    with open(output_file, 'wb') as out:
        for seq_id, sequence in ep.read_sequences_lazy(sequence_file):
            # Process one at a time
            embedding = ep.embed_sequence(sequence)
            
            # Save immediately
            pickle.dump({seq_id: embedding}, out)
            
            # Clear cache
            ep.clear_cache()
```

## Analysis with Embeddings

### Similarity Analysis

```python
# Calculate embedding similarity
similarity_matrix = ep.calculate_similarity_matrix(
    embeddings,
    metric="cosine"  # or "euclidean", "manhattan"
)

# Find similar sequences
query_embedding = embeddings["QUERY_PROTEIN"]
similar_proteins = ep.find_similar_embeddings(
    query_embedding,
    embeddings,
    top_k=10
)

# Cluster by embeddings
clusters = ep.cluster_embeddings(
    embeddings,
    method="kmeans",
    n_clusters=5
)
```

### Dimensionality Reduction

```python
# Reduce dimensions for visualization
reduced_embeddings = ep.reduce_dimensions(
    embeddings,
    method="umap",  # or "tsne", "pca"
    n_components=2
)

# Visualize
ep.visualize_embeddings(
    reduced_embeddings,
    labels=sequence_labels,
    output="embedding_plot.png"
)
```

### Machine Learning Integration

```python
# Prepare for ML tasks
X, y = ep.prepare_ml_dataset(
    embeddings,
    labels=protein_labels
)

# Train classifier
from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier()
clf.fit(X, y)

# Feature importance from embeddings
feature_importance = ep.analyze_feature_importance(
    embeddings,
    labels,
    model=clf
)
```

## Advanced Features

### Fine-tuning Support

```python
# Load model for fine-tuning
model, tokenizer = ep.load_model_for_finetuning("esm2_t6_8M_UR50D")

# Prepare data for fine-tuning
training_data = ep.prepare_finetuning_data(
    sequences,
    labels,
    task="classification"
)

# Fine-tune (requires additional setup)
finetuned_model = ep.finetune_model(
    model,
    training_data,
    epochs=10,
    learning_rate=1e-5
)
```

### Custom Embedding Functions

```python
# Define custom embedding extraction
def custom_embedding_function(model_output, sequence_length):
    """Extract custom embedding from model output."""
    # Custom logic, e.g., weighted average
    weights = torch.linspace(0.5, 1.0, sequence_length)
    weighted_embedding = (model_output * weights.unsqueeze(-1)).mean(dim=1)
    return weighted_embedding

# Use custom function
embeddings = ep.embed_sequences(
    sequences,
    custom_fn=custom_embedding_function
)
```

### Caching and Persistence

```python
# Enable persistent caching
ep.enable_cache(cache_dir="embedding_cache")

# Load cached embeddings
cached = ep.load_cached_embeddings(["PROT1", "PROT2"])

# Save embeddings efficiently
ep.save_embeddings(
    embeddings,
    format="hdf5",  # or "pickle", "numpy"
    compress=True
)

# Load saved embeddings
loaded_embeddings = ep.load_embeddings("embeddings.h5")
```

## Integration with Other Processors

### Sequence Integration

```python
# Load sequences and generate embeddings
sp = SeqProcessor()
sequences = sp.load_sequences(["P12345", "Q9Y6K1", "O14786"])

# Generate embeddings
embeddings = ep.embed_sequences(sequences)

# Cluster sequences by embeddings
clusters = ep.cluster_embeddings(embeddings, n_clusters=3)

# Save clustered sequences
for cluster_id, proteins in clusters.items():
    cluster_seqs = {p: sequences[p] for p in proteins}
    sp.save_sequences(cluster_seqs, f"cluster_{cluster_id}.fasta")
```

### Structure Integration

```python
# Combine structure and sequence embeddings
cp = CifBaseProcessor()

# Extract sequences from structures
struct_sequences = cp.get_seq_dict()

# Generate embeddings
struct_embeddings = ep.embed_sequences(struct_sequences)

# Map embeddings to structure
for pdb_id, embedding in struct_embeddings.items():
    # Use embedding for structure analysis
    structure = cp.load_structure(pdb_id)
    # Annotate structure with embedding values...
```

### Property Prediction

```python
# Use embeddings for property prediction
pp = PropertyProcessor()

# Get proteins with known properties
training_proteins = pp.get_entities_with_property("thermal_stability")

# Generate embeddings
training_embeddings = ep.embed_sequences(training_proteins)

# Train property predictor
property_model = ep.train_property_predictor(
    training_embeddings,
    property_values,
    property_type="regression"
)

# Predict properties for new proteins
new_embeddings = ep.embed_sequences(test_proteins)
predictions = property_model.predict(new_embeddings)
```

## Best Practices

### 1. Model Selection

Choose models based on requirements:

```python
# For speed (small sequences)
fast_embedding = ep.embed_sequence(
    short_sequence,
    model="esm2_t6_8M_UR50D"
)

# For quality (important proteins)
quality_embedding = ep.embed_sequence(
    important_sequence,
    model="esm2_t33_650M_UR50D"
)

# For long sequences
long_seq_embedding = ep.embed_sequence(
    long_sequence,
    model="esm2_t12_35M_UR50D"  # Better length handling
)
```

### 2. Batch Size Optimization

```python
# Determine optimal batch size
optimal_batch_size = ep.determine_optimal_batch_size(
    sequence_lengths=[len(s) for s in sequences.values()],
    available_memory=8192  # MB
)

# Use optimal batch size
embeddings = ep.batch_embed_sequences(
    sequences,
    batch_size=optimal_batch_size
)
```

### 3. GPU Utilization

```python
# Check GPU availability
if ep.cuda_available():
    ep.set_device("cuda:0")
    
    # Enable mixed precision for faster processing
    ep.enable_mixed_precision()
    
    # Monitor GPU memory
    memory_usage = ep.get_gpu_memory_usage()
    print(f"GPU memory: {memory_usage}MB")
```

### 4. Error Handling

```python
# Handle sequence length limits
max_length = ep.get_model_max_length()

for seq_id, sequence in sequences.items():
    if len(sequence) > max_length:
        # Truncate or split
        chunks = ep.split_sequence(sequence, max_length)
        chunk_embeddings = []
        
        for chunk in chunks:
            emb = ep.embed_sequence(chunk)
            chunk_embeddings.append(emb)
        
        # Combine chunk embeddings
        final_embedding = ep.combine_chunk_embeddings(chunk_embeddings)
```

## Common Workflows

### Workflow 1: Functional Classification

```python
def classify_protein_function(sequences, known_functions):
    """Classify proteins by function using embeddings."""
    
    # Generate embeddings
    all_embeddings = ep.embed_sequences(sequences)
    
    # Separate training and test
    training_embeddings = {
        seq_id: emb 
        for seq_id, emb in all_embeddings.items() 
        if seq_id in known_functions
    }
    
    test_embeddings = {
        seq_id: emb 
        for seq_id, emb in all_embeddings.items() 
        if seq_id not in known_functions
    }
    
    # Train classifier
    X_train = np.array(list(training_embeddings.values()))
    y_train = [known_functions[sid] for sid in training_embeddings.keys()]
    
    from sklearn.svm import SVC
    classifier = SVC(kernel='rbf')
    classifier.fit(X_train, y_train)
    
    # Predict functions
    predictions = {}
    for seq_id, embedding in test_embeddings.items():
        pred_function = classifier.predict([embedding])[0]
        predictions[seq_id] = pred_function
    
    return predictions
```

### Workflow 2: Evolutionary Analysis

```python
def analyze_protein_evolution(sequences, species_info):
    """Analyze evolutionary relationships using embeddings."""
    
    # Generate embeddings
    embeddings = ep.embed_sequences(sequences)
    
    # Calculate distances
    distance_matrix = ep.calculate_distance_matrix(embeddings)
    
    # Hierarchical clustering
    from scipy.cluster.hierarchy import dendrogram, linkage
    linkage_matrix = linkage(distance_matrix, method='ward')
    
    # Group by species
    species_embeddings = {}
    for seq_id, embedding in embeddings.items():
        species = species_info[seq_id]
        if species not in species_embeddings:
            species_embeddings[species] = []
        species_embeddings[species].append(embedding)
    
    # Calculate species centroids
    species_centroids = {
        species: np.mean(embs, axis=0)
        for species, embs in species_embeddings.items()
    }
    
    return {
        "distance_matrix": distance_matrix,
        "linkage": linkage_matrix,
        "species_centroids": species_centroids
    }
```

### Workflow 3: Variant Effect Prediction

```python
def predict_variant_effects(wild_type_seq, variants):
    """Predict effects of sequence variants using embeddings."""
    
    # Embed wild-type
    wt_embedding = ep.embed_sequence(wild_type_seq)
    
    results = []
    for variant_name, variant_seq in variants.items():
        # Embed variant
        var_embedding = ep.embed_sequence(variant_seq)
        
        # Calculate changes
        embedding_change = var_embedding - wt_embedding
        distance = np.linalg.norm(embedding_change)
        
        # Get per-residue changes if needed
        wt_residue_emb = ep.embed_sequence(
            wild_type_seq,
            embedding_type="per_residue"
        )
        var_residue_emb = ep.embed_sequence(
            variant_seq,
            embedding_type="per_residue"
        )
        
        # Find most affected regions
        residue_changes = np.linalg.norm(
            var_residue_emb - wt_residue_emb,
            axis=1
        )
        
        results.append({
            "variant": variant_name,
            "embedding_distance": distance,
            "max_residue_change": np.max(residue_changes),
            "affected_region": np.argmax(residue_changes)
        })
    
    return pd.DataFrame(results)
```

## Troubleshooting

### Common Issues

1. **Memory errors with large batches**
```python
# Reduce batch size or use CPU
ep.set_device("cpu")
ep.batch_embed_sequences(sequences, batch_size=8)
```

2. **Slow embedding generation**
```python
# Use smaller model or enable caching
ep.enable_cache()
ep.set_model("esm2_t6_8M_UR50D")  # Faster model
```

3. **Model loading failures**
```python
# Clear cache and re-download
ep.clear_model_cache()
ep.download_model("esm2_t12_35M_UR50D", force=True)
```

## Summary

EmbeddingProcessor provides:
- State-of-the-art protein language models
- Multiple embedding types and extraction methods
- Efficient batch processing and caching
- Integration with ML frameworks
- Cross-processor compatibility
- Advanced analysis capabilities

The processor enables modern machine learning approaches for protein analysis, from functional classification to evolutionary studies.
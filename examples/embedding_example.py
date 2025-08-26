"""
Example demonstrating the new EmbeddingProcessor.

This script shows how to generate embeddings for protein sequences
using the BaseProcessor-compatible EmbeddingProcessor.
"""

import os
from pathlib import Path
from protos.processing.embedding import EmbeddingProcessor

# Example protein sequences
example_sequences = {
    "bacteriorhodopsin": "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD",
    "hemoglobin_alpha": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    "lysozyme": "KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV"
}


def main():
    """Run embedding generation example."""
    print("=== EmbeddingProcessor Example ===\n")
    
    # Initialize processor with a small model for testing
    print("1. Initializing EmbeddingProcessor...")
    processor = EmbeddingProcessor(
        name="example_embeddings",
        model_name="esm2_t6_8m",  # Small model for demo
        batch_size=2
    )
    
    # Check dependencies first
    deps = processor.check_dependencies()
    print(f"   Dependencies available: {deps}")
    
    if not deps['ready']:
        print("\n⚠️  Missing dependencies detected!")
        print("   This example requires PyTorch and Transformers to be installed.")
        print("   Install them with ONE of these commands:")
        print("     pip install -e '.[gpu]'       # GPU support (recommended)")
        print("     pip install -e '.[embedding]' # CPU-only support")
        print("     pip install torch transformers # Manual installation")
        print("\n   The EmbeddingProcessor can still be imported and configured,")
        print("   but embedding generation requires these dependencies.")
        return
    
    print(f"   Model: {processor.model_name}")
    print(f"   Device: {processor.device}")
    print(f"   Embedding dimension: {processor.get_embedding_dim()}\n")
    
    # List available models
    print("2. Available models:")
    for model in processor.list_available_models():
        print(f"   - {model}")
    print()
    
    # Generate embeddings
    print("3. Generating embeddings for example sequences...")
    
    # Mean embeddings (default)
    print("   a) Mean embeddings:")
    mean_embeddings = processor.embed_sequences(
        example_sequences,
        embedding_type="mean"
    )
    for seq_id, embedding in mean_embeddings.items():
        print(f"      {seq_id}: shape {embedding.shape}")
    
    # CLS embeddings
    print("\n   b) CLS token embeddings:")
    cls_embeddings = processor.embed_sequences(
        example_sequences,
        embedding_type="cls"
    )
    for seq_id, embedding in cls_embeddings.items():
        print(f"      {seq_id}: shape {embedding.shape}")
    
    # Per-residue embeddings
    print("\n   c) Per-residue embeddings:")
    residue_embeddings = processor.embed_sequences(
        example_sequences,
        embedding_type="per_residue"
    )
    for seq_id, embedding in residue_embeddings.items():
        seq_len = len(example_sequences[seq_id])
        print(f"      {seq_id}: shape {embedding.shape} (sequence length: {seq_len})")
    
    # Save embeddings as dataset
    print("\n4. Saving embeddings as dataset...")
    dataset_name = "example_protein_embeddings"
    processor.embed_sequences(
        example_sequences,
        embedding_type="mean",
        save_dataset=dataset_name
    )
    print(f"   Saved dataset: {dataset_name}")
    
    # List datasets
    print("\n5. Available datasets:")
    datasets = processor.list_datasets()
    for ds in datasets:
        info = processor.get_dataset_info(ds)
        print(f"   - {ds}: {info.get('num_sequences', 'N/A')} sequences")
    
    # Load saved embeddings
    print(f"\n6. Loading embeddings from dataset '{dataset_name}'...")
    loaded_embeddings = processor.load_embeddings(dataset_name)
    print(f"   Loaded {len(loaded_embeddings)} embeddings")
    
    # Single sequence example
    print("\n7. Single sequence embedding:")
    single_seq = "MTEYKLVVVGAGGVGKSALTVQFVQGIFVEYDPNIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDL"
    single_embedding = processor.embed_sequences(single_seq)
    print(f"   Input: {single_seq[:50]}...")
    print(f"   Embedding shape: {single_embedding.shape}")
    
    # Clear cache
    print("\n8. Clearing model cache...")
    processor.clear_cache()
    print("   Cache cleared")
    
    print("\n=== Example complete ===")


if __name__ == "__main__":
    main()
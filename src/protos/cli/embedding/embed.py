"""
CLI command for generating embeddings using the entity system.

This command now uses the modern EmbeddingProcessor with entity tracking.
"""

import argparse
import torch
import logging
from typing import List, Optional

from protos.processing.embedding.embedding_processor import EmbeddingProcessor
from protos.processing.sequence import SequenceProcessor
from protos.io.data_access import GlobalRegistry

logger = logging.getLogger(__name__)


def embed_sequences(
    name: str,
    max_length: int = 1000,
    model_name: str = 'esm2_t6_8M',
    batch_size: int = 32,
    overwrite: bool = False,
    sequence_ids: Optional[List[str]] = None
):
    """
    Generate embeddings for protein sequences using the entity system.
    
    Args:
        name: Name for the embedding processor instance
        max_length: Maximum sequence length to process
        model_name: Embedding model to use
        batch_size: Batch size for processing
        overwrite: Whether to overwrite existing embeddings
        sequence_ids: Specific sequence IDs to embed (if None, embed all registered sequences)
    """
    # Initialize processors
    seq_proc = SeqProcessor(name=f"{name}_seq")
    emb_proc = EmbeddingProcessor(name=name, model_name=model_name)
    
    # Get sequences to embed
    if sequence_ids:
        # Specific sequences requested
        sequences = {}
        for seq_id in sequence_ids:
            sequence = seq_proc.load_sequence_entity(seq_id)
            if sequence:
                sequences[seq_id] = sequence
            else:
                logger.warning(f"Sequence {seq_id} not found")
    else:
        # Get all registered sequences
        registry = GlobalRegistry()
        sequences = {}
        
        for entity_id, entity_data in registry.entity_registry.entities.items():
            if 'sequence' in entity_data.get('formats', {}):
                original_id = entity_data['original_id']
                sequence = seq_proc.load_sequence_entity(original_id)
                if sequence and len(sequence) <= max_length:
                    sequences[original_id] = sequence
    
    if not sequences:
        print("No sequences found to embed")
        return False
    
    print(f"Found {len(sequences)} sequences to embed (max length: {max_length})")
    
    # Check for existing embeddings
    if not overwrite:
        already_embedded = []
        for seq_id in list(sequences.keys()):
            if emb_proc.has_embedding(seq_id):
                already_embedded.append(seq_id)
                del sequences[seq_id]
        
        if already_embedded:
            print(f"Skipping {len(already_embedded)} already embedded sequences")
    
    if not sequences:
        print("All sequences already have embeddings")
        return True
    
    # Convert to format expected by embed_sequences
    sequence_list = list(sequences.values())
    id_list = list(sequences.keys())
    
    print(f"Embedding {len(sequences)} sequences with {model_name}...")
    
    # Generate embeddings
    try:
        embeddings = emb_proc.embed_sequences(
            sequences=sequence_list,
            sequence_ids=id_list,
            batch_size=batch_size,
            register_entities=True  # Automatically register in entity system
        )
        
        print(f"Successfully generated embeddings for {len(embeddings)} sequences")
        
        # Save the embedding dataset
        emb_proc.save_dataset(f"{name}_embeddings")
        print(f"Saved embeddings to dataset: {name}_embeddings")
        
        return True
        
    except Exception as e:
        logger.error(f"Error generating embeddings: {e}")
        print(f"Error: {e}")
        return False


def main():
    """Command-line entry point for the embedding tool."""
    parser = argparse.ArgumentParser(
        description='Generate embeddings for protein sequences using the entity system'
    )
    
    parser.add_argument(
        'name',
        help='Name for this embedding run'
    )
    
    parser.add_argument(
        '--sequences', '-s',
        nargs='+',
        help='Specific sequence IDs to embed (if not provided, embeds all registered sequences)'
    )
    
    parser.add_argument(
        '--max-length', '-l',
        type=int,
        default=1000,
        help='Maximum sequence length to process (default: 1000)'
    )
    
    parser.add_argument(
        '--model', '-m',
        default='esm2_t6_8M',
        choices=[
            'esm2_t6_8M',
            'esm2_t12_35M', 
            'esm2_t30_150M',
            'esm2_t33_650M',
            'ankh_base',
            'ankh_large'
        ],
        help='Embedding model to use (default: esm2_t6_8M)'
    )
    
    parser.add_argument(
        '--batch-size', '-b',
        type=int,
        default=32,
        help='Batch size for processing (default: 32)'
    )
    
    parser.add_argument(
        '--overwrite', '-o',
        action='store_true',
        help='Overwrite existing embeddings'
    )
    
    parser.add_argument(
        '--device', '-d',
        choices=['cpu', 'cuda', 'auto'],
        default='auto',
        help='Device to use (default: auto-detect)'
    )
    
    args = parser.parse_args()
    
    # Set device
    if args.device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    else:
        device = args.device
    
    print(f"Using device: {device}")
    
    # Generate embeddings
    success = embed_sequences(
        name=args.name,
        max_length=args.max_length,
        model_name=args.model,
        batch_size=args.batch_size,
        overwrite=args.overwrite,
        sequence_ids=args.sequences
    )
    
    if success:
        print("Embedding generation completed successfully!")
    else:
        print("Embedding generation failed")
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
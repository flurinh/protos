import argparse
import torch
from protos.embedding.emb_processor import EMBProcessor, EMBModelManager
from protos.io.fasta_utils import read_fasta


def embed_sequences(name, max_length=1000, emb_model='ankh_large', batch_size=15000, overwrite=False):
    """
    Generate embeddings for protein sequences in a FASTA file.
    
    Args:
        name (str): Name of the dataset (used to locate the FASTA file and name the output)
        max_length (int, optional): Maximum sequence length to process. Defaults to 1000.
        emb_model (str, optional): Embedding model to use ('ankh_large' or 'esm2'). Defaults to 'ankh_large'.
        batch_size (int, optional): Batch size for processing. Defaults to 15000.
        overwrite (bool, optional): Whether to overwrite existing embeddings. Defaults to False.
    """
    fasta_path = f'data/fasta/processed/{name}.fasta'
    seq_dict = read_fasta(fasta_path)

    # Filter sequences by maximum length
    seq_dict = {k: v for k, v in seq_dict.items() if len(v) < max_length}
    
    dataset_name = f"{name}_{emb_model}"

    # Initialize model manager and processor
    model_manager = EMBModelManager(emb_model)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Setting device to {device}")
    embp = EMBProcessor(model_manager, device=device)

    # Load existing dataset if available
    if dataset_name in embp.list_available_datasets():
        embp.load_dataset(dataset_name)
        print("Loaded existing dataset...")
        embedded_seq_ids = embp.get_embedded_seq_ids()  # Using the helper method
        # Exclude already embedded sequences
        seq_dict = {k: v for k, v in seq_dict.items() if k not in embedded_seq_ids}

    print("Starting embedding...")

    # Define batch processing parameters
    total_sequences = len(seq_dict)
    num_batches = (total_sequences + batch_size - 1) // batch_size  # Ceiling division

    # Convert sequence dictionary to a list for easy slicing
    seq_items = list(seq_dict.items())

    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, total_sequences)
        print(f"Processing batch {batch_num + 1}/{num_batches}: sequences {start_idx} to {end_idx - 1}")

        # Embed the current batch
        embp.emb_seq_dict(
            seq_dict=dict(seq_items[start_idx:end_idx]),
            overwrite=overwrite
        )

        # Save the dataset after each batch to prevent data loss
        embp.save_dataset(dataset_name)
        print(f"Saved embeddings for batch {batch_num + 1}/{num_batches}.")

    print("Finished embedding all sequences.")
    return True


def main():
    """Command-line entry point for the embedding tool."""
    parser = argparse.ArgumentParser(description='Generate embeddings for protein sequences')
    parser.add_argument('--name', type=str, required=True, help='Name of the dataset')
    parser.add_argument('--max_length', type=int, default=1000, help='Maximum sequence length to process')
    parser.add_argument('--emb_model', type=str, default='ankh_large', choices=['ankh_large', 'esm2'], 
                        help='Embedding model to use')
    parser.add_argument('--batch_size', type=int, default=15000, help='Batch size for processing')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing embeddings')
    
    args = parser.parse_args()
    
    embed_sequences(
        name=args.name,
        max_length=args.max_length,
        emb_model=args.emb_model,
        batch_size=args.batch_size,
        overwrite=args.overwrite
    )
    
    print("Done!")


if __name__ == '__main__':
    main()
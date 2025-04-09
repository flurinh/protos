import argparse
import pickle
import pandas as pd
import numpy as np
import os
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.property.property_processor import PropertyProcessor
from protos.embedding.emb_processor import EMBProcessor
from models.predictor import Predictor


def save_dict_to_pickle(data, file_path, mode='ab'):
    """
    Save a dictionary to a pickle file.
    
    Args:
        data (dict): Dictionary to save
        file_path (str): Path to save the pickle file
        mode (str, optional): File open mode. Defaults to 'ab'.
    """
    with open(file_path, mode) as f:
        pickle.dump(data, f)


def load_dict_from_pickle(file_path):
    """
    Load a dictionary from a pickle file.
    
    Args:
        file_path (str): Path to the pickle file
    
    Returns:
        dict: Loaded dictionary or empty dict if file doesn't exist
    """
    if os.path.exists(file_path):
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    return {}


def process_batch(predictor, grnp, embp, protein_family, batch_size, protein_ids):
    """
    Process a batch of proteins for prediction.
    
    Args:
        predictor (Predictor): Predictor instance
        grnp (GRNProcessor): GRN processor instance
        embp (EMBProcessor): Embedding processor instance
        protein_family (str): Protein family name
        batch_size (int): Batch size for processing
        protein_ids (list): List of protein IDs to process (optional)
    
    Returns:
        tuple: (predictions, attention, global_embedding)
    """
    predictions, attn, global_embedding = predictor.predict(
        grnp,
        embp,
        protein_family=protein_family,
        batch_size=batch_size
    )
    return predictions, attn, global_embedding


def predict(dataset, run_id, protein_family, batch_size=8, save_attention=True):
    """
    Run predictions on a dataset.
    
    Args:
        dataset (str): Name of the dataset
        run_id (str): Run ID for the model
        protein_family (str): Protein family name
        batch_size (int, optional): Batch size for processing. Defaults to 8.
        save_attention (bool, optional): Whether to save attention maps. Defaults to True.
    
    Returns:
        pd.DataFrame: DataFrame with predictions
    """
    # Initialize processors to get all protein IDs
    grnp = GRNProcessor(f'{dataset}', 'data/grn/datasets/')
    grnp_ids = set(grnp.data.index.tolist())
    embp = EMBProcessor()
    emb_model = 'ankh_large'
    embp.load_dataset(dataset + '_' + emb_model)
    embp_ids = set(embp.emb_dict.keys())
    protein_ids = list(grnp_ids.intersection(embp_ids))  # valid ones
    print(f"Found {len(protein_ids)} valid protein ids (with preprocssed data).")

    # Define stats for normalization
    stats = {
        'Lmax_11_cis': (344.0, 611.0),
        'lmax_11_cis': (403.0, 611.0),
        'protonated_11_cis': (0.0, 1.0),
        'lmax_at': (436.0, 622.0),
        'Lmax_at': (385.0, 622.0),
        'lmax': (436.0, 622.0),
        'protonated': (0.0, 0.0),
        'lmax_uv': (344.0, 398.0),
        'protonated_at': (0.0, 1.0)
    }

    # Create attention folder if it doesn't exist
    if save_attention:
        attention_folder = 'data/attention'
        os.makedirs(attention_folder, exist_ok=True)
        attn_file = os.path.join(attention_folder, f'{dataset}_{run_id}_attn.pkl')
        pooled_embeddings_file = os.path.join(attention_folder, f'{dataset}_{run_id}_global_embedding.pkl')

    # Initialize a Predictor
    predictor = Predictor(run_id, stats)
    all_predictions = []

    # Process data in batches
    for start_idx in range(0, len(protein_ids), batch_size):
        end_idx = min(start_idx + batch_size, len(protein_ids))
        print(f"Processing batch {start_idx//batch_size + 1}/{(len(protein_ids)-1)//batch_size + 1}: proteins {start_idx} to {end_idx-1}")

        batch_protein_ids = protein_ids[start_idx:end_idx]

        # Reload and filter data for the current batch
        grnp = GRNProcessor(f'{dataset}', 'data/grn/datasets/')
        grnp.filter_by_ids(batch_protein_ids)

        embp = EMBProcessor()
        embp.load_dataset(dataset + '_' + emb_model)
        embp.filter_by_ids(batch_protein_ids)
        
        # Get predictions for this batch
        predictions, attn, global_embedding = process_batch(
            predictor, grnp, embp, protein_family, batch_size//2, []
        )
        
        all_predictions.append(predictions)
        
        # Save attention and embeddings if requested
        if save_attention:
            save_dict_to_pickle(attn, attn_file)
            save_dict_to_pickle(global_embedding, pooled_embeddings_file)

    # Combine all predictions
    combined_predictions = pd.concat(all_predictions, ignore_index=True)
    
    # Load the PropertyProcessor to update with new predictions
    pp = PropertyProcessor()
    pp.load_dataset(dataset)
    
    # Update with new predictions
    pp.properties.set_index('protein_id', inplace=True)
    combined_predictions.set_index('protein_id', inplace=True)
    pp.properties.update(combined_predictions)
    pp.properties.reset_index(inplace=True)
    
    # Save updated properties
    pp.compile_metadata()
    pp.save_dataset(dataset)
    
    print(f"Prediction completed for {len(protein_ids)} proteins in dataset: {dataset}")
    if save_attention:
        print(f"Attention saved to: {attn_file}")
        print(f"Pooled embeddings saved to: {pooled_embeddings_file}")
    
    return combined_predictions


def main():
    """Command-line entry point for prediction."""
    parser = argparse.ArgumentParser(description='Process GRN data and make predictions.')
    parser.add_argument('--dataset', type=str, required=True, 
                      help='Dataset name')
    parser.add_argument('--run_id', type=str, required=True, 
                      help='Run ID')
    parser.add_argument('--protein_family', type=str, required=True, 
                      help='Protein family')
    parser.add_argument('--batch_size', type=int, default=8, 
                      help='Batch size for processing')
    parser.add_argument('--no-save-attention', action='store_true',
                      help='Disable saving attention maps (saves disk space)')
    
    args = parser.parse_args()
    
    predict(
        dataset=args.dataset,
        run_id=args.run_id,
        protein_family=args.protein_family,
        batch_size=args.batch_size,
        save_attention=not args.no_save_attention
    )
    
    print("Done! Predictions completed successfully.")


if __name__ == '__main__':
    main()
"""
Example demonstrating how to use EmbeddingProcessor with GRN mapping.

This shows the proper way to extract residue embeddings and map them to:
- Sequence positions
- GRN assignments
- Structure information
"""

import os
from pathlib import Path
import torch

# Set up paths
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

from protos.processing.embedding.embedding_processor import EmbeddingProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.structure.struct_base_processor import CifBaseProcessor


def main():
    # Initialize processors
    embedding_processor = EmbeddingProcessor(
        name="embedding_example",
        model_name="esm2_t6_8m",  # Use small model for demonstration
        device="cpu"
    )
    
    grn_processor = GRNBaseProcessor(name="grn_example")
    structure_processor = CifBaseProcessor(name="structure_example")
    
    # Example sequence (bacteriorhodopsin fragment)
    sequence = "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"[:50]
    protein_id = "1AP9_A"
    
    print(f"Processing sequence of length {len(sequence)}")
    print(f"Sequence: {sequence}\n")
    
    # Step 1: Generate per-residue embeddings
    print("1. Generating per-residue embeddings...")
    full_embeddings = embedding_processor.embed_sequences(sequence, embedding_type="per_residue")
    print(f"   Full embeddings shape: {full_embeddings.shape} (includes CLS and EOS tokens)")
    
    # Step 2: Extract residue-only embeddings
    residue_embeddings = embedding_processor.get_residue_embeddings(full_embeddings)
    print(f"   Residue embeddings shape: {residue_embeddings.shape} (one per residue)")
    
    # Step 3: Simulate GRN assignments (in practice, this would come from GRNBaseProcessor)
    # For demonstration, we'll create mock GRN assignments
    print("\n2. Creating residue mapping with GRN assignments...")
    
    # Mock GRN assignments (position -> GRN)
    # In practice: grn_assignments = grn_processor.assign_grns(sequence, protein_id)
    grn_assignments = {}
    for i in range(1, 11):  # First 10 residues
        grn_assignments[i] = f"1.{50 + i - 1}"  # 1.50, 1.51, etc.
    for i in range(11, min(len(sequence) + 1, 21)):  # Next 10 residues
        grn_assignments[i] = f"2.{50 + i - 11}"  # 2.50, 2.51, etc.
    
    # Step 4: Create comprehensive residue mapping
    residue_data = []
    for i, aa in enumerate(sequence):
        position = i + 1  # 1-indexed position
        
        residue_info = {
            'position': position,
            'amino_acid': aa,
            'embedding': residue_embeddings[i],  # 0-indexed in tensor
            'grn': grn_assignments.get(position, 'unassigned'),
            'embedding_dim': residue_embeddings[i].shape[0]
        }
        residue_data.append(residue_info)
    
    # Step 5: Demonstrate usage
    print(f"\n3. Residue mapping created for {len(residue_data)} residues:")
    print("\nFirst 5 residues:")
    for res in residue_data[:5]:
        print(f"   Position {res['position']:3d}: {res['amino_acid']} -> GRN {res['grn']:8s} | Embedding shape: ({res['embedding_dim']},)")
    
    # Step 6: Show how to access embeddings by GRN
    print("\n4. Creating GRN -> Embedding mapping...")
    grn_embedding_map = {}
    for res in residue_data:
        if res['grn'] != 'unassigned':
            grn_embedding_map[res['grn']] = {
                'embedding': res['embedding'],
                'amino_acid': res['amino_acid'],
                'position': res['position']
            }
    
    print(f"   Created mapping for {len(grn_embedding_map)} GRN positions")
    print("\n   Sample GRN lookups:")
    for grn in ['1.50', '1.55', '2.50']:
        if grn in grn_embedding_map:
            info = grn_embedding_map[grn]
            print(f"   GRN {grn}: {info['amino_acid']} at position {info['position']}")
    
    # Step 7: Demonstrate batch processing
    print("\n5. Batch processing multiple sequences...")
    sequences = {
        "protein1": sequence[:20],
        "protein2": "ACDEFGHIKLMNPQRSTVWY",
        "protein3": "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELR"[:20]
    }
    
    batch_embeddings = embedding_processor.embed_sequences(sequences, embedding_type="per_residue")
    
    print("\n   Batch results:")
    for seq_id, embeddings in batch_embeddings.items():
        residue_embs = embedding_processor.get_residue_embeddings(embeddings)
        print(f"   {seq_id}: {len(sequences[seq_id])} residues -> {residue_embs.shape[0]} embeddings")
    
    # Step 8: Show how to save embeddings for later use
    print("\n6. Saving embeddings for later use...")
    # Create a structured dataset
    embedding_dataset = {
        'sequences': sequences,
        'embeddings': {},
        'metadata': {
            'model': embedding_processor.model_name,
            'embedding_type': 'per_residue',
            'include_special_tokens': False
        }
    }
    
    for seq_id, full_embs in batch_embeddings.items():
        # Store only residue embeddings (no special tokens)
        embedding_dataset['embeddings'][seq_id] = embedding_processor.get_residue_embeddings(full_embs)
    
    # In practice, you would save this:
    # torch.save(embedding_dataset, 'embeddings_with_grn.pt')
    
    print("   Dataset prepared for saving")
    print(f"   Contains embeddings for {len(embedding_dataset['embeddings'])} sequences")
    
    # Clean up
    embedding_processor.clear_cache()
    print("\n7. Cleaned up - model cache cleared")


if __name__ == "__main__":
    # Check dependencies
    try:
        import torch
        import transformers
        print("Dependencies available - running example...")
        main()
    except ImportError as e:
        print("Missing dependencies. Install with:")
        print("  pip install -e '.[gpu]'  # For GPU support")
        print("  pip install -e '.[embedding]'  # For CPU only")
        print(f"\nError: {e}")
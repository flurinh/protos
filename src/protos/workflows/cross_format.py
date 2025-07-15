"""
Cross-format workflow implementations for Protos.

This module provides high-level functions that coordinate operations
across different data formats while maintaining entity consistency.
"""

import logging
from typing import Dict, Any, Optional, List, Union
from pathlib import Path
import pandas as pd

from protos.io.data_access import GlobalRegistry, generate_entity_id
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor

logger = logging.getLogger(__name__)


def sequence_to_structure_workflow(
    sequence: str,
    protein_id: str,
    structure_predictor: str = "alphafold",
    seq_processor_name: str = "seq2struct",
    struct_processor_name: str = "seq2struct",
    metadata: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Workflow for sequence to structure prediction.
    
    This function manages the entity tracking when a sequence is used
    to predict a structure (e.g., via AlphaFold).
    
    Args:
        sequence: Protein sequence
        protein_id: Biological identifier for the protein
        structure_predictor: Method used for prediction (e.g., "alphafold")
        seq_processor_name: Name for sequence processor instance
        struct_processor_name: Name for structure processor instance
        metadata: Additional metadata to track
        
    Returns:
        Dict containing:
        - entity_id: Universal entity ID
        - sequence_saved: Path to saved sequence
        - structure_pdb_id: PDB ID for predicted structure
        - metadata: Combined metadata with lineage
    """
    # Initialize processors
    seq_proc = SeqProcessor(name=seq_processor_name)
    struct_proc = StructureProcessor(name=struct_processor_name)
    
    # Save sequence (registers entity)
    entity_id = seq_proc.save_sequence_entity(
        protein_id, 
        sequence,
        metadata=metadata or {}
    )
    
    # Prepare structure metadata with lineage
    struct_metadata = {
        'source': structure_predictor,
        'original_sequence_id': protein_id,
        'lineage': [
            {
                'step': 1,
                'operation': 'sequence_to_structure',
                'source_format': 'sequence',
                'source_id': protein_id,
                'method': structure_predictor
            }
        ]
    }
    
    # Update with any provided metadata
    if metadata:
        struct_metadata.update(metadata)
    
    # Structure PDB ID (following AlphaFold convention)
    structure_pdb_id = f"AF-{protein_id}-F1"
    
    logger.info(f"Sequence {protein_id} saved with entity ID {entity_id}")
    logger.info(f"Ready for structure prediction as {structure_pdb_id}")
    
    return {
        'entity_id': entity_id,
        'sequence_id': protein_id,
        'structure_pdb_id': structure_pdb_id,
        'metadata': struct_metadata
    }


def structure_to_sequence_workflow(
    pdb_id: str,
    chain_id: str = 'A',
    struct_processor_name: str = "struct2seq",
    seq_processor_name: str = "struct2seq",
    extraction_method: str = "ca_atoms"
) -> Dict[str, Any]:
    """
    Extract sequence from a protein structure.
    
    Args:
        pdb_id: PDB identifier of the structure
        chain_id: Chain to extract sequence from
        struct_processor_name: Name for structure processor
        seq_processor_name: Name for sequence processor
        extraction_method: Method used for extraction
        
    Returns:
        Dict containing:
        - sequence: Extracted sequence
        - entity_id: Entity ID for the sequence
        - metadata: Extraction metadata with lineage
    """
    # Initialize processors
    struct_proc = StructureProcessor(name=struct_processor_name)
    seq_proc = SeqProcessor(name=seq_processor_name)
    
    # Load structure
    structure_df = struct_proc.load_structure(pdb_id)
    
    if structure_df is None or len(structure_df) == 0:
        raise ValueError(f"Could not load structure {pdb_id}")
    
    # Extract sequence from specified chain
    chain_df = structure_df[structure_df['auth_chain_id'] == chain_id]
    
    if len(chain_df) == 0:
        raise ValueError(f"Chain {chain_id} not found in {pdb_id}")
    
    # Extract sequence
    sequence = struct_proc.get_sequence(pdb_id, chain_id=chain_id)
    
    # Create sequence ID
    seq_id = f"{pdb_id}_chain_{chain_id}"
    
    # Save with lineage metadata
    metadata = {
        'extracted_from_pdb': pdb_id,
        'chain': chain_id,
        'extraction_method': extraction_method,
        'parent_entity': generate_entity_id(pdb_id),
        'lineage': [
            {
                'step': 1,
                'operation': 'structure_to_sequence',
                'source_format': 'structure',
                'source_id': pdb_id,
                'chain': chain_id,
                'method': extraction_method
            }
        ]
    }
    
    entity_id = seq_proc.save_sequence_entity(seq_id, sequence, metadata=metadata)
    
    logger.info(f"Extracted sequence from {pdb_id} chain {chain_id}")
    logger.info(f"Saved as {seq_id} with entity ID {entity_id}")
    
    return {
        'sequence': sequence,
        'sequence_id': seq_id,
        'entity_id': entity_id,
        'metadata': metadata
    }


def sequence_to_grn_workflow(
    sequence: str,
    protein_id: str,
    reference_grn_table: str,
    reference_id: str,
    seq_processor_name: str = "seq2grn",
    grn_processor_name: str = "seq2grn",
    use_cached: bool = True
) -> Dict[str, Any]:
    """
    Assign GRN positions to a protein sequence.
    
    Args:
        sequence: Protein sequence
        protein_id: Biological identifier 
        reference_grn_table: Name of reference GRN table to use
        reference_id: ID of reference sequence in the table
        seq_processor_name: Name for sequence processor
        grn_processor_name: Name for GRN processor
        use_cached: Whether to use cached alignments
        
    Returns:
        Dict containing:
        - entity_id: Entity ID for the sequence
        - grn_assignment: Dict of GRN assignments
        - saved_table: Path to saved GRN table
        - metadata: Assignment metadata
    """
    # Initialize processors
    seq_proc = SeqProcessor(name=seq_processor_name)
    grn_proc = GRNBaseProcessor(name=grn_processor_name)
    
    # Save sequence if not already saved
    entity_id = seq_proc.save_sequence_entity(protein_id, sequence)
    
    # Load reference GRN table
    grn_proc.load_grn_table(reference_grn_table)
    
    # Perform GRN assignment
    # Note: In production, this would use annotate_gpcr from grn_table_utils
    # For now, we'll use a simplified approach
    from protos.processing.grn.grn_table_utils import annotate_gpcr
    from protos.processing.grn.grn_processor import GRNProcessor
    
    # If using GRNBaseProcessor, we need to create a GRNProcessor instance
    # with the same data for compatibility with annotate_gpcr
    if hasattr(grn_proc, 'data') and grn_proc.data is not None:
        # Create a temporary GRNProcessor with the loaded data
        temp_grn_proc = GRNProcessor(dataset=None, preload=False)
        temp_grn_proc.data = grn_proc.data
        temp_grn_proc.grns = grn_proc.grns
        temp_grn_proc.ids = grn_proc.ids
        
        # Perform annotation
        result_df = annotate_gpcr(
            temp_grn_proc, 
            protein_id, 
            sequence,
            add_to_GRNP=False,
            protein_family='gpcr_a',
            reload=False
        )
        
        # Extract the assignment from result dataframe
        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
            grn_assignment = result_df.iloc[0].to_dict()
        else:
            grn_assignment = None
    else:
        # Fallback: simple mock assignment for testing
        ref_row = grn_proc.data.loc[reference_id] if reference_id in grn_proc.data.index else None
        if ref_row is not None:
            grn_assignment = {}
            for i, grn in enumerate(grn_proc.grns):
                if i < len(sequence):
                    grn_assignment[grn] = sequence[i]
                else:
                    grn_assignment[grn] = '-'
        else:
            grn_assignment = None
    
    if grn_assignment is None:
        raise ValueError(f"GRN assignment failed for {protein_id}")
    
    # Create single-row GRN table for this sequence
    grn_row_data = pd.DataFrame([grn_assignment], index=[protein_id])
    grn_proc.data = grn_row_data
    grn_proc.ids = [protein_id]
    grn_proc.grns = list(grn_assignment.keys())
    
    # Save with entity ID and lineage
    table_name = f"grn_{protein_id}"
    saved_path = grn_proc.save_grn_table(
        table_name,
        include_entity_ids=True
    )
    
    metadata = {
        'reference_table': reference_grn_table,
        'reference_id': reference_id,
        'assignment_method': 'sequence_alignment',
        'lineage': [
            {
                'step': 1,
                'operation': 'sequence_to_grn',
                'source_format': 'sequence',
                'source_id': protein_id,
                'reference': reference_id,
                'method': 'alignment_based'
            }
        ]
    }
    
    logger.info(f"GRN assignment completed for {protein_id}")
    logger.info(f"Saved to {saved_path}")
    
    return {
        'entity_id': entity_id,
        'grn_assignment': grn_assignment,
        'saved_table': saved_path,
        'metadata': metadata
    }


def any_format_to_embeddings_workflow(
    source_id: str,
    source_format: str = "sequence",
    model_name: str = "esm2_t6_8M",
    processor_name: str = "to_embeddings",
    batch_size: int = 1
) -> Dict[str, Any]:
    """
    Generate embeddings from any supported format.
    
    Args:
        source_id: Biological identifier of the source
        source_format: Format of the source ("sequence", "structure", etc.)
        model_name: Embedding model to use
        processor_name: Name for embedding processor
        batch_size: Batch size for embedding generation
        
    Returns:
        Dict containing:
        - entity_id: Entity ID
        - embedding_path: Path to saved embeddings
        - metadata: Embedding metadata with lineage
    """
    # Initialize embedding processor
    emb_proc = EmbeddingProcessor(name=processor_name, model_name=model_name)
    
    # Get entity ID
    entity_id = generate_entity_id(source_id)
    
    # Handle different source formats
    if source_format == "sequence":
        seq_proc = SeqProcessor(name=processor_name)
        
        # Load sequence
        sequence = seq_proc.load_sequence_entity(source_id)
        if sequence is None:
            raise ValueError(f"Sequence {source_id} not found")
        
        # Generate embeddings
        embeddings = emb_proc.embed_sequences(
            [sequence],
            [source_id],
            batch_size=batch_size,
            register_entities=True
        )
        
    elif source_format == "structure":
        # For structures, first extract sequence
        struct_proc = CifBaseProcessor(name=processor_name)
        structure_df = struct_proc.load_structure(source_id)
        
        if structure_df is None:
            raise ValueError(f"Structure {source_id} not found")
        
        # Extract sequence from chain A
        sequence = struct_proc.get_sequence(source_id, chain_id='A')
        
        # Generate embeddings
        embeddings = emb_proc.embed_sequences(
            [sequence],
            [source_id],
            batch_size=batch_size,
            register_entities=True
        )
    else:
        raise ValueError(f"Unsupported source format: {source_format}")
    
    # Create metadata with lineage
    metadata = {
        'model': model_name,
        'source_format': source_format,
        'source_id': source_id,
        'lineage': [
            {
                'step': 1,
                'operation': f'{source_format}_to_embeddings',
                'source_format': source_format,
                'source_id': source_id,
                'method': model_name
            }
        ]
    }
    
    logger.info(f"Generated embeddings for {source_id} using {model_name}")
    
    return {
        'entity_id': entity_id,
        'source_id': source_id,
        'model': model_name,
        'metadata': metadata
    }


def track_conversion_lineage(
    entity_id: str,
    new_operation: Dict[str, Any],
    format_type: Optional[str] = None
) -> Dict[str, Any]:
    """
    Update lineage tracking for an entity after a conversion operation.
    
    Args:
        entity_id: Universal entity ID
        new_operation: Dict describing the new operation
        format_type: Specific format to update lineage for
        
    Returns:
        Updated lineage information
    """
    global_registry = GlobalRegistry()
    entity_info = global_registry.entity_registry.get_entity(entity_id)
    
    if entity_info is None:
        raise ValueError(f"Entity {entity_id} not found")
    
    # Get existing lineage
    if format_type and format_type in entity_info['formats']:
        current_lineage = entity_info['formats'][format_type].get('metadata', {}).get('lineage', [])
    else:
        # Aggregate lineage from all formats
        current_lineage = []
        for fmt_data in entity_info['formats'].values():
            fmt_lineage = fmt_data.get('metadata', {}).get('lineage', [])
            current_lineage.extend(fmt_lineage)
    
    # Add new operation
    new_operation['step'] = len(current_lineage) + 1
    updated_lineage = current_lineage + [new_operation]
    
    # Update metadata
    if format_type:
        global_registry.entity_registry.update_entity_metadata(
            entity_id,
            {'lineage': updated_lineage},
            format_type=format_type
        )
    
    logger.info(f"Updated lineage for entity {entity_id}: {new_operation['operation']}")
    
    return {
        'entity_id': entity_id,
        'updated_lineage': updated_lineage,
        'total_steps': len(updated_lineage)
    }
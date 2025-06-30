"""Sequence processor for managing sequence data and alignments."""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
import pandas as pd
import logging

from protos.core.base_processor import BaseProcessor
from protos.io.fasta_utils import read_fasta, write_fasta
from protos.io.file_utils import save_json, load_json
from protos.io.data_access import generate_entity_id
from .seq_alignment import (
    init_aligner, align_blosum62, format_alignment,
    mmseqs2_align, mmseqs2_align2, msa_blosum62,
    get_best_alignment, calc_alignment_score_restricted_area
)
from .seq_mutation_utils import (
    parse_mutation_str, apply_mutations_to_seq,
    generate_mutation_combinations, generate_rn_site_mutations
)

logger = logging.getLogger(__name__)


class SeqProcessor(BaseProcessor):
    """
    Processor for handling sequence data, alignments, and mutations.
    
    This processor manages:
    - Sequence loading from FASTA files
    - Sequence alignments (BioPython and MMseqs2)
    - Multiple sequence alignments
    - Sequence mutations and variant generation
    - Sequence metadata and annotations
    """
    
    def __init__(self, name: str = "seq_processor", 
                 processor_data_dir: str = "sequence",
                 **kwargs):
        """
        Initialize the sequence processor.
        
        Args:
            name: Processor name
            processor_data_dir: Subdirectory for sequence data
            **kwargs: Additional arguments for BaseProcessor
        """
        super().__init__(
            name=name,
            processor_data_dir=processor_data_dir,
            **kwargs
        )
        
        # Set up data directories
        self.data_dirs = {
            'fasta': Path(self.data_path) / 'fasta',
            'alignments': Path(self.data_path) / 'alignments',
            'metadata': Path(self.data_path) / 'metadata'
        }
        
        # Ensure directories exist
        for dir_path in self.data_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize aligner
        self.aligner = None
        self._init_aligner()
        
        # Cache for loaded sequences
        self.sequences = {}
        self.alignments = {}
        self.metadata = {}
        
    def _init_aligner(self, open_gap_score: int = -10):
        """Initialize BioPython aligner with BLOSUM62 matrix."""
        try:
            self.aligner = init_aligner(open_gap_score)
            logger.info("Initialized BioPython aligner with BLOSUM62")
        except Exception as e:
            logger.warning(f"Failed to initialize aligner: {e}")
            self.aligner = None
    
    def load_sequences(self, fasta_file: str, dataset_name: Optional[str] = None, 
                      register_entities: bool = True) -> Dict[str, str]:
        """
        Load sequences from a FASTA file with entity support.
        
        Args:
            fasta_file: Path to FASTA file or name in fasta/ directory
            dataset_name: Optional name to store sequences under
            register_entities: Whether to register each sequence as an entity
            
        Returns:
            Dictionary of sequence_id -> sequence
        """
        # Resolve path
        if not os.path.isabs(fasta_file):
            # Check in fasta directory
            fasta_path = self.data_dirs['fasta'] / fasta_file
            if not fasta_path.exists() and not fasta_file.endswith('.fasta'):
                fasta_path = self.data_dirs['fasta'] / f"{fasta_file}.fasta"
        else:
            fasta_path = Path(fasta_file)
            
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
            
        # Load sequences
        sequences = read_fasta(str(fasta_path))
        
        # Register entities if requested
        if register_entities:
            try:
                from protos.io.data_access import GlobalRegistry
                global_registry = GlobalRegistry()
                
                # Determine if this FASTA file should be registered as a dataset
                is_multi_sequence = len(sequences) > 1
                fasta_filename = fasta_path.stem
                
                # Register individual sequence entities
                for seq_id, sequence in sequences.items():
                    # Generate entity ID from sequence ID
                    entity_id = generate_entity_id(seq_id)
                    
                    # For single-sequence files, don't include dataset
                    # For multi-sequence files, use the filename as dataset
                    entity_datasets = []
                    if is_multi_sequence:
                        entity_datasets.append(fasta_filename)
                    if dataset_name:
                        entity_datasets.append(dataset_name)
                    
                    # Register entity
                    global_registry.entity_registry.register_entity(
                        entity_id=entity_id,
                        entity_type="sequence",
                        original_id=seq_id,
                        file_path=str(fasta_path),
                        metadata={
                            "length": len(sequence),
                            "fasta_file": fasta_filename,
                            "is_multi_sequence_file": is_multi_sequence
                        },
                        datasets=entity_datasets
                    )
                
                # Register the FASTA file as a dataset if it contains multiple sequences
                if is_multi_sequence and self.dataset_manager:
                    try:
                        self.dataset_manager.register_dataset(
                            dataset_id=fasta_filename,
                            name=fasta_filename,
                            description=f"Multi-sequence FASTA file containing {len(sequences)} sequences",
                            metadata={
                                "type": "multi_sequence_fasta",
                                "sequence_count": len(sequences),
                                "sequence_ids": list(sequences.keys()),
                                "file_path": str(fasta_path),
                                "file_size": fasta_path.stat().st_size
                            }
                        )
                        logger.info(f"Registered FASTA file '{fasta_filename}' as dataset")
                    except:
                        pass
                
                logger.info(f"Registered {len(sequences)} sequence entities from {fasta_filename}")
            except Exception as e:
                logger.warning(f"Could not register sequence entities: {e}")
        
        # Store in cache
        if dataset_name:
            self.sequences[dataset_name] = sequences
        else:
            self.sequences.update(sequences)
            
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
        return sequences
    
    def save_sequences(self, sequences: Dict[str, str], output_file: str,
                      dataset_name: Optional[str] = None):
        """
        Save sequences to a FASTA file.
        
        Args:
            sequences: Dictionary of sequence_id -> sequence
            output_file: Output filename
            dataset_name: Optional dataset name for registration
        """
        # Resolve output path
        if not os.path.isabs(output_file):
            output_path = self.data_dirs['fasta'] / output_file
            if not output_file.endswith('.fasta'):
                output_path = self.data_dirs['fasta'] / f"{output_file}.fasta"
        else:
            output_path = Path(output_file)
            
        # Ensure directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save sequences
        write_fasta(sequences, str(output_path))
        
        # Register dataset if name provided
        if dataset_name:
            self.register_dataset(dataset_name, str(output_path), "fasta")
            
        logger.info(f"Saved {len(sequences)} sequences to {output_path}")
    
    def align_sequences(self, seq1: str, seq2: str, 
                       seq1_id: str = "seq1", seq2_id: str = "seq2",
                       store_alignment: bool = True) -> Tuple[float, List[str]]:
        """
        Align two sequences using BioPython.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            seq1_id: ID for first sequence
            seq2_id: ID for second sequence
            store_alignment: Whether to store in alignment cache
            
        Returns:
            Tuple of (alignment score, formatted alignment)
        """
        if not self.aligner:
            raise RuntimeError("Aligner not initialized")
            
        # Perform alignment
        alignment = align_blosum62(seq1, seq2, self.aligner)
        formatted = format_alignment(alignment)
        
        # Store if requested
        if store_alignment:
            key = f"{seq1_id}_vs_{seq2_id}"
            self.alignments[key] = {
                'seq1_id': seq1_id,
                'seq2_id': seq2_id,
                'score': alignment.score,
                'alignment': formatted
            }
            
        return alignment.score, formatted
    
    def find_best_match(self, query_seq: str, reference_seqs: Dict[str, str],
                       use_mmseqs: bool = True) -> Tuple[str, float, List[str]]:
        """
        Find best matching sequence from references.
        
        Args:
            query_seq: Query sequence
            reference_seqs: Dictionary of reference sequences
            use_mmseqs: Whether to use MMseqs2 for fast search
            
        Returns:
            Tuple of (best_match_id, score, alignment)
        """
        if use_mmseqs:
            try:
                # Use MMseqs2 for fast search
                hits = mmseqs2_align(query_seq, reference_seqs)
                if hits is not None and not hits.empty:
                    best_hit = hits.iloc[0]
                    best_id = best_hit['target_id']
                    
                    # Get detailed alignment with BioPython
                    score, alignment = self.align_sequences(
                        query_seq, reference_seqs[best_id],
                        "query", best_id, store_alignment=False
                    )
                    
                    return best_id, score, alignment
            except Exception as e:
                logger.warning(f"MMseqs2 search failed: {e}, falling back to BioPython")
        
        # Fall back to BioPython exhaustive search
        best_id = None
        best_score = float('-inf')
        best_alignment = None
        
        for ref_id, ref_seq in reference_seqs.items():
            try:
                score, alignment = self.align_sequences(
                    query_seq, ref_seq,
                    "query", ref_id, store_alignment=False
                )
                if score > best_score:
                    best_score = score
                    best_id = ref_id
                    best_alignment = alignment
            except Exception as e:
                logger.warning(f"Failed to align with {ref_id}: {e}")
                continue
                
        return best_id, best_score, best_alignment
    
    def multiple_sequence_alignment(self, sequences: Dict[str, str],
                                  reference_seqs: Optional[Dict[str, str]] = None,
                                  use_mmseqs: bool = True) -> Dict[str, Tuple[str, float, List[str]]]:
        """
        Perform multiple sequence alignment against references.
        
        Args:
            sequences: Query sequences
            reference_seqs: Reference sequences (if None, align all-vs-all)
            use_mmseqs: Whether to use MMseqs2
            
        Returns:
            Dictionary of query_id -> (best_ref_id, score, alignment)
        """
        if reference_seqs is None:
            # All-vs-all alignment
            reference_seqs = sequences
            
        if use_mmseqs and len(sequences) > 1 and len(reference_seqs) > 1:
            try:
                # Use MMseqs2 for batch search
                hits = mmseqs2_align2(sequences, reference_seqs)
                if hits is not None and not hits.empty:
                    # Process results
                    results = {}
                    for query_id in sequences:
                        query_hits = hits[hits['query_id'] == query_id]
                        if not query_hits.empty:
                            best_hit = query_hits.iloc[0]
                            ref_id = best_hit['target_id']
                            
                            # Get detailed alignment
                            score, alignment = self.align_sequences(
                                sequences[query_id], reference_seqs[ref_id],
                                query_id, ref_id, store_alignment=True
                            )
                            results[query_id] = (ref_id, score, alignment)
                    
                    return results
            except Exception as e:
                logger.warning(f"MMseqs2 batch search failed: {e}")
        
        # Fall back to BioPython
        results = msa_blosum62(sequences, reference_seqs, self.aligner)
        
        # Convert to expected format
        formatted_results = {}
        for query_id, (ref_id, score, alignment) in results.items():
            formatted_results[query_id] = (ref_id, score, alignment)
            
        return formatted_results
    
    def mutate_sequence(self, sequence: str, mutations: List[str],
                       sequence_id: str = "mutant") -> str:
        """
        Apply mutations to a sequence.
        
        Args:
            sequence: Original sequence
            mutations: List of mutations in format "A123G"
            sequence_id: ID for the mutated sequence
            
        Returns:
            Mutated sequence
        """
        mutated = apply_mutations_to_seq(sequence, mutations)
        
        # Store if ID provided
        if sequence_id:
            self.sequences[sequence_id] = mutated
            
        return mutated
    
    def generate_variants(self, sequence: str, positions: List[int],
                         possible_aas: List[List[str]],
                         base_id: str = "variant") -> Dict[str, str]:
        """
        Generate all possible sequence variants at specified positions.
        
        Args:
            sequence: Base sequence
            positions: List of positions to mutate (1-indexed)
            possible_aas: List of possible amino acids for each position
            base_id: Base ID for variant naming
            
        Returns:
            Dictionary of variant_id -> sequence
        """
        combinations = generate_mutation_combinations(positions, possible_aas, sequence)
        
        variants = {}
        for i, mutations in enumerate(combinations):
            variant_id = f"{base_id}_{i+1}"
            if mutations:  # Only if there are actual mutations
                variant_seq = self.mutate_sequence(sequence, mutations)
                variants[variant_id] = variant_seq
                
        return variants
    
    def save_alignment(self, alignment_data: Dict, output_file: str):
        """Save alignment data to file."""
        output_path = self.data_dirs['alignments'] / output_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        save_json(alignment_data, str(output_path))
        logger.info(f"Saved alignment to {output_path}")
    
    def load_alignment(self, alignment_file: str) -> Dict:
        """Load alignment data from file."""
        alignment_path = self.data_dirs['alignments'] / alignment_file
        if not alignment_path.exists():
            raise FileNotFoundError(f"Alignment file not found: {alignment_path}")
            
        return load_json(str(alignment_path))
    
    def get_sequence_metadata(self, sequence_ids: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Get metadata for sequences.
        
        Args:
            sequence_ids: List of sequence IDs (None for all)
            
        Returns:
            DataFrame with sequence metadata
        """
        if sequence_ids is None:
            sequence_ids = list(self.sequences.keys())
            
        metadata = []
        for seq_id in sequence_ids:
            if seq_id in self.sequences:
                seq = self.sequences[seq_id]
                metadata.append({
                    'sequence_id': seq_id,
                    'length': len(seq),
                    'molecular_weight': self._calculate_mw(seq),
                    'isoelectric_point': self._calculate_pi(seq),
                    'amino_acid_composition': self._get_aa_composition(seq)
                })
                
        return pd.DataFrame(metadata)
    
    def _calculate_mw(self, sequence: str) -> float:
        """Calculate molecular weight of sequence."""
        # Simplified MW calculation
        mw_dict = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
        }
        
        mw = sum(mw_dict.get(aa, 0) for aa in sequence)
        # Subtract water for peptide bonds
        mw -= 18.0 * (len(sequence) - 1)
        return round(mw, 1)
    
    def _calculate_pi(self, sequence: str) -> float:
        """Calculate isoelectric point of sequence."""
        # Simplified pI calculation - just a placeholder
        # In reality, this requires iterative pH calculation
        basic_count = sum(1 for aa in sequence if aa in 'RKH')
        acidic_count = sum(1 for aa in sequence if aa in 'DE')
        
        if acidic_count > basic_count:
            return 4.5  # Acidic
        elif basic_count > acidic_count:
            return 9.5  # Basic
        else:
            return 7.0  # Neutral
    
    def _get_aa_composition(self, sequence: str) -> Dict[str, float]:
        """Get amino acid composition of sequence."""
        composition = {}
        seq_len = len(sequence)
        
        for aa in 'ARNDCQEGHILKMFPSTWYV':
            count = sequence.count(aa)
            composition[aa] = round(count / seq_len * 100, 1) if seq_len > 0 else 0
            
        return composition
    
    # Entity support methods
    def load_sequence_entity(self, identifier: str) -> Optional[str]:
        """
        Load a single sequence entity.
        
        Args:
            identifier: Sequence identifier (name or entity hash)
            
        Returns:
            Sequence string or None if not found
        """
        # Resolve identifier
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            entity_id = global_registry.entity_registry.resolve_identifier(identifier, format_type="sequence")
            
            # Get original ID
            original_id = global_registry.entity_registry.get_original_id(entity_id)
            if not original_id:
                original_id = identifier
        except:
            # Fallback
            original_id = identifier
            entity_id = generate_entity_id(identifier)
        
        # Check cache first
        for dataset, sequences in self.sequences.items():
            if original_id in sequences:
                return sequences[original_id]
        
        # Try to load from entity registry
        try:
            global_registry = GlobalRegistry()
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info and 'sequence' in entity_info.get('formats', {}):
                file_path = entity_info['formats']['sequence'].get('file_path')
                if file_path and os.path.exists(file_path):
                    # Load the FASTA file
                    sequences = read_fasta(file_path)
                    if original_id in sequences:
                        return sequences[original_id]
        except:
            pass
        
        logger.warning(f"Sequence entity not found: {identifier}")
        return None
    
    def save_sequence_entity(self, seq_id: str, sequence: str, 
                           datasets: Optional[List[str]] = None,
                           metadata: Optional[Dict[str, any]] = None) -> str:
        """
        Save a single sequence as an entity.
        
        Args:
            seq_id: Sequence identifier
            sequence: Sequence string
            datasets: Optional dataset IDs to associate with
            metadata: Additional metadata
            
        Returns:
            Entity ID
        """
        # Save sequence to individual file
        output_file = f"{seq_id}.fasta"
        output_path = self.data_dirs['fasta'] / output_file
        
        # Save single sequence
        write_fasta({seq_id: sequence}, str(output_path))
        
        # Generate entity ID
        entity_id = generate_entity_id(seq_id)
        
        # Register entity
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            
            # Merge metadata
            entity_metadata = {
                "length": len(sequence),
                "file": output_file
            }
            if metadata:
                entity_metadata.update(metadata)
            
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="sequence",
                original_id=seq_id,
                file_path=str(output_path),
                metadata=entity_metadata,
                datasets=datasets or []
            )
            logger.info(f"Saved sequence entity {seq_id} -> {entity_id}")
        except Exception as e:
            logger.warning(f"Could not register sequence entity: {e}")
        
        return entity_id
    
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all sequence entities (individual sequences within FASTA files).
        
        Each sequence within a FASTA file is an entity, not the file itself.
        
        Args:
            dataset: Optional dataset (FASTA filename) to filter by
            
        Returns:
            List of sequence IDs (not hash IDs!)
        """
        if dataset:
            # If dataset specified, load sequences from that FASTA file
            fasta_path = self.data_dirs['fasta'] / f"{dataset}.fasta"
            if not fasta_path.exists():
                fasta_path = self.data_dirs['fasta'] / f"{dataset}.fa"
            
            if fasta_path.exists():
                sequences = self.load_sequences(str(fasta_path))
                return list(sequences.keys())
            return []
        
        # List all sequences from all FASTA files
        all_sequences = []
        sequence_dir = self.data_dirs['fasta']
        
        if sequence_dir.exists():
            for fasta_file in sequence_dir.glob("*.fasta"):
                try:
                    sequences = self.load_sequences(str(fasta_file), register_entities=False)
                    all_sequences.extend(sequences.keys())
                except:
                    pass
            
            for fasta_file in sequence_dir.glob("*.fa"):
                try:
                    sequences = self.load_sequences(str(fasta_file), register_entities=False)
                    all_sequences.extend(sequences.keys())
                except:
                    pass
        
        # Also check registry for any registered sequences
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            entity_ids = global_registry.entity_registry.list_entities(format_type="sequence")
            
            for entity_id in entity_ids:
                original_id = global_registry.entity_registry.get_original_id(entity_id)
                if original_id and original_id not in all_sequences:
                    all_sequences.append(original_id)
        except:
            pass
            
        return list(set(all_sequences))  # Remove duplicates
    
    def list_sequence_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all sequence entities (backward compatibility).
        
        Deprecated: Use list_entities() instead.
        
        Args:
            dataset: Optional dataset to filter by
            
        Returns:
            List of sequence IDs
        """
        return self.list_entities(dataset=dataset)
    
    def get_sequence(self, identifier: str) -> Optional[str]:
        """
        Get a single sequence by identifier.
        
        This method allows accessing individual sequences that have been
        registered as entities, supporting the user's requirement that
        "we can access them individually".
        
        Args:
            identifier: Sequence identifier (ID or entity hash)
            
        Returns:
            Sequence string or None if not found
        """
        # Try load_sequence_entity first (handles entity resolution)
        sequence = self.load_sequence_entity(identifier)
        if sequence:
            return sequence
            
        # Fallback: check cached sequences
        if identifier in self.sequences:
            return self.sequences[identifier]
            
        # Try to find in all FASTA files
        for fasta_file in self.data_dirs['fasta'].glob("*.fasta"):
            try:
                sequences = read_fasta(str(fasta_file))
                if identifier in sequences:
                    return sequences[identifier]
            except:
                pass
                
        for fasta_file in self.data_dirs['fasta'].glob("*.fa"):
            try:
                sequences = read_fasta(str(fasta_file))
                if identifier in sequences:
                    return sequences[identifier]
            except:
                pass
        
        return None
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """
        List available sequence datasets.
        
        Datasets include:
        1. FASTA files with multiple sequences (>1 sequence)
        2. User-defined datasets from JSON files that group sequences
        
        FASTA files with only one sequence are NOT considered datasets.
        
        Returns:
            List of dataset information dictionaries
        """
        datasets = []
        sequence_dir = self.data_dirs['fasta']
        
        if sequence_dir.exists():
            # List all FASTA files with multiple sequences
            for fasta_file in sequence_dir.glob("*.fasta"):
                try:
                    # Quick scan to count sequences
                    with open(fasta_file, 'r') as f:
                        seq_count = sum(1 for line in f if line.startswith('>'))
                    
                    # Only consider files with multiple sequences as datasets
                    if seq_count > 1:
                        datasets.append({
                            'id': fasta_file.stem,
                            'type': 'multi_sequence_fasta',
                            'format': 'fasta',
                            'sequence_count': seq_count,
                            'file_path': str(fasta_file),
                            'file_size': fasta_file.stat().st_size,
                            'description': f'FASTA file containing {seq_count} sequences'
                        })
                except:
                    pass
            
            for fasta_file in sequence_dir.glob("*.fa"):
                try:
                    with open(fasta_file, 'r') as f:
                        seq_count = sum(1 for line in f if line.startswith('>'))
                    
                    # Only consider files with multiple sequences as datasets
                    if seq_count > 1:
                        datasets.append({
                            'id': fasta_file.stem,
                            'type': 'multi_sequence_fasta', 
                            'format': 'fasta',
                            'sequence_count': seq_count,
                            'file_path': str(fasta_file),
                            'file_size': fasta_file.stat().st_size,
                            'description': f'FASTA file containing {seq_count} sequences'
                        })
                except:
                    pass
        
        # Check for user-defined datasets in JSON files
        datasets_json = Path(self.data_path) / 'datasets.json'
        if datasets_json.exists():
            try:
                with open(datasets_json, 'r') as f:
                    user_datasets = json.load(f)
                    for dataset_id, dataset_info in user_datasets.items():
                        # User-defined datasets can group sequences from multiple sources
                        datasets.append({
                            'id': dataset_id,
                            'type': 'user_defined_dataset',
                            'description': dataset_info.get('description', ''),
                            'sequence_ids': dataset_info.get('sequence_ids', []),
                            'sequence_count': len(dataset_info.get('sequence_ids', [])),
                            'source': 'datasets.json'
                        })
            except:
                pass
        
        # Check dataset manager for managed datasets
        if self.dataset_manager:
            managed = self.dataset_manager.list_datasets()
            for ds in managed:
                if isinstance(ds, dict) and ds not in datasets:
                    datasets.append(ds)
        
        return datasets
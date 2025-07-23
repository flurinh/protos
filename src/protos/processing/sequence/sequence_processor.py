"""
Sequence processor for managing sequence data and alignments.

UPDATED: Follows DATA_MANAGEMENT_UNIFIED.md principles
- NO custom path handling - ProtosPaths manages ALL paths
- Zero configuration required
- Human-readable names for all operations
- Implements abstract methods from BaseProcessor
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
import pandas as pd
import logging

from protos.core.base_processor import BaseProcessor
from protos.io.fasta_utils import read_fasta, write_fasta
from protos.io.file_utils import save_json, load_json
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


class SequenceProcessor(BaseProcessor):
    """
    Processor for handling sequence data, alignments, and mutations.
    
    This processor manages:
    - Sequence loading from FASTA files
    - Sequence alignments (BioPython and MMseqs2)
    - Multiple sequence alignments
    - Sequence mutations and variant generation
    - Sequence metadata and annotations
    
    Key changes from legacy version:
    - NO path parameters in constructor
    - ALL paths managed by ProtosPaths
    - Implements abstract load_entity() and save_entity()
    - Automatic entity registration for sequences
    """
    
    def __init__(self, name: str = "seq_processor", paths=None):
        """
        Initialize the sequence processor.
        
        Args:
            name: Processor name
            paths: ProtosPaths instance (created if not provided)
        """
        # Initialize BaseProcessor with ProtosPaths
        super().__init__(name=name, paths=paths)
        
        # Initialize aligner
        self.aligner = None
        self._init_aligner()
        
        # Cache for loaded sequences
        self.sequences = {}
        self.alignments = {}
        self.metadata = {}
        
    # Path properties using ProtosPaths
    @property
    def path_fasta_dir(self):
        """Get path to FASTA files directory."""
        return self.get_subdirectory_path('fasta_dir')
        
    @property
    def path_alignments_dir(self):
        """Get path to alignments directory."""
        return self.get_subdirectory_path('alignment_dir')
        
    @property
    def path_metadata_dir(self):
        """Get path to metadata directory."""
        return self.get_subdirectory_path('metadata_dir')
        
    def _init_aligner(self, open_gap_score: int = -10):
        """Initialize BioPython aligner with BLOSUM62 matrix."""
        try:
            self.aligner = init_aligner(open_gap_score)
            logger.info("Initialized BioPython aligner with BLOSUM62")
        except Exception as e:
            logger.warning(f"Failed to initialize aligner: {e}")
            self.aligner = None
    
    def load_entity(self, name: str) -> Optional[Union[str, Dict[str, str]]]:
        """
        Load a sequence entity by name.
        
        This method implements the abstract load_entity from BaseProcessor.
        It can load either:
        1. A single sequence from a single-sequence FASTA file
        2. Multiple sequences from a multi-sequence FASTA file
        
        Args:
            name: Human-readable name (filename without extension or sequence ID)
            
        Returns:
            - For single-sequence files: The sequence string
            - For multi-sequence files: Dictionary of seq_id -> sequence
            - None if not found
        """
        # First, check if this is a FASTA file in the fasta directory
        fasta_path = self.path_fasta_dir / f"{name}.fasta"
        if not fasta_path.exists():
            fasta_path = self.path_fasta_dir / f"{name}.fa"
            
        if fasta_path.exists():
            sequences = read_fasta(str(fasta_path))
            if len(sequences) == 1:
                # Single sequence file - return just the sequence
                return list(sequences.values())[0]
            else:
                # Multi-sequence file - return dictionary
                return sequences
                
        # Check if it's a sequence ID within a multi-sequence file
        # This requires scanning through FASTA files
        for fasta_file in self.path_fasta_dir.glob("*.fasta"):
            try:
                sequences = read_fasta(str(fasta_file))
                if name in sequences:
                    return sequences[name]
            except:
                pass
                
        for fasta_file in self.path_fasta_dir.glob("*.fa"):
            try:
                sequences = read_fasta(str(fasta_file))
                if name in sequences:
                    return sequences[name]
            except:
                pass
                
        return None
        
    def save_entity(self, name: str, data: Union[str, Dict[str, str]]):
        """
        Save a sequence entity.
        
        This method implements the abstract save_entity from BaseProcessor.
        
        Args:
            name: Human-readable name for the entity
            data: Either a single sequence string or dict of seq_id -> sequence
        """
        # Sanitize filename
        safe_name = self._sanitize_filename(name)
        output_path = self.path_fasta_dir / f"{safe_name}.fasta"
        
        # Convert single sequence to dict format
        if isinstance(data, str):
            sequences = {name: data}
        else:
            sequences = data
            
        # Save sequences
        write_fasta(sequences, str(output_path))
        
        # Register the entity with the entity registry
        is_multi = len(sequences) > 1
        
        # Build metadata
        entity_metadata = {
            "file_path": str(output_path),
            "file_size": output_path.stat().st_size if output_path.exists() else 0
        }
        
        if is_multi:
            # Multi-sequence file metadata
            entity_metadata.update({
                "sequence_count": len(sequences),
                "sequence_ids": list(sequences.keys()),
                "type": "multi_sequence"
            })
        else:
            # Single sequence metadata
            seq_id = list(sequences.keys())[0]
            sequence = list(sequences.values())[0]
            entity_metadata.update({
                "sequence_id": seq_id,
                "length": len(sequence),
                "type": "single_sequence"
            })
        
        # Register with entity registry
        # Use relative path from data root
        rel_path = output_path.relative_to(self.paths.data_root)
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,  # 'sequence'
            file_path=str(rel_path),
            metadata=entity_metadata
        )
            
        logger.info(f"Saved {'dataset' if is_multi else 'sequence'} entity: {name}")
        
    def _sanitize_filename(self, name: str) -> str:
        """Sanitize a name for use as a filename."""
        # Replace problematic characters
        safe_name = name.replace('/', '_').replace('\\', '_').replace(':', '_')
        safe_name = safe_name.replace('|', '_').replace('*', '_').replace('?', '_')
        safe_name = safe_name.replace('<', '_').replace('>', '_').replace('"', '_')
        return safe_name
    
    def load_sequence(self, identifier: str, use_cache: bool = True) -> Optional[Union[str, Dict[str, str]]]:
        """
        Load a sequence or sequences from FASTA file.
        
        This is the primary method for loading sequences, similar to load_structure
        in CifBaseProcessor. It wraps load_entity with sequence-specific logic.
        
        Args:
            identifier: Sequence identifier (can be filename without extension or sequence ID)
            use_cache: Whether to use cached versions (for consistency with structure processor)
            
        Returns:
            - For single sequences: The sequence string
            - For multi-sequence files: Dictionary of seq_id -> sequence
            - None if not found
        """
        # Simply delegate to load_entity which handles all the logic
        return self.load_entity(identifier)
    
    def load_sequences(self, fasta_file: str, dataset_name: Optional[str] = None, 
                      register_entities: bool = True) -> Dict[str, str]:
        """
        Load sequences from a FASTA file with entity support.
        
        Args:
            fasta_file: Name of FASTA file (without path) in fasta/ directory
            dataset_name: Optional name to store sequences under
            register_entities: Whether to register each sequence as an entity
            
        Returns:
            Dictionary of sequence_id -> sequence
        """
        # Load using load_entity
        result = self.load_entity(fasta_file)
        
        if result is None:
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
            
        # Ensure we have a dictionary
        if isinstance(result, str):
            # Single sequence - use the filename as key
            sequences = {fasta_file: result}
        else:
            sequences = result
            
        # Register entities if requested
        if register_entities and len(sequences) > 1:
            # For multi-sequence files, register individual sequences
            for seq_id, sequence in sequences.items():
                # Each sequence in a multi-sequence file can be accessed individually
                # Create a relative path for this specific sequence
                rel_path = Path("sequence/fasta") / f"{fasta_file}"
                self.entity_registry.register_entity(
                    name=seq_id,
                    format_type=self.processor_type,
                    file_path=str(rel_path),
                    metadata={
                        "length": len(sequence),
                        "source_file": fasta_file,
                        "is_part_of_dataset": True,
                        "type": "sequence_in_dataset"
                    }
                )
            logger.info(f"Registered {len(sequences)} sequence entities from {fasta_file}")
        
        # Store in cache
        if dataset_name:
            self.sequences[dataset_name] = sequences
        else:
            self.sequences.update(sequences)
            
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_file}")
        return sequences
    
    def save_sequences(self, sequences: Dict[str, str], output_file: str,
                      dataset_name: Optional[str] = None):
        """
        Save sequences to a FASTA file.
        
        Args:
            sequences: Dictionary of sequence_id -> sequence
            output_file: Output filename (without path)
            dataset_name: Optional dataset name for registration
        """
        # Remove extension if provided
        if output_file.endswith('.fasta') or output_file.endswith('.fa'):
            name = output_file.rsplit('.', 1)[0]
        else:
            name = output_file
            
        # Save using save_entity
        self.save_entity(name, sequences)
        
        # Register dataset if name provided and multi-sequence
        if dataset_name and len(sequences) > 1:
            dataset_info = {
                "name": dataset_name,
                "entities": list(sequences.keys()),
                "metadata": {
                    "description": f"Sequence dataset with {len(sequences)} sequences",
                    "source_file": f"{name}.fasta"
                }
            }
            self.create_dataset(dataset_name, list(sequences.keys()), dataset_info["metadata"])
            
        logger.info(f"Saved {len(sequences)} sequences to {name}.fasta")
    
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
        output_path = self.path_alignments_dir / output_file
        # ProtosPaths handles directory creation
        
        save_json(alignment_data, str(output_path))
        logger.info(f"Saved alignment to {output_path}")
    
    def load_alignment(self, alignment_file: str) -> Dict:
        """Load alignment data from file."""
        alignment_path = self.path_alignments_dir / alignment_file
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
    
    def get_sequence(self, identifier: str) -> Optional[str]:
        """
        Get a single sequence by identifier.
        
        Args:
            identifier: Sequence identifier
            
        Returns:
            Sequence string or None if not found
        """
        # Try cache first
        if identifier in self.sequences:
            return self.sequences[identifier]
            
        # Try loading as entity
        result = self.load_entity(identifier)
        if isinstance(result, str):
            return result
        elif isinstance(result, dict) and identifier in result:
            return result[identifier]
            
        return None
    
    def list_entities(self, entity_type: Optional[str] = None) -> List[str]:
        """
        List all available sequence entities.
        
        Args:
            entity_type: Optional filter by type ('sequence' or 'sequence_dataset')
            
        Returns:
            List of entity names
        """
        entities = []
        
        # List FASTA files
        for fasta_file in self.path_fasta_dir.glob("*.fasta"):
            entities.append(fasta_file.stem)
            
        for fasta_file in self.path_fasta_dir.glob("*.fa"):
            entities.append(fasta_file.stem)
            
        # Also check entity registry
        registered = self.entity_registry.list_entities()
        entities.extend(registered)
        
        # Remove duplicates
        return list(set(entities))
    
    def list_datasets(self) -> List[str]:
        """
        List available sequence datasets (multi-sequence FASTA files).
        
        Returns:
            List of dataset names
        """
        datasets = []
        
        # Check for multi-sequence FASTA files
        for fasta_file in self.path_fasta_dir.glob("*.fasta"):
            try:
                # Quick scan to count sequences
                with open(fasta_file, 'r') as f:
                    seq_count = sum(1 for line in f if line.startswith('>'))
                
                # Only consider files with multiple sequences as datasets
                if seq_count > 1:
                    datasets.append(fasta_file.stem)
            except:
                pass
                
        for fasta_file in self.path_fasta_dir.glob("*.fa"):
            try:
                with open(fasta_file, 'r') as f:
                    seq_count = sum(1 for line in f if line.startswith('>'))
                
                if seq_count > 1:
                    datasets.append(fasta_file.stem)
            except:
                pass
                
        # Also check dataset manager
        if self.dataset_manager:
            managed_datasets = BaseProcessor.list_datasets(self)
            datasets.extend(managed_datasets)
            
        return list(set(datasets))
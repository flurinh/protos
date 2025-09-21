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
import math
import itertools
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any, Iterable
import pandas as pd
import numpy as np
import logging

from protos.io.core.base_processor import BaseProcessor
from protos.io.formats.fasta_utils import read_fasta, write_fasta
from protos.io.utils.file_utils import save_json, load_json
from protos.analysis.sequence import (
    init_biopython_aligner,
    perform_pairwise_alignment,
    format_alignment_human,
    run_mmseqs_alignment,
    MMseqsUnavailableError,
)
from protos.processing.sequence.seq_alignment import (
    mmseqs2_align,
    mmseqs2_align2,
    msa_blosum62,
    calc_alignment_score_restricted_area,
    get_best_alignment,
    get_score,
    check_chain_similarity,
)
from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine
from protos.processing.grn import GRNProcessor
from protos.processing.grn.grn_utils import get_seq, GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL
from .seq_mutation_utils import (
    parse_mutation_str, apply_mutations_to_seq,
    generate_mutation_combinations, generate_rn_site_mutations
)

logger = logging.getLogger(__name__)

STANDARD_AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")


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
    
    def __init__(self, name: str = "seq_processor"):
        """
        Initialize the sequence processor.
        
        Args:
            name: Processor name
        """
        # Initialize BaseProcessor with ProtosPaths
        super().__init__(name=name)
        
        # Initialize aligner
        self.aligner = init_biopython_aligner()
        self._alignment_engine = None
        
        # Cache for loaded sequences
        self.sequences = {}
        self.sequence_store: Dict[str, str] = {}
        self.alignments = {}
        self.metadata = {}
        
    # Path properties using ProtosPaths
    @property
    def path_entity_dir(self):
        """Get path to single-sequence entity directory."""
        return self.get_subdirectory_path('entity_fasta_dir')

    @property
    def path_fasta_dir(self):
        """Get path to dataset FASTA files directory."""
        return self.get_subdirectory_path('dataset_fasta_dir')

    @property
    def path_alignments_dir(self):
        """Get path to alignments directory."""
        return self.get_subdirectory_path('alignment_dir')
        
    @property
    def path_pairwise_alignments_dir(self):
        """Get path to pairwise alignments subdirectory."""
        return self.get_subdirectory_path('pairwise_alignment_dir')
        
    @property
    def path_multiple_alignments_dir(self):
        """Get path to multiple alignments subdirectory."""
        return self.get_subdirectory_path('multiple_alignment_dir')
        
    @property
    def path_mmseqs_alignments_dir(self):
        """Get path to MMseqs2 alignments subdirectory."""
        return self.get_subdirectory_path('mmseqs_alignment_dir')
        
    @property
    def path_databases_dir(self):
        """Get path to MMseqs2 databases directory."""
        return self.get_subdirectory_path('databases_dir')
        
    @property
    def path_metadata_dir(self):
        """Get path to metadata directory."""
        return self.get_subdirectory_path('metadata_dir')

    @property
    def alignment_engine(self):
        if self._alignment_engine is None:
            from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

            self._alignment_engine = SequenceAlignmentEngine()
        return self._alignment_engine
        
    def _init_aligner(self, open_gap_score: int = -10):
        """Deprecated helper; aligner is initialized via analysis utilities."""
        logger.warning("_init_aligner is deprecated; aligner initialized in __init__")
    
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
        # Look for a canonical entity FASTA file
        entity_path = Path(self.path_entity_dir) / f"{name}.fasta"
        if not entity_path.exists():
            entity_path = Path(self.path_entity_dir) / f"{name}.fa"

        if entity_path.exists():
            sequences = read_fasta(str(entity_path))
            # Entity files should contain exactly one sequence
            return list(sequences.values())[0] if len(sequences) == 1 else sequences

        return None
        
    def save_entity(
        self,
        name: str,
        data: Union[str, Dict[str, str]],
        metadata: Optional[Dict[str, Any]] = None,
    ):
        """
        Save a sequence entity.
        
        This method implements the abstract save_entity from BaseProcessor.
        
        Args:
            name: Human-readable name for the entity
            data: Either a single sequence string or dict of seq_id -> sequence
        """
        # Sanitize filename
        safe_name = self._sanitize_filename(name)
        output_path = Path(self.path_entity_dir) / f"{safe_name}.fasta"
        
        # Convert single sequence to dict format
        if isinstance(data, str):
            sequences = {name: data}
        else:
            if len(data) != 1:
                raise ValueError("save_entity expects a single sequence. Use save_sequences for multi-entry data.")
            sequences = data
            
        # Save sequences
        write_fasta(sequences, str(output_path))

        # Build metadata for single entity
        seq_id = list(sequences.keys())[0]
        sequence = list(sequences.values())[0]
        entity_metadata = {
            "file_path": str(output_path.relative_to(self.paths.data_root)),
            "length": len(sequence),
            "sequence_id": seq_id,
            "type": "single_sequence"
        }
        if metadata:
            entity_metadata.update(metadata)

        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=entity_metadata["file_path"],
            metadata=entity_metadata
        )

        logger.info(f"Saved sequence entity: {name}")
        
    def _sanitize_filename(self, name: str) -> str:
        """Sanitize a name for use as a filename."""
        # Replace problematic characters
        safe_name = name.replace('/', '_').replace('\\', '_').replace(':', '_')
        safe_name = safe_name.replace('|', '_').replace('*', '_').replace('?', '_')
        safe_name = safe_name.replace('<', '_').replace('>', '_').replace('"', '_')
        return safe_name
    
    def export_entity(
        self,
        name: str,
        out_path: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
        sequence_ids: Optional[List[str]] = None,
    ) -> Path:
        exporter = self._get_exporter()
        return exporter.export_entity(
            name,
            out_path,
            format=format,
            overwrite=overwrite,
            sequence_ids=sequence_ids,
        )

    def list_source_structures(
        self,
        sequence_ids: Optional[Union[str, Iterable[str]]] = None,
        *,
        rel_type: str = 'derived_from',
        direction: str = 'outgoing',
        include_unregistered: bool = False,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """List structure entities related to the provided sequences."""

        if sequence_ids is None:
            if include_unregistered:
                sequence_ids = self.list_entities()
            else:
                sequence_ids = list(self.sequence_store.keys()) or self.list_entities()

        return self.resolve_related_entities(
            sequence_ids,
            rel_type=rel_type,
            direction=direction,
            format_type='structure',
        )

    def list_dataset_source_structures(
        self,
        dataset_name: str,
        *,
        rel_type: str = 'derived_from',
        direction: str = 'outgoing',
    ) -> Dict[str, List[Dict[str, Any]]]:
        """List source structures for every sequence entity in a dataset."""

        return self.resolve_related_entities_for_dataset(
            dataset_name,
            rel_type=rel_type,
            direction=direction,
            format_type='structure',
        )

    # ---------- Internal Helpers ----------

    def _cache_sequences(
        self,
        sequences: Dict[str, str],
        *,
        dataset_name: Optional[str] = None,
    ) -> None:
        """Store sequences in local caches for quick reuse."""

        if dataset_name:
            existing = self.sequences.get(dataset_name, {})
            if isinstance(existing, dict):
                existing.update(sequences)
                self.sequences[dataset_name] = existing
            else:
                self.sequences[dataset_name] = sequences
        else:
            self.sequences.update(sequences)

        self.sequence_store.update(sequences)

    def _resolve_sequence(self, sequence_id: str) -> str:
        """Return sequence string for given identifier, loading if necessary."""

        if sequence_id in self.sequence_store:
            return self.sequence_store[sequence_id]

        # Check cached datasets
        for cache in self.sequences.values():
            if isinstance(cache, dict) and sequence_id in cache:
                seq = cache[sequence_id]
                self.sequence_store[sequence_id] = seq
                return seq

        data = self.load_entity(sequence_id)
        if isinstance(data, str):
            self.sequence_store[sequence_id] = data
            return data
        if isinstance(data, dict) and sequence_id in data:
            seq = data[sequence_id]
            self.sequence_store[sequence_id] = seq
            return seq

        raise ValueError(f"Sequence '{sequence_id}' not found in cache or registry")

    def _collect_sequences(
        self,
        sequences: Optional[Dict[str, str]] = None,
        dataset_name: Optional[str] = None,
    ) -> Dict[str, str]:
        """Normalize sequence inputs to a dictionary mapping id→sequence."""

        if sequences is not None:
            return OrderedDict(sequences)

        if dataset_name:
            loaded = self.load_dataset(dataset_name)
            self._cache_sequences(loaded, dataset_name=dataset_name)
            return OrderedDict(loaded)

        if self.sequence_store:
            return OrderedDict(sorted(self.sequence_store.items()))

        raise ValueError("No sequences available; load data or provide an explicit mapping")
    
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
        dataset_path = Path(self.path_fasta_dir) / fasta_file
        if not dataset_path.exists():
            # allow filename without extension
            dataset_path = Path(self.path_fasta_dir) / f"{fasta_file}.fasta"
        if not dataset_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        sequences = read_fasta(str(dataset_path))

        # Register entities if requested
        if register_entities:
            for seq_id, sequence in sequences.items():
                self.save_entity(seq_id, sequence)
            logger.info(f"Materialized {len(sequences)} sequence entities from {fasta_file}")

        dataset_key = dataset_name or self._sanitize_filename(dataset_path.stem)
        self._cache_sequences(sequences, dataset_name=dataset_key)

        dataset_metadata = {
            "artifact_path": str(dataset_path.relative_to(self.paths.data_root)),
            "sequence_count": len(sequences),
            "materialized": register_entities,
        }

        entity_names = list(sequences.keys()) if register_entities else []
        self.create_dataset(dataset_key, entity_names, dataset_metadata)

        logger.info(f"Loaded {len(sequences)} sequences from {dataset_path.name}")
        return sequences
    
    def save_sequences(
        self,
        sequences: Dict[str, str],
        output_file: str,
        dataset_name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        materialize_entities: bool = False,
    ) -> Path:
        """
        Save sequences to a FASTA file.

        Args:
            sequences: Dictionary of sequence_id -> sequence
            output_file: Output filename (without path)
            dataset_name: Optional dataset name for registration
            metadata: Optional metadata to attach to dataset record
            materialize_entities: Whether to save each sequence as an individual entity
        """
        if not sequences:
            raise ValueError("No sequences provided to save")

        if output_file.endswith(('.fasta', '.fa')):
            base_name = output_file.rsplit('.', 1)[0]
        else:
            base_name = output_file

        dataset_key = self._sanitize_filename(dataset_name or base_name)
        file_stem = self._sanitize_filename(base_name)

        dataset_path = Path(self.path_fasta_dir) / f"{file_stem}.fasta"
        write_fasta(sequences, str(dataset_path))

        if materialize_entities:
            for seq_id, sequence in sequences.items():
                self.save_entity(seq_id, sequence)

        dataset_metadata = {
            "artifact_path": str(dataset_path.relative_to(self.paths.data_root)),
            "sequence_count": len(sequences),
            "materialized": materialize_entities,
            "sequence_ids": list(sequences.keys()),
        }
        if metadata:
            dataset_metadata.update(metadata)

        entity_names = list(sequences.keys()) if materialize_entities else []
        self.create_dataset(dataset_key, entity_names, dataset_metadata)

        self._cache_sequences(sequences, dataset_name=dataset_key)

        logger.info("Saved %d sequences to %s", len(sequences), dataset_path)
        return dataset_path

    def load_dataset(self, dataset_name: str) -> Dict[str, str]:
        info = self.get_dataset_info(dataset_name)
        metadata = info.get('metadata', {}) if info else {}
        artifact_rel = metadata.get('artifact_path')

        sequences: Dict[str, str]
        if artifact_rel:
            dataset_path = Path(self.paths.data_root) / artifact_rel
            if not dataset_path.exists():
                raise FileNotFoundError(f"Dataset artifact missing at {dataset_path}")
            sequences = read_fasta(str(dataset_path))
        else:
            loaded = BaseProcessor.load_dataset(self, dataset_name)
            sequences = {}
            for key, value in loaded.items():
                if isinstance(value, str):
                    sequences[key] = value
                elif isinstance(value, dict):
                    sequences.update(value)

        self._cache_sequences(sequences, dataset_name=dataset_name)
        return sequences
    
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
        result = perform_pairwise_alignment(seq1, seq2, self.aligner, seq1_id, seq2_id)
        formatted = result.alignment_lines

        # Store if requested
        if store_alignment:
            key = f"{seq1_id}_vs_{seq2_id}"
            alignment_data = {
                'seq1_id': seq1_id,
                'seq2_id': seq2_id,
                'score': result.score,
                'alignment': formatted,
                'timestamp': pd.Timestamp.now().isoformat()
            }
            self.alignments[key] = alignment_data
            
            # Also save to file for persistence
            self.save_alignment(alignment_data, f"{key}.json", alignment_type="pairwise")
            
        return result.score, formatted

    def align_sequences_mmseqs(
        self,
        sequences: Dict[str, str],
    ) -> List[str]:
        try:
            return self.alignment_engine.align_pairwise_mmseqs(sequences)
        except MMseqsUnavailableError as exc:
            self.logger.warning("MMseqs unavailable: %s", exc)
            raise

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
            self.sequence_store[sequence_id] = mutated
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
                
        if variants:
            self._cache_sequences(variants)

        return variants

    def create_mutant_library(
        self,
        base_sequence_id: str,
        mutation_map: Dict[int, Iterable[str]],
        *,
        base_name: Optional[str] = None,
        include_wildtype: bool = True,
        limit: Optional[int] = None,
        zero_indexed: bool = False,
        register: bool = False,
        dataset_name: Optional[str] = None,
        materialize_entities: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
        return_metadata: bool = False,
    ) -> Union[Dict[str, str], Tuple[Dict[str, str], pd.DataFrame]]:
        """Generate a combinatorial mutant library for a base sequence.

        Args:
            base_sequence_id: Identifier of the sequence to mutate (must be registered).
            mutation_map: Mapping of residue positions (1-indexed by default) to
                iterable of replacement amino acids.
            base_name: Optional prefix for generated variant IDs.
            include_wildtype: Include the original sequence as ``<base_name>_WT``.
            limit: Maximum number of mutated variants to generate (``None`` = all).
            zero_indexed: Treat keys in ``mutation_map`` as zero-indexed positions.
            register: When True, persist the library via ``save_sequences``.
            dataset_name: Dataset name to use when registering (defaults to
                ``<base_name>_library``).
            materialize_entities: Save each variant as an entity when registering.
            metadata: Optional metadata dict to attach to the dataset when registering.
            return_metadata: If True, also return a DataFrame describing each variant.

        Returns:
            Either the variant mapping or a tuple of (variants, metadata_df).
        """

        if not mutation_map:
            raise ValueError("mutation_map must contain at least one position")

        base_sequence = self._resolve_sequence(base_sequence_id)
        base_length = len(base_sequence)

        if base_length == 0:
            raise ValueError(f"Base sequence '{base_sequence_id}' is empty")

        base_label = base_name or f"{self._sanitize_filename(base_sequence_id)}_mut"
        dataset_key = self._sanitize_filename(dataset_name or f"{base_label}_library")

        # Normalize positions
        positions = sorted(mutation_map.keys())
        if zero_indexed:
            positions = [pos + 1 for pos in positions]

        # Validate positions within sequence length
        for pos in positions:
            if pos < 1 or pos > base_length:
                raise ValueError(f"Mutation position {pos} out of range for sequence length {base_length}")

        residue_options: List[List[str]] = []
        for original_key in sorted(mutation_map.keys()):
            options = mutation_map[original_key]
            normalized: List[str] = []
            for aa in options:
                aa_upper = aa.upper()
                if aa_upper in {'X', '*'}:
                    normalized.extend(STANDARD_AMINO_ACIDS)
                else:
                    normalized.append(aa_upper)
            residue_options.append(sorted(set(normalized)))

        variants: "OrderedDict[str, str]" = OrderedDict()
        records: List[Dict[str, Any]] = []

        if include_wildtype:
            wt_id = f"{base_label}_WT"
            variants[wt_id] = base_sequence
            records.append({
                'variant_id': wt_id,
                'base_sequence_id': base_sequence_id,
                'mutations': [],
                'mutation_count': 0,
            })

        mutated_count = 0
        product_iter = itertools.product(*residue_options)

        for combo in product_iter:
            if limit is not None and mutated_count >= limit:
                break

            sequence_list = list(base_sequence)
            mutation_tags: List[str] = []

            for pos, aa in zip(positions, combo):
                idx = pos - 1
                original_aa = base_sequence[idx]
                if aa == original_aa:
                    continue
                sequence_list[idx] = aa
                mutation_tags.append(f"{original_aa}{pos}{aa}")

            if not mutation_tags:
                # Skip duplicates that are identical to the wildtype
                continue

            mutation_signature = '_'.join(mutation_tags)
            variant_id = f"{base_label}__{mutation_signature}"

            mutated_sequence = ''.join(sequence_list)
            variants[variant_id] = mutated_sequence
            records.append({
                'variant_id': variant_id,
                'base_sequence_id': base_sequence_id,
                'mutations': mutation_tags,
                'mutation_count': len(mutation_tags),
            })

            mutated_count += 1

        if not variants:
            raise ValueError("No variants generated; adjust mutation_map or parameters")

        self._cache_sequences(dict(variants), dataset_name=dataset_key)

        dataset_path: Optional[Path] = None
        if register:
            dataset_metadata = {
                'base_sequence_id': base_sequence_id,
                'positions': positions,
                'variant_count': len(variants),
            }
            if metadata:
                dataset_metadata.update(metadata)

            dataset_path = self.save_sequences(
                dict(variants),
                output_file=dataset_key,
                dataset_name=dataset_key,
                metadata=dataset_metadata,
                materialize_entities=materialize_entities,
            )

        result = dict(variants)

        if return_metadata:
            metadata_df = pd.DataFrame(records)
            if register:
                return result, metadata_df, dataset_path
            return result, metadata_df

        if register:
            return result, dataset_path

        return result

    def compute_conservation(
        self,
        sequences: Optional[Dict[str, str]] = None,
        *,
        dataset_name: Optional[str] = None,
        alphabet: Optional[Iterable[str]] = None,
        ignore_gaps: bool = True,
        gap_characters: Optional[Iterable[str]] = None,
        pseudocount: float = 0.0,
        normalize_entropy: bool = True,
        store_result: bool = False,
        result_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Compute per-position conservation metrics for a sequence collection."""

        seq_map = self._collect_sequences(sequences, dataset_name)
        if not seq_map:
            raise ValueError("No sequences available for conservation analysis")

        lengths = {len(seq) for seq in seq_map.values()}
        if len(lengths) != 1:
            raise ValueError("Sequences must be aligned (equal length) for conservation analysis")

        seq_length = lengths.pop()
        matrix = pd.DataFrame(
            [list(seq) for seq in seq_map.values()],
            index=list(seq_map.keys()),
        )

        gap_set = set(gap_characters or ['-', '.'])

        if alphabet is None:
            observed = set()
            for seq in seq_map.values():
                observed.update(ch for ch in seq if ch not in gap_set)
            alphabet_set = sorted(observed) if observed else STANDARD_AMINO_ACIDS
        else:
            alphabet_set = sorted({aa.upper() for aa in alphabet})

        records = []
        max_entropy = math.log2(len(alphabet_set)) if alphabet_set else 1.0

        for position in range(seq_length):
            column = matrix.iloc[:, position]
            if ignore_gaps:
                column = column[~column.isin(gap_set)]

            total = len(column)
            if total == 0:
                records.append({
                    'position': position + 1,
                    'consensus': '',
                    'consensus_frequency': 0.0,
                    'entropy': 0.0,
                    'normalized_entropy': 0.0,
                    'counts': {},
                })
                continue

            counts_series = column.value_counts().astype(float)
            all_counts = {aa: counts_series.get(aa, 0.0) for aa in alphabet_set}

            if pseudocount > 0:
                all_counts = {aa: count + pseudocount for aa, count in all_counts.items()}

            count_values = np.array(list(all_counts.values()), dtype=float)
            freq_values = count_values / count_values.sum()

            consensus_idx = int(np.argmax(freq_values))
            consensus_residue = alphabet_set[consensus_idx]
            consensus_frequency = float(freq_values[consensus_idx])

            nonzero = freq_values[freq_values > 0]
            entropy = float(-(nonzero * np.log2(nonzero)).sum()) if nonzero.size else 0.0
            normalized = entropy / max_entropy if (normalize_entropy and max_entropy > 0) else entropy

            records.append({
                'position': position + 1,
                'consensus': consensus_residue,
                'consensus_frequency': consensus_frequency,
                'entropy': entropy,
                'normalized_entropy': normalized,
                'counts': {aa: int(counts_series.get(aa, 0)) for aa in alphabet_set},
            })

        result_df = pd.DataFrame(records)

        if store_result:
            key = result_name or dataset_name or 'conservation'
            self.metadata[key] = result_df

        return result_df

    def compute_linkage(
        self,
        sequences: Optional[Dict[str, str]] = None,
        *,
        dataset_name: Optional[str] = None,
        ignore_gaps: bool = True,
        gap_characters: Optional[Iterable[str]] = None,
        min_observations: int = 5,
        normalize: bool = True,
        top_k: Optional[int] = None,
        store_result: bool = False,
        result_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Compute pairwise residue linkage using mutual information."""

        seq_map = self._collect_sequences(sequences, dataset_name)
        if not seq_map:
            raise ValueError("No sequences available for linkage analysis")

        lengths = {len(seq) for seq in seq_map.values()}
        if len(lengths) != 1:
            raise ValueError("Sequences must be aligned (equal length) for linkage analysis")

        seq_length = lengths.pop()
        if seq_length < 2:
            raise ValueError("Linkage analysis requires sequences of length >= 2")

        seq_array = np.array([list(seq) for seq in seq_map.values()])
        gap_set = set(gap_characters or ['-', '.'])

        rows: List[Dict[str, Any]] = []

        for i in range(seq_length):
            col_i = seq_array[:, i]
            for j in range(i + 1, seq_length):
                col_j = seq_array[:, j]

                pair_counts: Counter = Counter()
                for a, b in zip(col_i, col_j):
                    if ignore_gaps and (a in gap_set or b in gap_set):
                        continue
                    pair_counts[(a, b)] += 1

                total = sum(pair_counts.values())
                if total < min_observations:
                    continue

                px: Counter = Counter()
                py: Counter = Counter()
                for (a, b), count in pair_counts.items():
                    px[a] += count
                    py[b] += count

                mi = 0.0
                for (a, b), count in pair_counts.items():
                    pxy = count / total
                    px_prob = px[a] / total
                    py_prob = py[b] / total
                    if pxy <= 0 or px_prob <= 0 or py_prob <= 0:
                        continue
                    mi += pxy * math.log2(pxy / (px_prob * py_prob))

                hx = -sum(
                    (count / total) * math.log2(count / total)
                    for count in px.values()
                    if count > 0
                )
                hy = -sum(
                    (count / total) * math.log2(count / total)
                    for count in py.values()
                    if count > 0
                )

                if normalize and min(hx, hy) > 0:
                    normalized_mi = mi / min(hx, hy)
                else:
                    normalized_mi = mi

                top_pair, top_count = max(pair_counts.items(), key=lambda item: item[1])

                rows.append({
                    'pos_i': i + 1,
                    'pos_j': j + 1,
                    'mutual_information': mi,
                    'normalized_mi': normalized_mi,
                    'top_pair': f"{top_pair[0]}-{top_pair[1]}",
                    'top_pair_fraction': top_count / total,
                    'observations': total,
                })

        if not rows:
            return pd.DataFrame(
                columns=[
                    'pos_i',
                    'pos_j',
                    'mutual_information',
                    'normalized_mi',
                    'top_pair',
                    'top_pair_fraction',
                    'observations',
                ]
            )

        linkage_df = pd.DataFrame(rows)
        linkage_df = linkage_df.sort_values('mutual_information', ascending=False)

        if top_k is not None:
            linkage_df = linkage_df.head(top_k)

        if store_result:
            key = result_name or dataset_name or 'sequence_linkage'
            self.metadata[key] = linkage_df

        return linkage_df

    def annotate_with_grn_reference(
        self,
        dataset_name: str,
        reference_table: str,
        *,
        output_table_name: Optional[str] = None,
        allow_create: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> pd.DataFrame:
        """Align sequences to a GRN reference table and record the resulting GRN annotations."""

        sequences = self.load_dataset(dataset_name)
        if not sequences:
            raise ValueError(f"Dataset '{dataset_name}' contained no sequences")

        grn_proc = GRNProcessor()
        reference_df = grn_proc.load_reference_table(reference_table)
        if reference_df.empty:
            raise ValueError(f"Reference table '{reference_table}' is empty")

        ref_sequences: Dict[str, str] = {}
        ref_valid_columns: Dict[str, List[str]] = {}
        for ref_id in reference_df.index:
            seq = get_seq(ref_id, reference_df)
            if not seq:
                continue
            ref_sequences[ref_id] = seq
            valid_cols = [
                col
                for col, val in reference_df.loc[ref_id].items()
                if val not in {GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL, '-', 'X'}
            ]
            ref_valid_columns[ref_id] = valid_cols

        if not ref_sequences:
            raise ValueError("Reference table did not contain valid sequences")

        engine = SequenceAlignmentEngine()
        table_rows: Dict[str, pd.Series] = {}

        for seq_id, sequence in sequences.items():
            best_ref = None
            best_result = None
            best_score = float('-inf')

            for ref_id, ref_seq in ref_sequences.items():
                result = engine.align_pairwise(ref_id, ref_seq, seq_id, sequence)
                score = result.score / max(len(sequence), 1)
                if score > best_score:
                    best_score = score
                    best_ref = ref_id
                    best_result = result

            if best_ref is None or best_result is None:
                continue

            reference_row = reference_df.loc[best_ref]
            valid_columns = ref_valid_columns.get(best_ref, [])
            if not valid_columns:
                continue

            row = pd.Series('-', index=reference_row.index, dtype=object)
            ref_chars = best_result.alignment[0]
            query_chars = best_result.alignment[2]

            ref_positions = [
                (col, reference_row[col])
                for col in reference_row.index
                if reference_row[col] not in {GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL, '-', 'X'}
            ]

            ref_pointer = 0
            query_position = 0

            for ref_char, query_char in zip(ref_chars, query_chars):
                if ref_char != '-':
                    if ref_pointer >= len(ref_positions):
                        break
                    grn_col, _ = ref_positions[ref_pointer]
                    ref_pointer += 1
                    if query_char != '-':
                        query_position += 1
                        row[grn_col] = f"{query_char}{query_position}"
                    else:
                        row[grn_col] = '-'
                else:
                    if query_char != '-':
                        query_position += 1

            table_rows[seq_id] = row

        if not table_rows:
            raise ValueError("No sequences could be annotated against the reference table")

        table_df = pd.DataFrame(table_rows).T
        table_df.index.name = 'entity_name'

        output_name = output_table_name or f"{dataset_name}_{reference_table}_grn"
        table_metadata = {
            'dataset_name': dataset_name,
            'reference': reference_table,
            'best_match_count': len(table_df),
        }
        if metadata:
            table_metadata.update(metadata)

        grn_proc.record_table(
            output_name,
            table_df,
            metadata=table_metadata,
            allow_create=allow_create,
        )

        return table_df

    def save_alignment(self, alignment_data: Dict, output_file: str, 
                      alignment_type: str = "pairwise"):
        """
        Save alignment data to appropriate subdirectory.
        
        Args:
            alignment_data: Alignment data to save
            output_file: Output filename
            alignment_type: Type of alignment ("pairwise", "multiple", "mmseqs")
        """
        # Determine appropriate directory based on alignment type
        if alignment_type == "pairwise":
            output_path = self.path_pairwise_alignments_dir / output_file
        elif alignment_type == "multiple":
            output_path = self.path_multiple_alignments_dir / output_file
        elif alignment_type == "mmseqs":
            output_path = self.path_mmseqs_alignments_dir / output_file
        else:
            # Default to general alignments directory
            output_path = self.path_alignments_dir / output_file
        
        # ProtosPaths handles directory creation
        save_json(alignment_data, str(output_path))
        logger.info(f"Saved {alignment_type} alignment to {output_path}")
    
    def load_alignment(self, alignment_file: str, alignment_type: Optional[str] = None) -> Dict:
        """
        Load alignment data from file.
        
        Args:
            alignment_file: Alignment filename
            alignment_type: Type of alignment (if None, searches all subdirectories)
            
        Returns:
            Alignment data dictionary
        """
        # If alignment type specified, look in specific directory
        if alignment_type:
            if alignment_type == "pairwise":
                alignment_path = self.path_pairwise_alignments_dir / alignment_file
            elif alignment_type == "multiple":
                alignment_path = self.path_multiple_alignments_dir / alignment_file
            elif alignment_type == "mmseqs":
                alignment_path = self.path_mmseqs_alignments_dir / alignment_file
            else:
                alignment_path = self.path_alignments_dir / alignment_file
                
            if alignment_path.exists():
                return load_json(str(alignment_path))
            else:
                raise FileNotFoundError(f"Alignment file not found: {alignment_path}")
        
        # Otherwise, search all alignment directories
        search_dirs = [
            self.path_alignments_dir,
            self.path_pairwise_alignments_dir,
            self.path_multiple_alignments_dir,
            self.path_mmseqs_alignments_dir
        ]
        
        for dir_path in search_dirs:
            alignment_path = dir_path / alignment_file
            if alignment_path.exists():
                return load_json(str(alignment_path))
                
        raise FileNotFoundError(f"Alignment file not found in any alignment directory: {alignment_file}")
    
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
    
    def create_mmseqs_database(self, sequences: Dict[str, str], db_name: str) -> str:
        """
        Create MMseqs2 database from sequences.
        
        Args:
            sequences: Dictionary of sequence_id -> sequence
            db_name: Name for the database
            
        Returns:
            Path to the created database
        """
        # Save sequences to temporary FASTA
        temp_fasta = self.path_databases_dir / f"{db_name}_temp.fasta"
        write_fasta(sequences, str(temp_fasta))
        
        # Create MMseqs2 database
        db_path = self.path_databases_dir / db_name
        
        try:
            # Import mmseqs_helper if available
            from .mmseqs_helper import create_database
            create_database(str(temp_fasta), str(db_path))
            logger.info(f"Created MMseqs2 database: {db_path}")
        except ImportError:
            logger.warning("MMseqs2 helper not available. Database creation skipped.")
        finally:
            # Clean up temporary file
            if temp_fasta.exists():
                temp_fasta.unlink()
                
        return str(db_path)
    
    def list_mmseqs_databases(self) -> List[str]:
        """List available MMseqs2 databases."""
        db_dir = self.path_databases_dir
        if not db_dir.exists():
            return []
            
        # MMseqs2 databases have multiple files, look for .dbtype files
        db_files = list(db_dir.glob("*.dbtype"))
        return [f.stem for f in db_files]
    
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

    def export_dataset(
        self,
        dataset_name: str,
        output_dir: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
        name_pattern: Optional[str] = None,
        sequence_ids: Optional[List[str]] = None,
    ) -> Dict[str, Path]:
        exporter = self._get_exporter()
        return exporter.export_dataset(
            dataset_name,
            output_dir,
            format=format,
            overwrite=overwrite,
            name_pattern=name_pattern,
            sequence_ids=sequence_ids,
        )

    def _get_exporter(self):
        from protos.io.export.sequence_exporter import SequenceExporter

        return SequenceExporter(self)
    
    def list_entities(self, entity_type: Optional[str] = None) -> List[str]:
        """
        List all available sequence entities.
        
        Args:
            entity_type: Optional filter by type ('sequence' or 'sequence_dataset')
            
        Returns:
            List of entity names
        """
        return BaseProcessor.list_entities(self, entity_type)
    
    def list_datasets(self) -> List[str]:
        """
        List available sequence datasets (multi-sequence FASTA files).
        
        Returns:
            List of dataset names
        """
        return BaseProcessor.list_datasets(self)
from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

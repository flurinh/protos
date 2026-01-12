"""Loader for sequence data from NCBI (Entrez and BLAST).

This module provides:
- NCBILoader: BaseLoader implementation for NCBI sequences
- Fetching sequences by accession via Entrez
- BLAST homology search
- Batch sequence downloads
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from protos.io.core.base_loader import BaseLoader
from protos.io.formats.fasta_utils import write_fasta
from protos.io.ingest.utils.ncbi_utils import (
    BlastResult,
    fetch_sequence_entrez,
    fetch_sequences_batch,
    parse_fasta_text,
    run_blast_search,
    search_entrez,
)


class NCBILoader(BaseLoader):
    """Loader for fetching sequences from NCBI.

    Supports:
    - Fetching sequences by NCBI accession (via Entrez)
    - Running BLAST searches
    - Batch sequence downloads
    - Integration with SequenceProcessor for registration

    Example usage:
        loader = NCBILoader()

        # Fetch single sequence
        loader.download_and_register("NP_000530.1", name="RHO_HUMAN")

        # Batch fetch
        loader.download_batch(
            ["P02699", "P08100", "P15409"],
            dataset_name="rhodopsins"
        )

        # BLAST search
        results = loader.blast_search(
            sequence="MNGTEGPNFYVPFSN...",
            database="swissprot"
        )
    """

    loader_type = "ncbi"

    def __init__(
        self,
        name: str = "ncbi_loader",
        *,
        processor: Optional[Any] = None,
        api_key: Optional[str] = None,
    ) -> None:
        """Initialize NCBILoader.

        Args:
            name: Loader instance name
            processor: Optional SequenceProcessor for registration
            api_key: Optional NCBI API key for higher rate limits
        """
        super().__init__(name=name)
        self._processor = processor
        self.api_key = api_key

    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        """Parse an NCBI identifier.

        Supported formats:
        - Direct accession: "NP_000530.1", "P02699"
        - Prefixed: "ncbi:NP_000530.1", "refseq:NP_000530"
        - Multiple: "ncbi:P02699,P08100,P15409"
        - Local path: "/path/to/file.fasta"

        Args:
            identifier: Identifier string to parse

        Returns:
            Dictionary with source, ids, and name
        """
        path = Path(identifier)
        if path.exists():
            return {
                "source": "local",
                "path": path,
                "name": path.stem,
            }

        token = identifier.strip()

        # Handle prefixed identifiers
        if ":" in token:
            prefix, ids_str = token.split(":", 1)
            prefix = prefix.lower()

            if prefix in ("ncbi", "refseq", "genbank", "protein"):
                ids = [part.strip() for part in ids_str.split(",") if part.strip()]
                if not ids:
                    raise ValueError("No NCBI identifiers provided")
                return {"source": "ncbi", "ids": ids, "name": ids[0]}

        # Try to detect if it looks like an NCBI accession
        if self._looks_like_ncbi_accession(token):
            return {"source": "ncbi", "ids": [token], "name": token}

        raise ValueError(f"Unsupported NCBI identifier: {identifier}")

    def fetch_entity(
        self,
        identifier: str,
        source: Optional[str] = None,
        db: str = "protein",
        **kwargs,
    ) -> Optional[Path]:
        """Fetch a sequence from NCBI.

        Args:
            identifier: NCBI accession or identifier string
            source: Source type (auto-detected if None)
            db: NCBI database (protein, nucleotide)

        Returns:
            Path to downloaded FASTA file, or None if failed
        """
        info = self.parse_identifier(identifier)
        source = info["source"]

        if source == "local":
            return info["path"]

        if source == "ncbi":
            return self._fetch_from_ncbi(info["ids"], info.get("name"), db=db)

        raise ValueError(f"Unsupported source for NCBI loader: {source}")

    def download_and_register(
        self,
        identifier: str,
        name: Optional[str] = None,
        source: str = "ncbi",
        materialize_entities: bool = True,
        metadata: Optional[Dict[str, Any]] = None,
        db: str = "protein",
        **kwargs,
    ) -> Optional[str]:
        """Download sequence from NCBI and register with SequenceProcessor.

        Args:
            identifier: NCBI accession or identifier
            name: Optional name for the sequence/dataset
            source: Source type
            materialize_entities: Whether to save individual entity files
            metadata: Optional metadata dictionary
            db: NCBI database

        Returns:
            Registered entity/dataset name, or None if failed
        """
        from protos.io.formats.fasta_utils import read_fasta

        info = self.parse_identifier(identifier)
        source_path = self.fetch_entity(identifier, source=info["source"], db=db)

        if source_path is None:
            return None

        sequences = read_fasta(str(source_path))
        if not sequences:
            return None

        processor = self._get_processor()

        if len(sequences) == 1 and materialize_entities:
            seq_id, sequence = next(iter(sequences.items()))
            entity_name = name or seq_id
            processor.save_entity(entity_name, sequence, metadata=metadata)
            saved_name = entity_name
        else:
            dataset_name = name or source_path.stem
            processor.save_sequences(
                sequences,
                output_file=dataset_name,
                dataset_name=dataset_name,
                metadata=metadata,
                materialize_entities=materialize_entities,
            )
            saved_name = dataset_name

        # Cleanup temp file
        if source_path.exists():
            try:
                source_path.unlink()
            except OSError:
                pass

        return saved_name

    def download_batch(
        self,
        accessions: List[str],
        dataset_name: str,
        db: str = "protein",
        metadata: Optional[Dict[str, Any]] = None,
        materialize_entities: bool = True,
    ) -> Dict[str, Any]:
        """Download multiple sequences from NCBI and register as dataset.

        Args:
            accessions: List of NCBI accessions
            dataset_name: Name for the dataset
            db: NCBI database
            metadata: Optional metadata
            materialize_entities: Whether to save individual entities

        Returns:
            Dictionary with download results
        """
        results = {
            "success": [],
            "failed": [],
            "dataset": dataset_name,
        }

        # Fetch sequences
        sequences = {}
        for acc in accessions:
            content = fetch_sequence_entrez(acc, db=db, api_key=self.api_key)
            if content:
                parsed = parse_fasta_text(content)
                sequences.update(parsed)
                results["success"].append(acc)
            else:
                results["failed"].append(acc)

        if not sequences:
            return results

        # Register with processor
        processor = self._get_processor()
        ds_metadata = metadata.copy() if metadata else {}
        ds_metadata["source"] = "ncbi"
        ds_metadata["database"] = db
        ds_metadata["accessions"] = results["success"]

        processor.save_sequences(
            sequences,
            output_file=dataset_name,
            dataset_name=dataset_name,
            metadata=ds_metadata,
            materialize_entities=materialize_entities,
        )

        return results

    def blast_search(
        self,
        sequence: str,
        query_id: str = "query",
        program: str = "blastp",
        database: str = "swissprot",
        hitlist_size: int = 50,
        expect: float = 10.0,
        max_wait: int = 300,
    ) -> Optional[BlastResult]:
        """Run a BLAST search against NCBI databases.

        Args:
            sequence: Query sequence (amino acid or nucleotide)
            query_id: Identifier for the query
            program: BLAST program (blastp, blastn, etc.)
            database: NCBI database (nr, refseq_protein, swissprot, pdb)
            hitlist_size: Maximum number of hits
            expect: E-value threshold
            max_wait: Maximum seconds to wait for results

        Returns:
            BlastResult object with parsed hits, or None if failed
        """
        return run_blast_search(
            sequence=sequence,
            query_id=query_id,
            program=program,
            database=database,
            hitlist_size=hitlist_size,
            expect=expect,
            api_key=self.api_key,
            max_wait=max_wait,
        )

    def blast_search_and_fetch(
        self,
        sequence: str,
        query_id: str = "query",
        database: str = "swissprot",
        hitlist_size: int = 20,
        dataset_name: Optional[str] = None,
        min_identity: float = 30.0,
        min_coverage: float = 50.0,
    ) -> Optional[Dict[str, Any]]:
        """Run BLAST search and fetch hit sequences.

        Args:
            sequence: Query sequence
            query_id: Query identifier
            database: NCBI database
            hitlist_size: Maximum hits
            dataset_name: Optional dataset name for registration
            min_identity: Minimum identity percentage to fetch
            min_coverage: Minimum coverage percentage to fetch

        Returns:
            Dictionary with blast results and fetched sequences
        """
        # Run BLAST
        result = self.blast_search(
            sequence=sequence,
            query_id=query_id,
            database=database,
            hitlist_size=hitlist_size,
        )

        if not result:
            return None

        # Filter hits
        filtered_hits = [
            hit for hit in result.hits
            if hit.identity_percent >= min_identity
            and hit.coverage_percent >= min_coverage
        ]

        if not filtered_hits:
            return {
                "blast_result": result,
                "sequences": {},
                "filtered_hits": 0,
            }

        # Fetch hit sequences
        accessions = [hit.hit_accession for hit in filtered_hits]
        sequences = {}

        for acc in accessions[:hitlist_size]:  # Limit fetches
            content = fetch_sequence_entrez(acc, api_key=self.api_key)
            if content:
                parsed = parse_fasta_text(content)
                sequences.update(parsed)

        # Register if dataset name provided
        if dataset_name and sequences:
            processor = self._get_processor()
            processor.save_sequences(
                sequences,
                output_file=dataset_name,
                dataset_name=dataset_name,
                metadata={
                    "source": "ncbi_blast",
                    "query_id": query_id,
                    "database": database,
                    "min_identity": min_identity,
                    "min_coverage": min_coverage,
                },
                materialize_entities=True,
            )

        return {
            "blast_result": result,
            "sequences": sequences,
            "filtered_hits": len(filtered_hits),
            "dataset": dataset_name,
        }

    def search_and_download(
        self,
        query: str,
        dataset_name: str,
        db: str = "protein",
        retmax: int = 100,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Search NCBI Entrez and download matching sequences.

        Args:
            query: Entrez search query (e.g., "rhodopsin[gene] AND human[organism]")
            dataset_name: Name for the dataset
            db: NCBI database
            retmax: Maximum number of results
            metadata: Optional metadata

        Returns:
            Dictionary with search and download results
        """
        # Search
        ids = search_entrez(query, db=db, retmax=retmax, api_key=self.api_key)

        if not ids:
            return {"success": [], "failed": [], "dataset": None}

        # Download
        return self.download_batch(
            accessions=ids,
            dataset_name=dataset_name,
            db=db,
            metadata=metadata,
        )

    def to_dataframe(self, result: BlastResult) -> pd.DataFrame:
        """Convert BLAST results to pandas DataFrame.

        Args:
            result: BlastResult from blast_search

        Returns:
            DataFrame with hit information
        """
        records = []
        for hit in result.hits:
            records.append({
                "query_id": result.query_id,
                "hit_id": hit.hit_id,
                "hit_accession": hit.hit_accession,
                "hit_def": hit.hit_def[:100],  # Truncate long descriptions
                "identity_percent": hit.identity_percent,
                "coverage_percent": hit.coverage_percent,
                "evalue": hit.hsp_evalue,
                "bit_score": hit.hsp_bit_score,
                "align_len": hit.hsp_align_len,
                "gaps": hit.hsp_gaps,
            })

        return pd.DataFrame(records)

    def _get_processor(self):
        """Get or create SequenceProcessor."""
        if self._processor is None:
            from protos.processing.sequence import SequenceProcessor
            self._processor = SequenceProcessor(name=f"{self.name}_processor")
        return self._processor

    def _fetch_from_ncbi(
        self,
        ids: List[str],
        target_name: Optional[str],
        db: str = "protein",
    ) -> Optional[Path]:
        """Fetch sequences from NCBI and save to temp file.

        Args:
            ids: List of NCBI accessions
            target_name: Target filename
            db: NCBI database

        Returns:
            Path to FASTA file, or None if all failed
        """
        processor = self._get_processor()
        safe_name = processor._sanitize_filename(target_name or ids[0])
        output_path = Path(processor.path_fasta_dir) / f"{safe_name}.fasta"

        sequences = {}
        for ncbi_id in ids:
            content = fetch_sequence_entrez(ncbi_id, db=db, api_key=self.api_key)
            if content:
                parsed = parse_fasta_text(content)
                sequences.update(parsed)

        if not sequences:
            return None

        write_fasta(sequences, str(output_path))
        return output_path

    @staticmethod
    def _looks_like_ncbi_accession(identifier: str) -> bool:
        """Check if identifier looks like an NCBI accession.

        Patterns:
        - RefSeq: NP_000530.1, NM_000001.2, XP_123456.1
        - GenBank: AAA12345.1
        - UniProt-like: P02699 (also valid for NCBI)
        """
        if not identifier or len(identifier) < 5:
            return False

        # RefSeq patterns
        if identifier[:2] in ("NP", "NM", "XP", "XM", "YP", "WP", "NC", "NG"):
            return True

        # GenBank protein
        if identifier[:3].isalpha() and identifier[3:8].isdigit():
            return True

        # UniProt-style (also works with NCBI)
        if identifier[0].isalpha() and len(identifier) >= 6:
            return True

        return False

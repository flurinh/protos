"""NCBI utilities for sequence fetching and BLAST searches.

This module provides functions for:
- Fetching sequences from NCBI Entrez (by accession)
- Running NCBI BLAST searches
- Parsing BLAST results
"""

from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# NCBI API endpoints
ENTREZ_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BLAST_BASE_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# Rate limiting: NCBI allows 3 requests/second without API key, 10 with key
DEFAULT_DELAY = 0.34  # ~3 requests per second

# Setup session with retries
retries = Retry(total=5, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


@dataclass
class BlastHit:
    """Represents a single BLAST hit."""
    hit_id: str
    hit_def: str
    hit_accession: str
    hit_len: int
    hsp_bit_score: float
    hsp_score: int
    hsp_evalue: float
    hsp_identity: int
    hsp_positive: int
    hsp_gaps: int
    hsp_align_len: int
    hsp_query_from: int
    hsp_query_to: int
    hsp_hit_from: int
    hsp_hit_to: int
    hsp_qseq: str
    hsp_hseq: str
    hsp_midline: str
    identity_percent: float
    coverage_percent: float


@dataclass
class BlastResult:
    """Container for BLAST search results."""
    query_id: str
    query_len: int
    database: str
    hits: List[BlastHit]
    request_id: str
    execution_time: float


def fetch_sequence_entrez(
    accession: str,
    db: str = "protein",
    rettype: str = "fasta",
    api_key: Optional[str] = None,
) -> Optional[str]:
    """Fetch a sequence from NCBI Entrez by accession.

    Args:
        accession: NCBI accession number (e.g., NP_000530.1, P02699)
        db: Database to search (protein, nucleotide, etc.)
        rettype: Return type (fasta, gb, etc.)
        api_key: Optional NCBI API key for higher rate limits

    Returns:
        Sequence in requested format, or None if not found
    """
    url = f"{ENTREZ_BASE_URL}/efetch.fcgi"
    params = {
        "db": db,
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()

        content = response.text.strip()
        if not content or "Error" in content[:100]:
            return None
        return content

    except requests.RequestException as e:
        print(f"Error fetching {accession}: {e}")
        return None


def fetch_sequences_batch(
    accessions: List[str],
    db: str = "protein",
    rettype: str = "fasta",
    api_key: Optional[str] = None,
    delay: float = DEFAULT_DELAY,
) -> Dict[str, str]:
    """Fetch multiple sequences from NCBI Entrez.

    Args:
        accessions: List of NCBI accession numbers
        db: Database to search
        rettype: Return type
        api_key: Optional NCBI API key
        delay: Delay between requests (seconds)

    Returns:
        Dictionary mapping accession to sequence (FASTA format)
    """
    results = {}

    for i, acc in enumerate(accessions):
        if i > 0:
            time.sleep(delay)

        content = fetch_sequence_entrez(acc, db=db, rettype=rettype, api_key=api_key)
        if content:
            results[acc] = content
        else:
            print(f"Warning: Could not fetch {acc}")

    return results


def search_entrez(
    query: str,
    db: str = "protein",
    retmax: int = 100,
    api_key: Optional[str] = None,
) -> List[str]:
    """Search NCBI Entrez and return accession IDs.

    Args:
        query: Search query (e.g., "rhodopsin[gene] AND human[organism]")
        db: Database to search
        retmax: Maximum number of results
        api_key: Optional NCBI API key

    Returns:
        List of accession IDs
    """
    url = f"{ENTREZ_BASE_URL}/esearch.fcgi"
    params = {
        "db": db,
        "term": query,
        "retmax": retmax,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        return data.get("esearchresult", {}).get("idlist", [])

    except requests.RequestException as e:
        print(f"Error searching Entrez: {e}")
        return []


def submit_blast(
    sequence: str,
    program: str = "blastp",
    database: str = "nr",
    hitlist_size: int = 50,
    expect: float = 10.0,
    api_key: Optional[str] = None,
) -> Optional[str]:
    """Submit a BLAST search to NCBI.

    Args:
        sequence: Query sequence (amino acid or nucleotide)
        program: BLAST program (blastp, blastn, blastx, tblastn, tblastx)
        database: Database to search (nr, refseq_protein, swissprot, pdb, etc.)
        hitlist_size: Maximum number of hits to return
        expect: E-value threshold
        api_key: Optional NCBI API key

    Returns:
        Request ID (RID) for retrieving results, or None if submission failed
    """
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "HITLIST_SIZE": hitlist_size,
        "EXPECT": expect,
        "FORMAT_TYPE": "XML",
    }
    if api_key:
        params["API_KEY"] = api_key

    try:
        response = session.post(BLAST_BASE_URL, data=params, timeout=60)
        response.raise_for_status()

        # Parse RID from response
        content = response.text
        for line in content.split("\n"):
            if line.strip().startswith("RID = "):
                return line.split("=")[1].strip()

        print("Warning: Could not parse RID from BLAST response")
        return None

    except requests.RequestException as e:
        print(f"Error submitting BLAST: {e}")
        return None


def check_blast_status(rid: str, api_key: Optional[str] = None) -> Tuple[str, Optional[str]]:
    """Check the status of a BLAST search.

    Args:
        rid: Request ID from submit_blast
        api_key: Optional NCBI API key

    Returns:
        Tuple of (status, message) where status is one of:
        - "WAITING": Still running
        - "READY": Results available
        - "UNKNOWN": Unknown RID
        - "ERROR": Error occurred
    """
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo",
    }
    if api_key:
        params["API_KEY"] = api_key

    try:
        response = session.get(BLAST_BASE_URL, params=params, timeout=30)
        response.raise_for_status()

        content = response.text

        if "Status=WAITING" in content:
            return "WAITING", None
        elif "Status=READY" in content:
            if "ThereAreHits=yes" in content:
                return "READY", "hits_found"
            else:
                return "READY", "no_hits"
        elif "Status=UNKNOWN" in content:
            return "UNKNOWN", "Invalid RID"
        else:
            return "ERROR", content[:200]

    except requests.RequestException as e:
        return "ERROR", str(e)


def get_blast_results(
    rid: str,
    api_key: Optional[str] = None,
) -> Optional[str]:
    """Retrieve BLAST results in XML format.

    Args:
        rid: Request ID from submit_blast
        api_key: Optional NCBI API key

    Returns:
        XML string with BLAST results, or None if error
    """
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
    }
    if api_key:
        params["API_KEY"] = api_key

    try:
        response = session.get(BLAST_BASE_URL, params=params, timeout=120)
        response.raise_for_status()
        return response.text

    except requests.RequestException as e:
        print(f"Error retrieving BLAST results: {e}")
        return None


def parse_blast_xml(xml_content: str, query_len: int = 0) -> List[BlastHit]:
    """Parse BLAST XML output into BlastHit objects.

    Args:
        xml_content: XML string from BLAST results
        query_len: Length of query sequence (for coverage calculation)

    Returns:
        List of BlastHit objects
    """
    hits = []

    try:
        root = ET.fromstring(xml_content)

        # Find all Hit elements
        for hit in root.iter("Hit"):
            hit_id = hit.findtext("Hit_id", "")
            hit_def = hit.findtext("Hit_def", "")
            hit_accession = hit.findtext("Hit_accession", "")
            hit_len = int(hit.findtext("Hit_len", "0"))

            # Get best HSP (first one)
            hsp = hit.find(".//Hsp")
            if hsp is None:
                continue

            hsp_bit_score = float(hsp.findtext("Hsp_bit-score", "0"))
            hsp_score = int(hsp.findtext("Hsp_score", "0"))
            hsp_evalue = float(hsp.findtext("Hsp_evalue", "1e10"))
            hsp_identity = int(hsp.findtext("Hsp_identity", "0"))
            hsp_positive = int(hsp.findtext("Hsp_positive", "0"))
            hsp_gaps = int(hsp.findtext("Hsp_gaps", "0"))
            hsp_align_len = int(hsp.findtext("Hsp_align-len", "1"))
            hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
            hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
            hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
            hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
            hsp_qseq = hsp.findtext("Hsp_qseq", "")
            hsp_hseq = hsp.findtext("Hsp_hseq", "")
            hsp_midline = hsp.findtext("Hsp_midline", "")

            # Calculate percentages
            identity_percent = (hsp_identity / hsp_align_len * 100) if hsp_align_len > 0 else 0
            coverage = hsp_query_to - hsp_query_from + 1
            coverage_percent = (coverage / query_len * 100) if query_len > 0 else 0

            hits.append(BlastHit(
                hit_id=hit_id,
                hit_def=hit_def,
                hit_accession=hit_accession,
                hit_len=hit_len,
                hsp_bit_score=hsp_bit_score,
                hsp_score=hsp_score,
                hsp_evalue=hsp_evalue,
                hsp_identity=hsp_identity,
                hsp_positive=hsp_positive,
                hsp_gaps=hsp_gaps,
                hsp_align_len=hsp_align_len,
                hsp_query_from=hsp_query_from,
                hsp_query_to=hsp_query_to,
                hsp_hit_from=hsp_hit_from,
                hsp_hit_to=hsp_hit_to,
                hsp_qseq=hsp_qseq,
                hsp_hseq=hsp_hseq,
                hsp_midline=hsp_midline,
                identity_percent=identity_percent,
                coverage_percent=coverage_percent,
            ))

    except ET.ParseError as e:
        print(f"Error parsing BLAST XML: {e}")

    return hits


def run_blast_search(
    sequence: str,
    query_id: str = "query",
    program: str = "blastp",
    database: str = "swissprot",
    hitlist_size: int = 50,
    expect: float = 10.0,
    api_key: Optional[str] = None,
    poll_interval: int = 10,
    max_wait: int = 300,
) -> Optional[BlastResult]:
    """Run a complete BLAST search and return parsed results.

    This function submits the search, polls for completion, and parses results.

    Args:
        sequence: Query sequence
        query_id: Identifier for the query
        program: BLAST program
        database: Database to search
        hitlist_size: Maximum hits
        expect: E-value threshold
        api_key: Optional NCBI API key
        poll_interval: Seconds between status checks
        max_wait: Maximum seconds to wait for results

    Returns:
        BlastResult object, or None if search failed
    """
    import time as time_module
    start_time = time_module.time()

    # Submit search
    print(f"Submitting BLAST search ({program} vs {database})...")
    rid = submit_blast(
        sequence=sequence,
        program=program,
        database=database,
        hitlist_size=hitlist_size,
        expect=expect,
        api_key=api_key,
    )

    if not rid:
        print("Failed to submit BLAST search")
        return None

    print(f"Request ID: {rid}")

    # Poll for completion
    elapsed = 0
    while elapsed < max_wait:
        time_module.sleep(poll_interval)
        elapsed = time_module.time() - start_time

        status, message = check_blast_status(rid, api_key=api_key)
        print(f"  Status: {status} ({elapsed:.0f}s elapsed)")

        if status == "READY":
            break
        elif status in ("UNKNOWN", "ERROR"):
            print(f"BLAST error: {message}")
            return None

    if elapsed >= max_wait:
        print(f"BLAST search timed out after {max_wait}s")
        return None

    # Get results
    print("Retrieving results...")
    xml_content = get_blast_results(rid, api_key=api_key)

    if not xml_content:
        print("Failed to retrieve BLAST results")
        return None

    # Parse results
    query_len = len(sequence)
    hits = parse_blast_xml(xml_content, query_len=query_len)

    execution_time = time_module.time() - start_time

    return BlastResult(
        query_id=query_id,
        query_len=query_len,
        database=database,
        hits=hits,
        request_id=rid,
        execution_time=execution_time,
    )


def parse_fasta_text(text: str) -> Dict[str, str]:
    """Parse FASTA format text into dictionary.

    Args:
        text: FASTA format text

    Returns:
        Dictionary mapping sequence ID to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []

    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(current_seq)
            # Parse header - take first word after >
            header = line[1:].strip()
            current_id = header.split()[0] if header else "unknown"
            current_seq = []
        else:
            current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences

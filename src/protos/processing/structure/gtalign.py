import numpy as np
import pandas as pd
import subprocess
import os
import re
import tempfile
import time
from typing import List, Dict, Union, Optional, Tuple, Any
from pathlib import Path


def format_alignment_result(alignment_data: Dict[str, Any], line_width: int = 80) -> str:
    """
    Format alignment results in a visually appealing way.

    Args:
        alignment_data: Dictionary containing alignment results from GTalign
        line_width: Number of characters per line for sequence display

    Returns:
        Formatted string representation of the alignment
    """
    # Check if we have valid alignment data
    if alignment_data is None:
        return "No alignment data provided"
    
    output = []

    # Add header with basic statistics
    output.append("=" * line_width)
    output.append(f"Alignment Results: {alignment_data.get('file', '')}")
    output.append("-" * line_width)
    output.append(f"TM-score: {alignment_data.get('tm_score', 0.0):.4f}")
    output.append(f"RMSD: {alignment_data.get('rmsd', 0.0):.4f} Å")
    output.append(f"Alignment Length: {alignment_data.get('alignment_length', 0)}")
    output.append(f"Query Length: {alignment_data.get('query_length', 0)}")
    output.append(f"Reference Length: {alignment_data.get('reference_length', 0)}")
    output.append("=" * line_width)

    # Format rotation matrix and translation vector
    if alignment_data.get('rotation_matrix') and alignment_data.get('translation_vector'):
        output.append("\nTransformation Matrix (Rotation [3,3] and Translation [3,1]):")
        rot_matrix = alignment_data['rotation_matrix']
        trans_vector = alignment_data['translation_vector']

        # Format each row of the transformation matrix
        for i in range(3):
            rot_str = " ".join([f"{rot_matrix[i][j]:10.6f}" for j in range(3)])
            output.append(f"  {rot_str}  |  {trans_vector[i]:10.6f}")

        output.append("")

    # Format sequence alignment in chunks
    query_seq = alignment_data.get('query_sequence_clean', '')
    ref_seq = alignment_data.get('reference_sequence_clean', '')
    match_line = alignment_data.get('match_line', '')

    # Ensure all three lines have the same length
    max_len = max(len(query_seq), len(ref_seq), len(match_line))
    query_seq = query_seq.ljust(max_len)
    ref_seq = ref_seq.ljust(max_len)
    match_line = match_line.ljust(max_len)

    output.append("Sequence Alignment:")
    output.append("-" * line_width)

    # Determine position numbering for query and reference
    query_positions = []
    ref_positions = []

    # If we have alignment blocks with positions, use them
    if alignment_data.get('alignment_blocks'):
        query_positions = [block['query_start'] for block in alignment_data['alignment_blocks']]
        ref_positions = [block['ref_start'] for block in alignment_data['alignment_blocks']]

    # Default to starting at position 1 if no blocks available
    query_pos = query_positions[0] if query_positions else 1
    ref_pos = ref_positions[0] if ref_positions else 1

    # Display alignment in chunks
    for i in range(0, max_len, line_width):
        chunk_end = min(i + line_width, max_len)
        query_chunk = query_seq[i:chunk_end]
        match_chunk = match_line[i:chunk_end]
        ref_chunk = ref_seq[i:chunk_end]

        # Calculate ending positions for this chunk
        query_ungapped = query_chunk.replace('-', '')
        ref_ungapped = ref_chunk.replace('-', '')
        query_chunk_end = query_pos + len(query_ungapped) - 1
        ref_chunk_end = ref_pos + len(ref_ungapped) - 1

        # Add alignment chunk to output
        output.append(f"Query: {query_pos:4d} {query_chunk} {query_chunk_end:4d}")
        output.append(f"Match:      {match_chunk}")
        output.append(f"Refn.: {ref_pos:4d} {ref_chunk} {ref_chunk_end:4d}")
        output.append("")

        # Update positions for next chunk
        query_pos = query_chunk_end + 1
        ref_pos = ref_chunk_end + 1

    # Add alignment path summary if available
    if alignment_data.get('alignment_path'):
        output.append("Alignment Path Summary:")
        output.append("-" * line_width)

        # Count aligned positions, gaps in query, and gaps in reference
        aligned_count = 0
        query_gaps = 0
        ref_gaps = 0

        for q_pos, r_pos in alignment_data['alignment_path']:
            if q_pos is None:
                query_gaps += 1
            elif r_pos is None:
                ref_gaps += 1
            else:
                aligned_count += 1

        total_positions = aligned_count + query_gaps + ref_gaps

        output.append(f"Aligned positions: {aligned_count} ({aligned_count / total_positions:.1%})")
        output.append(f"Gaps in query:     {query_gaps} ({query_gaps / total_positions:.1%})")
        output.append(f"Gaps in reference: {ref_gaps} ({ref_gaps / total_positions:.1%})")
        output.append(f"Total positions:   {total_positions}")

    return "\n".join(output)


def print_alignment_result(alignment_data: Dict[str, Any], line_width: int = 80) -> None:
    """
    Print formatted alignment results to the console.

    Args:
        alignment_data: Dictionary containing alignment results from GTalign
        line_width: Number of characters per line for sequence display
    """
    # Handle case where we're passed the full result dict instead of an alignment
    if 'alignments' in alignment_data:
        # We were passed the full result dict - find the first alignment with data
        for alignment in alignment_data.get('alignments', []):
            if alignment.get('tm_score') is not None:
                formatted_result = format_alignment_result(alignment, line_width)
                print(formatted_result)
                return
        print("No valid alignments found in result")
        return
    
    # Normal case - we were passed a single alignment
    formatted_result = format_alignment_result(alignment_data, line_width)
    print(formatted_result)


def create_alignment_path(query_seq, ref_seq, query_start, ref_start):
    """
    Create a position mapping between query and reference sequences.

    Args:
        query_seq: Query sequence with gaps (-)
        ref_seq: Reference sequence with gaps (-)
        query_start: Starting position of the query sequence (1-indexed)
        ref_start: Starting position of the reference sequence (1-indexed)

    Returns:
        List of tuples (query_pos, ref_pos) where:
        - query_pos is the position in the query (or None if gap)
        - ref_pos is the position in the reference (or None if gap)
    """
    if len(query_seq) != len(ref_seq):
        raise ValueError(f"Sequence lengths don't match: query={len(query_seq)}, ref={len(ref_seq)}")

    alignment_path = []
    query_pos = query_start
    ref_pos = ref_start

    for q_char, r_char in zip(query_seq, ref_seq):
        if q_char == '-' and r_char == '-':
            # Both are gaps (shouldn't happen in valid alignments)
            continue
        elif q_char == '-':
            # Gap in query
            alignment_path.append((None, ref_pos))
            ref_pos += 1
        elif r_char == '-':
            # Gap in reference
            alignment_path.append((query_pos, None))
            query_pos += 1
        else:
            # Match or mismatch
            alignment_path.append((query_pos, ref_pos))
            query_pos += 1
            ref_pos += 1

    return alignment_path


def create_full_alignment_path(alignment_data):
    """
    Create a complete alignment path from all alignment blocks.

    Args:
        alignment_data: Dictionary containing parsed alignment data

    Returns:
        List of tuples (query_pos, ref_pos) for the entire alignment
    """
    full_path = []

    # If we have alignment blocks, use them
    if alignment_data["alignment_blocks"]:
        for block in alignment_data["alignment_blocks"]:
            query_start = block["query_start"]
            ref_start = block["ref_start"]
            query_seq = block["query_seq"]
            ref_seq = block["ref_seq"]

            block_path = create_alignment_path(query_seq, ref_seq, query_start, ref_start)
            full_path.extend(block_path)
    # Otherwise try to use the full sequences
    elif alignment_data["query_sequence"] and alignment_data["reference_sequence"]:
        query_seq = alignment_data["query_sequence"]
        ref_seq = alignment_data["reference_sequence"]

        # Get the starting positions from positions list or default to 1
        query_start = alignment_data["query_positions"][0][0] if alignment_data["query_positions"] else 1
        ref_start = alignment_data["ref_positions"][0][0] if alignment_data["ref_positions"] else 1

        full_path = create_alignment_path(query_seq, ref_seq, query_start, ref_start)

    return full_path


def get_query_to_ref_mapping(alignment_path: List[Tuple[Optional[int], Optional[int]]]) -> Dict[int, Optional[int]]:
    """
    Create a mapping from query positions to reference positions.

    Args:
        alignment_path: List of (query_pos, ref_pos) tuples

    Returns:
        Dictionary mapping query positions to reference positions
        If a query position aligns to a gap, the value will be None
    """
    mapping = {}
    for q_pos, r_pos in alignment_path:
        if q_pos is not None:  # Skip gaps in query
            mapping[q_pos] = r_pos
    return mapping


def get_ref_to_query_mapping(alignment_path: List[Tuple[Optional[int], Optional[int]]]) -> Dict[int, Optional[int]]:
    """
    Create a mapping from reference positions to query positions.

    Args:
        alignment_path: List of (query_pos, ref_pos) tuples

    Returns:
        Dictionary mapping reference positions to query positions
        If a reference position aligns to a gap, the value will be None
    """
    mapping = {}
    for q_pos, r_pos in alignment_path:
        if r_pos is not None:  # Skip gaps in reference
            mapping[r_pos] = q_pos
    return mapping


def transform_coordinates(coordinates: np.ndarray, rotation_matrix: np.ndarray, 
                         translation_vector: np.ndarray) -> np.ndarray:
    """
    Apply rotation and translation to coordinates.
    
    Args:
        coordinates: Array of shape (n, 3) containing x, y, z coordinates
        rotation_matrix: 3x3 rotation matrix
        translation_vector: Translation vector of length 3
        
    Returns:
        Array of shape (n, 3) containing transformed coordinates
    """
    # Ensure arrays are numpy arrays with correct type
    coords = np.array(coordinates, dtype=float)
    rot_matrix = np.array(rotation_matrix, dtype=float)
    trans_vector = np.array(translation_vector, dtype=float)
    
    # Apply rotation (dot product) and translation
    transformed = np.dot(coords, rot_matrix.T) + trans_vector
    
    return transformed


def apply_alignment_transformation(coordinate_df: pd.DataFrame, rotation_matrix: np.ndarray, 
                                  translation_vector: np.ndarray) -> pd.DataFrame:
    """
    Apply rotation and translation to coordinate DataFrame.
    
    Args:
        coordinate_df: DataFrame with x, y, z columns
        rotation_matrix: 3x3 rotation matrix
        translation_vector: Translation vector of length 3
        
    Returns:
        DataFrame with transformed coordinates
    """
    # Extract coordinates and convert to numpy array
    coords = coordinate_df[['x', 'y', 'z']].values.astype(float)
    
    # Apply transformation
    transformed_coords = transform_coordinates(coords, rotation_matrix, translation_vector)
    
    # Create new DataFrame with transformed coordinates
    result_df = coordinate_df.copy()
    result_df[['x', 'y', 'z']] = transformed_coords
    
    return result_df


class GTalign:
    """Python interface for GTalign protein structure alignment tool."""

    def __init__(self, executable_path: str = "gtalign"):
        """
        Initialize GTalign interface.

        Args:
            executable_path: Path to the GTalign executable. Defaults to "gtalign",
                             which assumes it's in your PATH.
        """
        self.executable = executable_path
        self._check_executable()
        # Use a persistent output directory for GTalign operations.
        self.default_output_dir = os.path.join("data", "tempcif", "gtalign_output")
        os.makedirs(self.default_output_dir, exist_ok=True)

    def _check_executable(self):
        """Check if GTalign is available and executable."""
        try:
            result = subprocess.run([self.executable],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True,
                                    check=False)
            return True
        except FileNotFoundError:
            raise FileNotFoundError(
                f"GTalign executable not found at {self.executable}. Make sure it's installed and in your PATH."
            )

    def align(self,
              query_structures: Union[str, List[str]],
              reference_structures: Union[str, List[str]],
              output_dir: Optional[str] = None,
              tm_score_threshold: float = 0.5,
              speed: int = 9,
              depth: int = 2,
              sort_by: int = 0,
              max_hits: int = 100,
              referenced: bool = False,
              include_2tm_score: bool = False,
              verbose: bool = True,
              **kwargs) -> Dict:
        """
        Align query structures against reference structures.

        Args:
            query_structures: Path(s) to query structure file(s) or directory.
            reference_structures: Path(s) to reference structure file(s) or directory.
            output_dir: Directory to store output files. If None, a persistent subfolder is used.
            tm_score_threshold: Minimum TM-score to report (0-1).
            speed: Speed setting (0-13, higher is faster but potentially less accurate).
            depth: Search depth (0-3, lower is deeper).
            sort_by: Sorting criterion (0-8, see GTalign help).
            max_hits: Maximum number of hits to return.
            referenced: Generate transformation matrices for reference structures.
            include_2tm_score: Include secondary TM-score calculation.
            verbose: Print verbose output.
            **kwargs: Additional GTalign parameters.

        Returns:
            Dictionary with alignment results.
        """
        if output_dir is None:
            output_dir = os.path.join(self.default_output_dir, "alignment_output")
        os.makedirs(output_dir, exist_ok=True)

        if isinstance(query_structures, str):
            query_structures = [query_structures]
        if isinstance(reference_structures, str):
            reference_structures = [reference_structures]

        cmd = [self.executable]
        cmd.extend(["--qrs=" + ",".join(query_structures)])
        cmd.extend(["--rfs=" + ",".join(reference_structures)])
        cmd.extend(["-o", output_dir])
        cmd.extend(["-s", str(tm_score_threshold)])
        cmd.extend(["--speed=" + str(speed)])
        cmd.extend(["--depth=" + str(depth)])
        cmd.extend(["--sort=" + str(sort_by)])
        cmd.extend(["--outfmt=" + str(1)])
        cmd.extend(["--nhits=" + str(max_hits)])

        if referenced:
            cmd.append("--referenced")
        if include_2tm_score:
            cmd.append("--2tm-score")
        if verbose:
            cmd.append("-v")

        # Filter out unsupported options (save_aligned, save_transform, save_matrix)
        # Only pass known parameters to avoid errors
        supported_params = {
            "referenced", "2tm_score", "speed", "depth", "sort", "nhits", "query_format", 
            "reference_format", "chain", "outfmt"
        }
        
        for key, value in kwargs.items():
            if key in supported_params:
                param = key.replace('_', '-')
                if isinstance(value, bool):
                    if value:
                        cmd.append(f"--{param}")
                else:
                    cmd.append(f"--{param}={value}")

        result = subprocess.run(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)

        if result.returncode != 0:
            raise RuntimeError(f"GTalign failed with error:\n{result.stderr}")

        alignments = self._parse_results(output_dir)
        return {
            "command": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "alignments": alignments
        }

    def _parse_results(self, output_dir: str) -> List[Dict]:
        """Parse alignment results from output directory."""
        results = []
        
        # Get current time to use as a cutoff
        current_time = time.time()
        # Only consider files created in the last 60 seconds (from this run)
        cutoff_time = current_time - 60
        
        for root, _, files in os.walk(output_dir):
            for file in files:
                # Accept files ending in .aln or .out
                if file.endswith('.aln') or file.endswith('.out'):
                    alignment_file = os.path.join(root, file)
                    
                    # Check if file was recently created (in this run)
                    file_path = Path(alignment_file)
                    file_mod_time = file_path.stat().st_mtime
                    
                    if file_mod_time >= cutoff_time:
                        alignment_data = self._parse_alignment_file(alignment_file)
                        results.append(alignment_data)
        
        return results

    def _parse_alignment_file(self, alignment_file: str) -> Dict:
        """Parse a single alignment file from GTalign's output.

        This parser is designed for the GPU version output.
        It processes the header (extracting TM-score, RMSD, and aligned residues)
        and then collects alignment blocks (lines starting with 'Query:' and 'Refn:' or 'Reference:')
        to build the aligned query and reference sequences.
        It also extracts rotation matrix and translation vector if available.
        Returns a dictionary with the parsed alignment data.

        Args:
            alignment_file: Path to the GTalign output file

        Returns:
            Dictionary containing parsed alignment data including:
            - tm_score: TM-score of the alignment
            - rmsd: RMSD of the alignment
            - rotation_matrix: 3x3 rotation matrix as list of lists
            - translation_vector: 3-element translation vector
            - query_sequence: Aligned query sequence with gaps
            - reference_sequence: Aligned reference sequence with gaps
            - match_line: Alignment match symbols
            - alignment_blocks: Detailed information for each alignment block
            - alignment_path: Position mapping between query and reference [(q_pos, r_pos),...]
        """
        with open(alignment_file, 'r') as f:
            content = f.read()

        alignment_data = {
            "file": alignment_file,
            "query": "",
            "reference": "",
            "tm_score": 0.0,
            "rmsd": 0.0,
            "alignment_length": 0,
            "query_length": 0,
            "reference_length": 0,
            "query_sequence": "",  # Raw query sequence with gaps
            "reference_sequence": "",  # Raw reference sequence with gaps
            "query_sequence_clean": "",  # Cleaned query sequence (only amino acids and gaps)
            "reference_sequence_clean": "",  # Cleaned reference sequence (only amino acids and gaps)
            "match_line": "",  # Alignment match line
            "query_ss": "",  # Query secondary structure
            "ref_ss": "",  # Reference secondary structure
            "query_positions": [],  # List of (start, end) positions for each block
            "ref_positions": [],  # List of (start, end) positions for each block
            "alignment_blocks": [],  # Detailed info for each alignment block
            "alignment_path": [],  # Position mapping between query and reference [(q_pos, r_pos),...]
            "rotation_matrix": None,  # Will be 3x3 matrix as list of lists
            "translation_vector": None,  # Will be 3-element list
            "aligned_query_structure": None,  # Path to aligned query structure if available
            "aligned_reference_structure": None  # Path to aligned reference structure if available
        }

        # Extract query and reference structure paths
        query_path_pattern = r"Query\s+structure\s*:\s*(\S+)"
        ref_path_pattern = r"Refn\.\s+structure\s*:\s*(\S+)"
        
        query_path_match = re.search(query_path_pattern, content)
        ref_path_match = re.search(ref_path_pattern, content)
        
        if query_path_match:
            alignment_data["query"] = query_path_match.group(1)
        if ref_path_match:
            alignment_data["reference"] = ref_path_match.group(1)

        # Extract rotation matrix and translation vector
        # First find the section header for either Query or Reference
        rotation_section_pattern = r"Rotation\s*\[3,3\]\s*and\s*translation\s*\[3,1\]\s*for\s*(Query|Reference):"
        rotation_section_match = re.search(rotation_section_pattern, content, re.IGNORECASE)

        if rotation_section_match:
            # Get content after the header
            start_pos = rotation_section_match.end()
            remaining_content = content[start_pos:]

            # Look for lines that start with whitespace and contain numeric data
            matrix_lines = []
            for line in remaining_content.split('\n')[:10]:  # Check first 10 lines
                # Check if line starts with whitespace and contains numeric data
                if line.startswith(' ') and re.search(r'[-+]?\d+\.\d+', line):
                    matrix_lines.append(line)
                    if len(matrix_lines) == 3:  # We found all 3 matrix lines
                        break

            if len(matrix_lines) == 3:
                rotation_matrix = []
                translation_vector = []

                for line in matrix_lines:
                    # Find all numbers in the line
                    numbers = re.findall(r'[-+]?\d+\.\d+', line)
                    if len(numbers) >= 4:
                        rotation_matrix.append([float(numbers[0]), float(numbers[1]), float(numbers[2])])
                        translation_vector.append(float(numbers[3]))

                if len(rotation_matrix) == 3 and len(translation_vector) == 3:
                    alignment_data["rotation_matrix"] = rotation_matrix
                    alignment_data["translation_vector"] = translation_vector

        # Extract aligned structure paths if they exist
        aligned_query_pattern = r"Aligned\s+query\s+structure\s+saved\s+to:\s+(\S+)"
        aligned_ref_pattern = r"Aligned\s+reference\s+structure\s+saved\s+to:\s+(\S+)"

        query_match = re.search(aligned_query_pattern, content)
        ref_match = re.search(aligned_ref_pattern, content)

        if query_match:
            alignment_data["aligned_query_structure"] = query_match.group(1)
        if ref_match:
            alignment_data["aligned_reference_structure"] = ref_match.group(1)

        # Extract TM-score and RMSD
        tm_score_pattern = r"TM-score\s+\(Refn\./Query\)\s*=\s*([\d.]+)"
        rmsd_pattern = r"RMSD\s*=\s*([\d.]+)\s*A"

        tm_match = re.search(tm_score_pattern, content)
        rmsd_match = re.search(rmsd_pattern, content)

        if tm_match:
            alignment_data["tm_score"] = float(tm_match.group(1))
        if rmsd_match:
            alignment_data["rmsd"] = float(rmsd_match.group(1))

        # Extract alignment length
        alignment_length_pattern = r"Identities\s*=\s*\d+/(\d+)"
        alignment_length_match = re.search(alignment_length_pattern, content)
        if alignment_length_match:
            alignment_data["alignment_length"] = int(alignment_length_match.group(1))

        # Extract length information
        length_pattern = r"Length:\s*Refn\.\s*=\s*(\d+),\s*Query\s*=\s*(\d+)"
        length_match = re.search(length_pattern, content)
        if length_match:
            alignment_data["reference_length"] = int(length_match.group(1))
            alignment_data["query_length"] = int(length_match.group(2))

        # Extract sequence alignment - handle 5-line blocks
        alignment_blocks = []
        query_seq = ""
        ref_seq = ""
        match_line = ""
        query_ss = ""
        ref_ss = ""
        query_positions = []
        ref_positions = []

        lines = content.split("\n")
        i = 0

        # Find the start of the alignment section
        while i < len(lines):
            if "Length: Refn. =" in lines[i]:
                # We found the start of alignment section
                i += 1
                break
            i += 1

        # Skip past the TM-score and Identities lines
        while i < len(lines) and not lines[i].strip().startswith("struct"):
            i += 1

        # Process 5-line alignment blocks
        while i < len(lines):
            # Skip empty lines
            if not lines[i].strip():
                i += 1
                continue

            # If we've reached the rotation matrix section, we're done
            if "Rotation" in lines[i]:
                break

            # Look for a complete block that starts with "struct"
            if i + 4 < len(lines) and lines[i].strip().startswith("struct"):
                # Line 1: Query secondary structure
                query_ss_line = lines[i].strip()

                # Line 2: Query sequence
                query_line = lines[i + 1].strip()
                query_match = re.match(r'Query:\s+(\d+)\s+([A-Z\-]+)\s+(\d+)', query_line)

                # Line 3: Match line showing alignment quality
                match_line_text = lines[i + 2].strip()

                # Line 4: Reference sequence
                ref_line = lines[i + 3].strip()
                ref_match = re.match(r'(?:Refn\.|Reference):\s+(\d+)\s+([A-Z\-]+)\s+(\d+)', ref_line)

                # Line 5: Reference secondary structure
                ref_ss_line = lines[i + 4].strip()

                if query_match and ref_match:
                    # Extract positions
                    query_start = int(query_match.group(1))
                    query_end = int(query_match.group(3))
                    ref_start = int(ref_match.group(1))
                    ref_end = int(ref_match.group(3))

                    # Extract sequences
                    query_seq_part = query_match.group(2)
                    ref_seq_part = ref_match.group(2)

                    # Extract match line - start from the same position as query seq part
                    match_part = match_line_text[match_line_text.find('+'):]
                    if len(match_part) > len(query_seq_part):
                        match_part = match_part[:len(query_seq_part)]
                    elif len(match_part) < len(query_seq_part):
                        match_part = match_part.ljust(len(query_seq_part))

                    # Extract secondary structure annotations
                    # Find the position where the sequence part starts
                    query_ss_offset = len("struct")
                    query_ss_part = query_ss_line[query_ss_offset:].strip()
                    if len(query_ss_part) > len(query_seq_part):
                        query_ss_part = query_ss_part[:len(query_seq_part)]
                    elif len(query_ss_part) < len(query_seq_part):
                        query_ss_part = query_ss_part.ljust(len(query_seq_part))

                    ref_ss_offset = len("struct")
                    ref_ss_part = ref_ss_line[ref_ss_offset:].strip()
                    if len(ref_ss_part) > len(ref_seq_part):
                        ref_ss_part = ref_ss_part[:len(ref_seq_part)]
                    elif len(ref_ss_part) < len(ref_seq_part):
                        ref_ss_part = ref_ss_part.ljust(len(ref_seq_part))

                    # Store the block information
                    alignment_blocks.append({
                        'query_start': query_start,
                        'query_end': query_end,
                        'ref_start': ref_start,
                        'ref_end': ref_end,
                        'query_seq': query_seq_part,
                        'ref_seq': ref_seq_part,
                        'match_line': match_part,
                        'query_ss': query_ss_part,
                        'ref_ss': ref_ss_part
                    })

                    # Append to full sequences
                    query_seq += query_seq_part
                    ref_seq += ref_seq_part
                    match_line += match_part
                    query_ss += query_ss_part
                    ref_ss += ref_ss_part

                    # Store start/end positions
                    query_positions.append((query_start, query_end))
                    ref_positions.append((ref_start, ref_end))

                    # Move to the next potential block
                    i += 5
                else:
                    # If pattern doesn't match, move to next line
                    i += 1
            else:
                # Not a valid block start, move to next line
                i += 1

        # Store the raw sequences
        alignment_data["query_sequence"] = query_seq
        alignment_data["reference_sequence"] = ref_seq
        alignment_data["match_line"] = match_line

        # Also store cleaned sequences (only standard amino acids and gaps)
        alignment_data["query_sequence_clean"] = re.sub(r'[^A-Z\-]', '', query_seq)
        alignment_data["reference_sequence_clean"] = re.sub(r'[^A-Z\-]', '', ref_seq)

        # Add additional information
        alignment_data["query_ss"] = query_ss  # Query secondary structure
        alignment_data["ref_ss"] = ref_ss  # Reference secondary structure
        alignment_data["query_positions"] = query_positions  # List of (start, end) positions
        alignment_data["ref_positions"] = ref_positions  # List of (start, end) positions
        alignment_data["alignment_blocks"] = alignment_blocks  # List of block dictionaries

        # Generate alignment path (position mapping)
        alignment_path = []

        # Process each alignment block to build the position mapping
        for block in alignment_blocks:
            query_start = block["query_start"]
            ref_start = block["ref_start"]
            query_seq_part = block["query_seq"]
            ref_seq_part = block["ref_seq"]

            # Create position mapping for this block
            q_pos = query_start
            r_pos = ref_start

            for q_char, r_char in zip(query_seq_part, ref_seq_part):
                if q_char == '-' and r_char == '-':
                    # Both are gaps (shouldn't happen in valid alignments)
                    continue
                elif q_char == '-':
                    # Gap in query
                    alignment_path.append((None, r_pos))
                    r_pos += 1
                elif r_char == '-':
                    # Gap in reference
                    alignment_path.append((q_pos, None))
                    q_pos += 1
                else:
                    # Match or mismatch
                    alignment_path.append((q_pos, r_pos))
                    q_pos += 1
                    r_pos += 1

        # Store the alignment path
        alignment_data["alignment_path"] = alignment_path

        return alignment_data

    def cluster(self,
                structures: Union[str, List[str]],
                output_dir: Optional[str] = None,
                tm_score_threshold: float = 0.5,
                coverage: float = 0.7,
                speed: int = 13,
                algorithm: int = 0,
                verbose: bool = True,
                **kwargs) -> Dict:
        """
        Cluster protein structures.

        Args:
            structures: Path(s) to structure file(s) or directory.
            output_dir: Directory to store output files. If None, a persistent subfolder is used.
            tm_score_threshold: TM-score threshold for clustering (0-1).
            coverage: Length coverage threshold (0-1).
            speed: Speed setting (0-13, higher is faster).
            algorithm: Clustering algorithm (0: complete-linkage, 1: single-linkage).
            verbose: Print verbose output.
            **kwargs: Additional GTalign parameters.

        Returns:
            Dictionary with clustering results.
        """
        if output_dir is None:
            output_dir = os.path.join(self.default_output_dir, "cluster_output")
        os.makedirs(output_dir, exist_ok=True)
        if isinstance(structures, str):
            structures = [structures]
        for struct in structures:
            if not os.path.exists(struct):
                raise FileNotFoundError(f"Structure file or directory not found: {struct}")
        cmd = [self.executable]
        cmd.extend(["--cls=" + ",".join(structures)])
        cmd.extend(["-o", output_dir])
        cmd.extend(["--cls-threshold=" + str(tm_score_threshold)])
        cmd.extend(["--cls-coverage=" + str(coverage)])
        cmd.extend(["--speed=" + str(speed)])

        cmd.extend(["--cls-algorithm=" + str(algorithm)])
        if verbose:
            cmd.append("-v")
        for key, value in kwargs.items():
            param = key.replace('_', '-')
            if isinstance(value, bool):
                if value:
                    cmd.append(f"--{param}")
            else:
                cmd.append(f"--{param}={value}")
        result = subprocess.run(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
        if result.returncode != 0:
            raise RuntimeError(f"GTalign clustering failed with error:\n{result.stderr}")
        clusters = self._parse_clustering_results(output_dir)
        return {
            "command": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "clusters": clusters
        }

    def _parse_clustering_results(self, output_dir: str) -> List[Dict]:
        """Parse clustering results from output directory."""
        clusters = []
        cluster_file = os.path.join(output_dir, "clusters.txt")
        if os.path.exists(cluster_file):
            with open(cluster_file, 'r') as f:
                content = f.read()
            current_cluster = None
            for line in content.split('\n'):
                if line.startswith("Cluster"):
                    if current_cluster:
                        clusters.append(current_cluster)
                    try:
                        cluster_id = int(line.split()[1])
                        current_cluster = {"id": cluster_id, "members": []}
                    except (IndexError, ValueError):
                        pass
                elif line.strip() and current_cluster:
                    parts = line.strip().split()
                    if parts:
                        current_cluster["members"].append(parts[0])
            if current_cluster:
                clusters.append(current_cluster)
        return clusters

    def get_temp_dir(self) -> str:
        """Return the default persistent temp directory for GTalign operations."""
        return self.default_output_dir

    def save_structure_as_cif(self, structure_df: pd.DataFrame, output_path: str) -> str:
        """
        Save a transformed structure DataFrame as a CIF file.
        
        Args:
            structure_df: DataFrame containing transformed structure data
            output_path: Path to save the CIF file
            
        Returns:
            Path to the saved CIF file
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Write mmCIF format header
        with open(output_path, 'w') as f:
            f.write("data_aligned_structure\n")
            f.write("#\n")
            f.write("loop_\n")
            f.write("_atom_site.group_PDB\n")
            f.write("_atom_site.id\n")
            f.write("_atom_site.type_symbol\n")
            f.write("_atom_site.label_atom_id\n")
            f.write("_atom_site.label_alt_id\n")
            f.write("_atom_site.label_comp_id\n")
            f.write("_atom_site.label_asym_id\n")
            f.write("_atom_site.label_entity_id\n")
            f.write("_atom_site.label_seq_id\n")
            f.write("_atom_site.pdbx_PDB_ins_code\n")
            f.write("_atom_site.Cartn_x\n")
            f.write("_atom_site.Cartn_y\n")
            f.write("_atom_site.Cartn_z\n")
            f.write("_atom_site.occupancy\n")
            f.write("_atom_site.B_iso_or_equiv\n")
            f.write("_atom_site.pdbx_formal_charge\n")
            f.write("_atom_site.auth_seq_id\n")
            f.write("_atom_site.auth_comp_id\n")
            f.write("_atom_site.auth_asym_id\n")
            f.write("_atom_site.auth_atom_id\n")
            f.write("_atom_site.pdbx_PDB_model_num\n")
            
            # Write atom data
            col_map = {
                'group_PDB': 'record_name',
                'id': 'atom_id',
                'type_symbol': 'element_symbol',
                'label_atom_id': 'atom_name',
                'label_alt_id': 'alt_loc',
                'label_comp_id': 'residue_name',
                'label_asym_id': 'chain_id',
                'label_entity_id': 'entity_id',
                'label_seq_id': 'residue_number',
                'pdbx_PDB_ins_code': 'insertion_code',
                'Cartn_x': 'x',
                'Cartn_y': 'y',
                'Cartn_z': 'z',
                'occupancy': 'occupancy',
                'B_iso_or_equiv': 'b_factor',
                'pdbx_formal_charge': 'charge',
                'auth_seq_id': 'residue_number',
                'auth_comp_id': 'residue_name',
                'auth_asym_id': 'chain_id',
                'auth_atom_id': 'atom_name',
                'pdbx_PDB_model_num': 'model_num'
            }
            
            # Sort the DataFrame by chain, residue number, and atom name for clean output
            if 'residue_number' in structure_df.columns and 'chain_id' in structure_df.columns:
                structure_df = structure_df.sort_values(['chain_id', 'residue_number', 'atom_name'])
            
            for _, row in structure_df.iterrows():
                # Prepare a list of values for this atom
                values = []
                for cif_col, df_col in col_map.items():
                    if df_col in row.index:
                        # Format differently based on column
                        if df_col in ['x', 'y', 'z']:
                            # Format coordinates with 3 decimal places
                            values.append(f"{row[df_col]:.3f}")
                        elif df_col in ['occupancy', 'b_factor']:
                            # Format B-factor and occupancy with 2 decimal places
                            values.append(f"{row[df_col]:.2f}")
                        else:
                            values.append(str(row[df_col]))
                    else:
                        # Use placeholders for missing columns
                        if cif_col in ['Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv']:
                            values.append("0.00")
                        elif cif_col in ['pdbx_formal_charge', 'label_entity_id', 'pdbx_PDB_model_num']:
                            values.append("1")
                        elif cif_col == 'group_PDB' and 'record_name' not in row.index:
                            values.append("ATOM")
                        else:
                            values.append(".")
                            
                # Write the atom line
                f.write(" ".join(values) + "\n")
                
        return output_path
    
    def align_from_processor(self, 
                            cif_processor, 
                            reference_pdb_id: str, 
                            target_pdb_ids: List[str],
                            chain_id: str = 'A', 
                            output_dir: Optional[str] = None,
                            tm_score_threshold: float = 0.5,
                            speed: int = 9,
                            referenced: bool = True,
                            save_transformed: bool = False,
                            save_matrix: bool = False,
                            **kwargs) -> Dict:
        """
        Comprehensive helper to align structures using GTalign directly from CifProcessor.
        
        Args:
            cif_processor: CifProcessor instance with loaded structures
            reference_pdb_id: PDB ID of reference structure
            target_pdb_ids: List of PDB IDs to align
            chain_id: Chain ID to use (default: 'A')
            output_dir: Output directory
            tm_score_threshold: TM-score threshold
            speed: Speed setting (0-13, higher is faster)
            referenced: Generate transformation matrices for reference structures
            save_transformed: Whether to save transformed structures as CIF files
            save_matrix: Whether to save transformation matrices
            **kwargs: Additional parameters for GTalign.align
            
        Returns:
            Dictionary with alignment results and transformed structures
        """
        # Create a unique timestamp for this run
        run_timestamp = int(time.time())
        
        # Create temp directory for intermediate files with unique timestamp
        temp_dir = tempfile.mkdtemp() if output_dir is None else os.path.join(output_dir, f'temp_{run_timestamp}')
        os.makedirs(temp_dir, exist_ok=True)
        
        # Save reference structure
        ref_file = cif_processor.save_temp_cif(reference_pdb_id, suffix=f"chain{chain_id}", output_dir=temp_dir)
        
        # Save target structures
        target_files = []
        target_file_map = {}  # Map PDB IDs to their temporary filenames
        for pdb_id in target_pdb_ids:
            target_file = cif_processor.save_temp_cif(pdb_id, suffix=f"chain{chain_id}", output_dir=temp_dir)
            target_files.append(str(target_file))
            target_file_map[pdb_id] = str(target_file)
        
        # Create a dedicated output directory for this run
        align_output_dir = os.path.join(temp_dir, f'align_output_{run_timestamp}')
        os.makedirs(align_output_dir, exist_ok=True)
        
        # Create directory for transformed structures if needed
        transformed_dir = os.path.join(align_output_dir, "transformed_structures")
        if save_transformed:
            os.makedirs(transformed_dir, exist_ok=True)
            
        # Create directory for transformation matrices if needed
        matrix_dir = os.path.join(align_output_dir, "matrices")
        if save_matrix:
            os.makedirs(matrix_dir, exist_ok=True)
        
        print(f"Saving reference and target structures to: {temp_dir}")
        print(f"Saving alignment output to: {align_output_dir}")
        
        # Run alignment
        result = self.align(
            query_structures=str(ref_file),
            reference_structures=target_files,
            output_dir=align_output_dir,
            tm_score_threshold=tm_score_threshold,
            speed=speed,
            referenced=referenced,
            **kwargs
        )
        
        # Debug information
        print(f"Target file map: {target_file_map}")
        print(f"Alignments found: {len(result['alignments'])}")
        
        # Apply transformations to structures
        aligned_structures = {}
        transformed_structure_paths = {}
        
        # Get the first alignment with rotation matrix and translation vector
        for alignment in result['alignments']:
            if alignment.get('rotation_matrix') is not None and alignment.get('translation_vector') is not None:
                rot_matrix = alignment.get('rotation_matrix')
                trans_vector = alignment.get('translation_vector')
                print(f"Using transformation from alignment: {os.path.basename(alignment['file'])}")
                
                # Apply this transformation to all target structures
                for pdb_id in target_pdb_ids:
                    # Get the structure data
                    structure = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id].copy()
                    
                    # Force coordinate values to float for transformation
                    structure[['x', 'y', 'z']] = structure[['x', 'y', 'z']].astype(float)
                    
                    # Apply transformation
                    transformed_structure = apply_alignment_transformation(
                        structure, rot_matrix, trans_vector
                    )
                    
                    # Store aligned structure
                    aligned_structures[pdb_id] = transformed_structure
                    
                    # Save transformed structure if requested
                    if save_transformed:
                        transformed_file_path = os.path.join(
                            transformed_dir, f"{pdb_id}_aligned_{run_timestamp}.cif"
                        )
                        self.save_structure_as_cif(transformed_structure, transformed_file_path)
                        transformed_structure_paths[pdb_id] = transformed_file_path
                        print(f"Saved transformed structure for {pdb_id} to {transformed_file_path}")
                
                # Save transformation matrix if requested
                if save_matrix:
                    matrix_file_path = os.path.join(matrix_dir, f"transformation_matrix_{run_timestamp}.txt")
                    with open(matrix_file_path, 'w') as f:
                        f.write("# Transformation Matrix\n")
                        f.write("# Rotation matrix [3,3] and translation vector [3,1]\n")
                        f.write("# Format: rotation matrix (3x3) followed by translation vector (3x1)\n")
                        f.write("\n")
                        for i in range(3):
                            rot_str = " ".join([f"{rot_matrix[i][j]:10.6f}" for j in range(3)])
                            f.write(f"{rot_str}  {trans_vector[i]:10.6f}\n")
                    print(f"Saved transformation matrix to {matrix_file_path}")
                
                break
        
        # Add aligned structures and paths to result
        result['aligned_structures'] = aligned_structures
        
        if save_transformed:
            result['transformed_structure_paths'] = transformed_structure_paths
        
        return result

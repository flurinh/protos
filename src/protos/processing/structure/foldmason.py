import subprocess
import shlex
import os
import shutil
import re
import shlex
import platform
from pathlib import Path
from typing import List, Dict, Optional, Any, Union, Sequence


def is_windows_path(path: Union[str, Path]) -> bool:
    """Check if a path is a Windows-style path."""
    path_str = str(path)
    # Check for drive letter pattern (C:, D:, etc.)
    if len(path_str) > 1 and path_str[1] == ':':
        return True
    # Check for UNC paths (\\server\share)
    if path_str.startswith('\\\\'):
        return True
    # Check for backslashes in the path
    if '\\' in path_str:
        return True
    return False


def to_wsl_path(win_path):
    """Convert a Windows path to WSL path."""
    import re

    # Check if the path starts with a drive letter (e.g., 'C:')
    match = re.match(r'^([A-Za-z]:)', win_path)

    if match:
        # Extract the drive letter (e.g., 'C:')
        drive = match.group(1)
        # Get the rest of the path after the drive letter
        rest = win_path[len(drive):]
        # Split the rest by any sequence of \ or /, filtering out empty strings
        components = [comp for comp in re.split(r'[\\/]+', rest) if comp]
        # Convert drive to WSL format (e.g., 'C:' -> '/mnt/c/') and join components
        wsl_path = '/mnt/' + drive[0].lower() + '/' + '/'.join(components)
    else:
        # No drive letter, split the entire path by \ or /, filtering out empty strings
        components = [comp for comp in re.split(r'[\\/]+', win_path) if comp]
        # Join components with forward slashes
        wsl_path = '/' + '/'.join(components)

    # Clean up any double slashes
    while '//' in wsl_path:
        wsl_path = wsl_path.replace('//', '/')

    return wsl_path


class FoldMason:
    """Python interface for FoldMason protein structure alignment and analysis tool."""

    def __init__(self, use_wsl: bool = True, debug: bool = True):
        """
        Initialize FoldMason interface.

        Args:
            use_wsl: If True, commands will be run via WSL.
            debug: If True, print debug information.
        """
        self.use_wsl = use_wsl
        self.debug = debug

        # Detect if we're on Windows
        self.is_windows = platform.system() == 'Windows'

    def run_command(self, command, debug=False):
        """Run a command via WSL if on Windows, with proper path conversion."""
        if self.debug or debug:
            print(f"DEBUG: Original command: {command}")

        if not self.use_wsl or not self.is_windows:
            # If not using WSL or not on Windows, just run the command directly
            result = subprocess.run(
                command, shell=True, capture_output=True, text=True, check=True
            )
            return result.stdout.strip()

        # Simple approach: just split by spaces
        parts = command.split()
        converted_parts = []

        for part in parts:
            # Check if this part looks like a Windows path
            # Simple heuristic: contains a colon and backslash
            if ':' in part and ('\\' in part or part.startswith('\\\\')):
                try:
                    # Convert Windows path to WSL path
                    wsl_path = to_wsl_path(part)
                    converted_parts.append(wsl_path)
                    if self.debug or debug:
                        print(f"DEBUG: Converted path: {part} -> {wsl_path}")
                except Exception as e:
                    if self.debug or debug:
                        print(f"DEBUG: Path conversion failed: {e}")
                    raise ValueError(f"Failed to convert path '{part}': {e}")
            else:
                # Not a path, keep as is
                converted_parts.append(part)

        # Reconstruct the command with converted paths
        wsl_command = " ".join(shlex.quote(part) for part in converted_parts)

        # Create the full WSL command
        full_command = f"wsl bash -lic {shlex.quote(wsl_command)}"

        if self.debug or debug:
            print(f"DEBUG: Running WSL command: {full_command}")

        try:
            result = subprocess.run(
                full_command, shell=True, capture_output=True, text=True, check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError as e:
            if self.debug or debug:
                print(f"DEBUG: Command failed: {e.stderr.strip()}")
            raise RuntimeError(f"FoldMason command failed: {e.stderr.strip()}")

    def easy_msa(self,
                 input_files: Sequence[Union[str, Path]],
                 output_prefix: Union[str, Path],
                 tmp_folder: Union[str, Path],
                 report_mode: int = 1,
                 precluster: bool = False) -> Dict[str, Any]:
        """
        Run the 'easy-msa' command to generate a multiple sequence alignment from structure files.

        Args:
            input_files: List of input structure files (PDB or mmCIF)
            output_prefix: Output file prefix
            tmp_folder: Directory for temporary files
            report_mode: Report mode (0, 1, or 2)
            precluster: Whether to pre-cluster structures before MSA

        Returns:
            Dictionary with command output and file paths
        """
        # Convert paths to Path objects and then to strings
        input_paths = [str(Path(f).resolve()) for f in input_files]
        output_prefix_path = Path(output_prefix).resolve()
        tmp_folder_path = Path(tmp_folder).resolve()

        # Ensure tmp folder exists
        tmp_folder_path.mkdir(parents=True, exist_ok=True)

        # Build command with positional arguments
        input_files_str = ' '.join(input_paths)
        command = f"foldmason easy-msa {input_files_str} {output_prefix_path} {tmp_folder_path}"

        # Add options
        if report_mode != 0:
            command += f" --report-mode {report_mode}"
        if precluster:
            command += " --precluster"

        if self.debug:
            print(f"DEBUG: Running easy-msa command")
            print(f"DEBUG: Input files: {input_paths}")
            print(f"DEBUG: Output prefix: {output_prefix_path}")
            print(f"DEBUG: Temp folder: {tmp_folder_path}")

        output = self.run_command(command)

        # Expected output files based on report_mode
        expected_files = {
            "alignment": f"{output_prefix_path}_alignment.fasta"
        }

        if report_mode == 1:
            expected_files["report"] = f"{output_prefix_path}_report.html"
        elif report_mode == 2:
            expected_files["report"] = f"{output_prefix_path}_report.json"

        return {
            "command": command,
            "stdout": output,
            "files": expected_files
        }

    def convertalis(self,
                    input_db: Union[str, Path],
                    output_format: str,
                    output_file: Union[str, Path],
                    extra_args: str = "") -> Dict[str, Any]:
        """
        Run the 'convertalis' command to convert alignment database to different formats.

        Args:
            input_db: Path to the alignment DB
            output_format: Desired output format (e.g., 'BLAST-tab', 'SAM')
            output_file: Path for the converted output
            extra_args: Any extra command-line arguments

        Returns:
            Dictionary with command output and file paths
        """
        input_db_path = Path(input_db).resolve()
        output_file_path = Path(output_file).resolve()

        cmd = f"foldmason convertalis {input_db_path} {output_format} {output_file_path} {extra_args}".strip()
        output = self.run_command(cmd)

        return {
            "command": cmd,
            "stdout": output,
            "output_file": str(output_file_path)
        }

    def structuremsa(self,
                     structure_db: Union[str, Path],
                     output_prefix: Union[str, Path]) -> Dict[str, Any]:
        """
        Run the 'structuremsa' command to create a structure-based MSA.

        Args:
            structure_db: Path to the structure database
            output_prefix: Output file prefix

        Returns:
            Dictionary with command output and file paths
        """
        structure_db_path = Path(structure_db).resolve()
        output_prefix_path = Path(output_prefix).resolve()

        cmd = f"foldmason structuremsa {structure_db_path} {output_prefix_path}"
        output = self.run_command(cmd)

        expected_output_files = {
            "alignment": f"{output_prefix_path}.fasta",
            "tree": f"{output_prefix_path}.dnd"
        }

        return {
            "command": cmd,
            "stdout": output,
            "files": expected_output_files
        }

    def structuremsacluster(self,
                            structure_db: Union[str, Path],
                            output_prefix: Union[str, Path],
                            cluster_threshold: float = 0.5) -> Dict[str, Any]:
        """
        Run the 'structuremsacluster' command for structure-based clustering.

        Args:
            structure_db: Path to the structure database
            output_prefix: Output file prefix
            cluster_threshold: Clustering threshold (default 0.5)

        Returns:
            Dictionary with command output and file paths
        """
        structure_db_path = Path(structure_db).resolve()
        output_prefix_path = Path(output_prefix).resolve()

        cmd = f"foldmason structuremsacluster {structure_db_path} {output_prefix_path} --cluster-threshold {cluster_threshold}"
        output = self.run_command(cmd)

        expected_output_files = {
            "alignment": f"{output_prefix_path}.fasta",
            "tree": f"{output_prefix_path}.dnd",
            "clusters": f"{output_prefix_path}_clusters.tsv"
        }

        return {
            "command": cmd,
            "stdout": output,
            "files": expected_output_files
        }

    def msa2lddt(self,
                 structure_db: Union[str, Path],
                 input_fasta: Union[str, Path]) -> Dict[str, Any]:
        """
        Run the 'msa2lddt' command to calculate the LDDT score of an MSA.

        Args:
            structure_db: Path to the structure database
            input_fasta: Path to the input FASTA alignment file

        Returns:
            Dictionary with command output and LDDT score
        """
        structure_db_path = Path(structure_db).resolve()
        input_fasta_path = Path(input_fasta).resolve()

        cmd = f"foldmason msa2lddt {structure_db_path} {input_fasta_path}"
        output = self.run_command(cmd)

        # Extract LDDT score from output if present
        lddt_score = None
        for line in output.splitlines():
            if "LDDT score:" in line:
                try:
                    lddt_score = float(line.split("LDDT score:")[1].strip())
                except (ValueError, IndexError):
                    pass

        return {
            "command": cmd,
            "stdout": output,
            "lddt_score": lddt_score
        }

    def msa2lddtreport(self,
                       structure_db: Union[str, Path],
                       input_fasta: Union[str, Path],
                       output_html: Union[str, Path],
                       guide_tree: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
        """
        Run the 'msa2lddtreport' command to calculate LDDT and generate an HTML report.

        Args:
            structure_db: Path to the structure database
            input_fasta: Path to the input FASTA alignment file
            output_html: Path for the HTML report
            guide_tree: Path to the guide tree file (optional)

        Returns:
            Dictionary with command output and file paths
        """
        structure_db_path = Path(structure_db).resolve()
        input_fasta_path = Path(input_fasta).resolve()
        output_html_path = Path(output_html).resolve()

        if guide_tree:
            guide_tree_path = Path(guide_tree).resolve()
            cmd = f"foldmason msa2lddtreport {structure_db_path} {input_fasta_path} {output_html_path} --guide-tree {guide_tree_path}"
        else:
            cmd = f"foldmason msa2lddtreport {structure_db_path} {input_fasta_path} {output_html_path}"

        output = self.run_command(cmd)

        return {
            "command": cmd,
            "stdout": output,
            "output_file": str(output_html_path)
        }

    def msa2lddtjson(self,
                     structure_db: Union[str, Path],
                     input_fasta: Union[str, Path],
                     output_json: Union[str, Path]) -> Dict[str, Any]:
        """
        Run the 'msa2lddtjson' command to calculate LDDT and generate a JSON report.

        Args:
            structure_db: Path to the structure database
            input_fasta: Path to the input FASTA alignment file
            output_json: Path for the JSON report

        Returns:
            Dictionary with command output and file paths
        """
        structure_db_path = Path(structure_db).resolve()
        input_fasta_path = Path(input_fasta).resolve()
        output_json_path = Path(output_json).resolve()

        cmd = f"foldmason msa2lddtjson {structure_db_path} {input_fasta_path} {output_json_path}"
        output = self.run_command(cmd)

        return {
            "command": cmd,
            "stdout": output,
            "output_file": str(output_json_path)
        }

    def refinemsa(self,
                  structure_db: Union[str, Path],
                  input_fasta: Union[str, Path],
                  output_fasta: Union[str, Path],
                  refine_iters: int = 1000) -> Dict[str, Any]:
        """
        Run the 'refinemsa' command to iteratively refine an MSA.

        Args:
            structure_db: Path to the structure database
            input_fasta: Path to the input FASTA alignment file
            output_fasta: Path where the refined alignment will be saved
            refine_iters: Number of refinement iterations (default: 1000)

        Returns:
            Dictionary with command output and file paths
        """
        structure_db_path = Path(structure_db).resolve()
        input_fasta_path = Path(input_fasta).resolve()
        output_fasta_path = Path(output_fasta).resolve()

        cmd = f"foldmason refinemsa {structure_db_path} {input_fasta_path} {output_fasta_path} --refine-iters {refine_iters}"
        output = self.run_command(cmd)

        return {
            "command": cmd,
            "stdout": output,
            "output_file": str(output_fasta_path)
        }
"""
Input Manager for safe user data registration.

This module provides a safe workflow for registering user-provided data files,
replacing the unsafe "drag-and-drop" approach. Files are placed in a designated
input folder, validated, and then registered with appropriate conflict resolution.
"""

import os
import shutil
import hashlib
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any, Tuple
from enum import Enum
from dataclasses import dataclass
import logging

from protos.io.paths import ProtosPaths
from protos.io.entity_registry import EntityRegistry
from protos.io.format_registry import format_registry, ProcessorType
from protos.io.conflict_resolver import ConflictResolutionStrategy


@dataclass
class InputFile:
    """Represents a file to be processed from the input folder."""
    path: Path
    file_type: str
    entity_name: str
    size: int = 0
    modified: datetime = None


@dataclass
class ProcessingResult:
    """Result of processing a single input file."""
    input_file: InputFile
    success: bool
    action_taken: str
    entity_name: str = None
    error_message: str = None
    details: Dict[str, Any] = None


class ProcessingReport:
    """Aggregates results from processing multiple files."""
    
    def __init__(self):
        self.processed = []
        self.rejected = []
        self.errors = []
        self.skipped = []
        self.total = 0
    
    def add_result(self, input_file: InputFile, result: ProcessingResult):
        """Add a processing result."""
        self.total += 1
        if result.success:
            if result.action_taken == "skipped":
                self.skipped.append(result)
            else:
                self.processed.append(result)
        else:
            if result.error_message:
                self.errors.append(result)
            else:
                self.rejected.append(result)
    
    def add_rejected(self, input_file: InputFile, errors: List[str]):
        """Add a rejected file."""
        self.total += 1
        result = ProcessingResult(
            input_file=input_file,
            success=False,
            action_taken="rejected",
            error_message="; ".join(errors)
        )
        self.rejected.append(result)
    
    def add_error(self, input_file: InputFile, error: str):
        """Add an error result."""
        self.total += 1
        result = ProcessingResult(
            input_file=input_file,
            success=False,
            action_taken="error",
            error_message=error
        )
        self.errors.append(result)
    
    def display(self):
        """Display report summary."""
        print(f"\nProcessing Report:")
        print(f"Total files: {self.total}")
        print(f"Successfully processed: {len(self.processed)}")
        print(f"Skipped (duplicates): {len(self.skipped)}")
        print(f"Rejected (validation): {len(self.rejected)}")
        print(f"Errors: {len(self.errors)}")
        
        if self.processed:
            print("\nProcessed files:")
            for result in self.processed:
                print(f"  ✓ {result.input_file.path.name} -> {result.entity_name} ({result.action_taken})")
        
        if self.skipped:
            print("\nSkipped files (already registered):")
            for result in self.skipped:
                print(f"  - {result.input_file.path.name} (duplicate of {result.entity_name})")
        
        if self.rejected:
            print("\nRejected files:")
            for result in self.rejected:
                print(f"  ✗ {result.input_file.path.name}: {result.error_message}")
        
        if self.errors:
            print("\nErrors:")
            for result in self.errors:
                print(f"  ✗ {result.input_file.path.name}: {result.error_message}")


class InputManager:
    """Manages user data input and registration."""
    
    def __init__(self, paths: ProtosPaths = None, entity_registry: EntityRegistry = None):
        """Initialize the input manager."""
        self.paths = paths or ProtosPaths()
        self.entity_registry = entity_registry or EntityRegistry(paths=self.paths)
        
        # Set up directories
        self.input_dir = self._get_input_directory()
        self.processed_dir = self._get_processed_directory()
        self.rejected_dir = self._get_rejected_directory()
        
        # Set up logging
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Create README in input directory
        self._create_readme()
    
    def _get_input_directory(self) -> Path:
        """Get user input directory."""
        input_dir = Path(self.paths.data_root) / "input"
        input_dir.mkdir(exist_ok=True)
        return input_dir
    
    def _get_processed_directory(self) -> Path:
        """Get processed input directory."""
        processed = Path(self.paths.data_root) / "input" / "processed"
        processed.mkdir(parents=True, exist_ok=True)
        return processed
    
    def _get_rejected_directory(self) -> Path:
        """Get rejected input directory."""
        rejected = Path(self.paths.data_root) / "input" / "rejected"
        rejected.mkdir(parents=True, exist_ok=True)
        return rejected
    
    def _create_readme(self):
        """Create README file in input directory with instructions."""
        readme_path = self.input_dir / "README.txt"
        if not readme_path.exists():
            readme_content = """
PROTOS INPUT FOLDER
==================

Place your data files here for registration into Protos.

Supported file types:
- Structure files: .cif, .pdb, .mmcif
- Sequence files: .fasta, .fa, .faa
- Property tables: .csv, .tsv
- Embeddings: .npy, .pt

Usage:
1. Copy your files to this folder
2. Run: protos register-input --dry-run  (to preview)
3. Run: protos register-input  (to register)

Files will be:
- Validated for format and integrity
- Registered in the appropriate processor
- Moved to processed/ if successful
- Moved to rejected/ with error log if failed

For more information:
protos register-input --help
"""
            readme_path.write_text(readme_content)
    
    def _detect_file_type(self, file_path: Path) -> Optional[str]:
        """Detect the file type based on extension and content."""
        # Use format registry to detect format
        fmt = format_registry.get_format_by_path(file_path)
        
        if fmt:
            # Map processor type to file type string for backward compatibility
            processor_map = {
                ProcessorType.STRUCTURE: 'structure',
                ProcessorType.SEQUENCE: 'sequence',
                ProcessorType.GRN: 'grn',
                ProcessorType.PROPERTY: 'property',
                ProcessorType.EMBEDDING: 'embedding',
                ProcessorType.LIGAND: 'ligand',
                ProcessorType.GRAPH: 'graph'
            }
            file_type = processor_map.get(fmt.processor)
            
            # Special handling for ambiguous formats
            if file_path.suffix.lower() in ['.csv', '.tsv']:
                # Could be property or GRN table - check content
                content_type = self._detect_table_type(file_path)
                if content_type:
                    file_type = content_type
            
            return file_type
        else:
            # Try content-based detection for unknown extensions
            return self._detect_by_content(file_path)
    
    def _detect_table_type(self, file_path: Path) -> Optional[str]:
        """Detect if a table file is property or GRN format."""
        try:
            import pandas as pd
            # Read first few rows
            sep = '\t' if file_path.suffix.lower() == '.tsv' else ','
            df = pd.read_csv(file_path, sep=sep, nrows=5)
            
            # Check for GRN-like columns (e.g., "1.50", "2.50", etc.)
            grn_pattern = r'^\d+\.\d+$'
            grn_cols = [col for col in df.columns if isinstance(col, str) and 
                       pd.Series([col]).str.match(grn_pattern).any()]
            
            if len(grn_cols) > 3:  # Multiple GRN columns suggest GRN table
                return 'grn'
            else:
                return 'property'
        except:
            return None
    
    def _detect_by_content(self, file_path: Path) -> Optional[str]:
        """Detect file type by examining content."""
        try:
            # Read first few lines
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                first_lines = [f.readline() for _ in range(10)]
                content = ''.join(first_lines)
            
            # Check for common patterns
            if content.startswith('>'):
                return 'sequence'  # FASTA format
            elif 'data_' in content or '_atom_site' in content:
                return 'structure'  # CIF format
            elif 'HEADER' in content or 'ATOM' in content:
                return 'structure'  # PDB format
            elif '\t' in first_lines[0] and len(first_lines[0].split('\t')) > 3:
                return 'property'  # TSV table
            elif ',' in first_lines[0] and len(first_lines[0].split(',')) > 3:
                return 'property'  # CSV table
            
        except Exception:
            pass
        
        return None
    
    def _extract_entity_name(self, file_path: Path) -> str:
        """Extract entity name from filename."""
        # Remove extension and clean up
        name = file_path.stem
        
        # Remove common suffixes
        for suffix in ['_processed', '_clean', '_final', '_v1', '_v2']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        
        return name
    
    def scan_input_folder(self) -> List[InputFile]:
        """Scan input folder for new files to process."""
        input_files = []
        
        # Skip hidden files and README
        skip_patterns = ['.', '_', 'README', 'processed', 'rejected']
        
        for file_path in self.input_dir.iterdir():
            if file_path.is_file():
                # Skip if matches skip pattern
                if any(file_path.name.startswith(p) for p in skip_patterns):
                    continue
                
                file_type = self._detect_file_type(file_path)
                if file_type:
                    stat = file_path.stat()
                    input_files.append(InputFile(
                        path=file_path,
                        file_type=file_type,
                        entity_name=self._extract_entity_name(file_path),
                        size=stat.st_size,
                        modified=datetime.fromtimestamp(stat.st_mtime)
                    ))
                else:
                    self.logger.warning(f"Unknown file type: {file_path.name}")
        
        return sorted(input_files, key=lambda x: x.path.name)
    
    def _compute_file_hash(self, file_path: Path, chunk_size: int = 8192) -> str:
        """Compute SHA256 hash of file content."""
        hasher = hashlib.sha256()
        with open(file_path, 'rb') as f:
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
        return hasher.hexdigest()
    
    def _get_processor_for_type(self, file_type: str):
        """Get the appropriate processor for a file type."""
        # Import here to avoid circular imports
        if file_type == 'structure':
            from protos.processing.structure.structure_processor import StructureProcessor
            return StructureProcessor(paths=self.paths)
        elif file_type == 'sequence':
            from protos.processing.sequence.sequence_processor import SequenceProcessor
            return SequenceProcessor(paths=self.paths)
        elif file_type == 'property':
            from protos.processing.property.property_processor import PropertyProcessor
            return PropertyProcessor(paths=self.paths)
        elif file_type == 'embedding':
            from protos.processing.embedding.embedding_processor import EmbeddingProcessor
            return EmbeddingProcessor(paths=self.paths)
        else:
            raise ValueError(f"Unknown file type: {file_type}")
    
    def process_input_files(self, 
                          conflict_strategy: ConflictResolutionStrategy = ConflictResolutionStrategy.SKIP,
                          dry_run: bool = False) -> ProcessingReport:
        """Process all files in input folder."""
        report = ProcessingReport()
        
        for input_file in self.scan_input_folder():
            try:
                # Get appropriate processor
                processor = self._get_processor_for_type(input_file.file_type)
                
                # Compute content hash for duplicate detection
                content_hash = self._compute_file_hash(input_file.path)
                
                # Check for duplicate content
                existing_entities = self.entity_registry.find_by_content_hash(content_hash)
                if existing_entities:
                    # Duplicate content found
                    result = ProcessingResult(
                        input_file=input_file,
                        success=True,
                        action_taken="skipped",
                        entity_name=existing_entities[0].original_id,
                        details={"reason": "duplicate_content", "existing": existing_entities[0].original_id}
                    )
                    report.add_result(input_file, result)
                    
                    if not dry_run and conflict_strategy == ConflictResolutionStrategy.SKIP:
                        self._move_to_processed(input_file, "duplicate")
                    continue
                
                # Check for name conflicts
                if self.entity_registry.entity_exists(input_file.entity_name, input_file.file_type):
                    if conflict_strategy == ConflictResolutionStrategy.SKIP:
                        result = ProcessingResult(
                            input_file=input_file,
                            success=True,
                            action_taken="skipped",
                            entity_name=input_file.entity_name,
                            details={"reason": "name_exists"}
                        )
                        report.add_result(input_file, result)
                        if not dry_run:
                            self._move_to_processed(input_file, "name_conflict")
                        continue
                    elif conflict_strategy == ConflictResolutionStrategy.VERSION:
                        # Generate versioned name
                        input_file.entity_name = self._generate_version_name(input_file.entity_name)
                
                # Validate file content
                validation_errors = self._validate_file(input_file, processor)
                if validation_errors:
                    self._move_to_rejected(input_file, validation_errors)
                    report.add_rejected(input_file, validation_errors)
                    continue
                
                # Register the file
                if not dry_run:
                    # Copy to appropriate directory
                    target_path = self._get_target_path(processor, input_file)
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(str(input_file.path), str(target_path))
                    
                    # Register with entity registry
                    self.entity_registry.register_entity(
                        name=input_file.entity_name,
                        format_type=input_file.file_type,
                        file_path=str(target_path.relative_to(self.paths.data_root)),
                        metadata={
                            "source": "user_input",
                            "original_name": input_file.path.name,
                            "registered": datetime.now().isoformat(),
                            "content_hash": content_hash
                        }
                    )
                    
                    # Move to processed
                    self._move_to_processed(input_file)
                
                result = ProcessingResult(
                    input_file=input_file,
                    success=True,
                    action_taken="registered" if not dry_run else "would_register",
                    entity_name=input_file.entity_name
                )
                report.add_result(input_file, result)
                
            except Exception as e:
                self.logger.exception(f"Error processing {input_file.path}")
                report.add_error(input_file, str(e))
        
        return report
    
    def _validate_file(self, input_file: InputFile, processor) -> List[str]:
        """Validate file before registration."""
        errors = []
        
        # Check file size
        if input_file.size == 0:
            errors.append("File is empty")
        elif input_file.size > 1_000_000_000:  # 1GB limit
            errors.append("File too large (>1GB)")
        
        # Type-specific validation
        if input_file.file_type == 'structure':
            errors.extend(self._validate_structure_file(input_file.path))
        elif input_file.file_type == 'sequence':
            errors.extend(self._validate_sequence_file(input_file.path))
        elif input_file.file_type == 'property':
            errors.extend(self._validate_property_file(input_file.path))
        
        return errors
    
    def _validate_structure_file(self, file_path: Path) -> List[str]:
        """Validate structure file format."""
        errors = []
        try:
            # Use format-specific validation from formats.py if available
            fmt = format_registry.get_format_by_path(file_path)
            
            if file_path.suffix.lower() in ['.cif', '.mmcif']:
                # Use CIF validation from cif_utils if available
                try:
                    from protos.io.cif_utils import validate_cif_file
                    is_valid, error_msg = validate_cif_file(str(file_path))
                    if not is_valid:
                        errors.append(error_msg)
                except ImportError:
                    # Fallback to basic validation
                    with open(file_path, 'r') as f:
                        first_line = f.readline().strip()
                        if not first_line.startswith('data_'):
                            errors.append("Invalid CIF format: missing data_ header")
            elif file_path.suffix.lower() == '.pdb':
                # Basic PDB validation
                with open(file_path, 'r') as f:
                    content = f.read(1000)
                    if 'ATOM' not in content and 'HEADER' not in content:
                        errors.append("Invalid PDB format: no ATOM or HEADER records found")
        except Exception as e:
            errors.append(f"Cannot read file: {e}")
        return errors
    
    def _validate_sequence_file(self, file_path: Path) -> List[str]:
        """Validate sequence file format."""
        errors = []
        try:
            # Use existing validation from file_utils
            from protos.io.file_utils import validate_fasta_format
            if not validate_fasta_format(str(file_path)):
                errors.append("Invalid FASTA format")
            else:
                # Additional validation for amino acid sequences
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    for i, line in enumerate(lines):
                        line = line.strip()
                        if line and not line.startswith('>'):
                            # Check for valid amino acids
                            if not all(c in 'ACDEFGHIKLMNPQRSTVWYX-' for c in line.upper()):
                                errors.append(f"Invalid amino acid characters in sequence at line {i+1}")
                                break
        except Exception as e:
            errors.append(f"Cannot read file: {e}")
        return errors
    
    def _validate_property_file(self, file_path: Path) -> List[str]:
        """Validate property table format."""
        errors = []
        try:
            import pandas as pd
            # Try to read as CSV/TSV
            sep = '\t' if file_path.suffix == '.tsv' else ','
            df = pd.read_csv(file_path, sep=sep, index_col=0, nrows=5)
            if df.empty:
                errors.append("Table is empty")
            elif len(df.columns) == 0:
                errors.append("No data columns found")
        except Exception as e:
            errors.append(f"Cannot parse table: {e}")
        return errors
    
    def _get_target_path(self, processor, input_file: InputFile) -> Path:
        """Get the target path for the registered file."""
        # Get processor-specific directory
        if input_file.file_type == 'structure':
            base_dir = Path(self.paths.get_processor_path('structure')) / 'mmcif'
        elif input_file.file_type == 'sequence':
            base_dir = Path(self.paths.get_processor_path('sequence')) / 'fasta'
        elif input_file.file_type == 'property':
            base_dir = Path(self.paths.get_processor_path('property')) / 'tables'
        elif input_file.file_type == 'embedding':
            base_dir = Path(self.paths.get_processor_path('embedding')) / 'embeddings'
        else:
            base_dir = Path(self.paths.get_processor_path(input_file.file_type))
        
        # Keep original extension
        target_name = f"{input_file.entity_name}{input_file.path.suffix}"
        return base_dir / target_name
    
    def _generate_version_name(self, base_name: str) -> str:
        """Generate a versioned name that doesn't conflict."""
        version = 2
        while True:
            versioned_name = f"{base_name}_v{version}"
            if not self.entity_registry.entity_exists(versioned_name):
                return versioned_name
            version += 1
    
    def _move_to_processed(self, input_file: InputFile, suffix: str = None):
        """Move successfully processed file to processed folder."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if suffix:
            processed_name = f"{timestamp}_{suffix}_{input_file.path.name}"
        else:
            processed_name = f"{timestamp}_{input_file.path.name}"
        processed_path = self.processed_dir / processed_name
        shutil.move(str(input_file.path), str(processed_path))
    
    def _move_to_rejected(self, input_file: InputFile, errors: List[str]):
        """Move rejected file with error log."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        rejected_name = f"{timestamp}_{input_file.path.name}"
        rejected_path = self.rejected_dir / rejected_name
        error_log = self.rejected_dir / f"{timestamp}_{input_file.path.stem}_errors.txt"
        
        # Move file
        shutil.move(str(input_file.path), str(rejected_path))
        
        # Write error log
        error_content = f"File: {input_file.path.name}\n"
        error_content += f"Type: {input_file.file_type}\n"
        error_content += f"Time: {datetime.now().isoformat()}\n"
        error_content += "\nErrors:\n"
        error_content += "\n".join(f"- {error}" for error in errors)
        error_log.write_text(error_content)
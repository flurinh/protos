"""
Simplified Input Manager for Protos.

This module provides a streamlined workflow for registering user data:
- Flat input folder structure
- Automatic format detection
- Delegates to appropriate loaders
- Handles dataset creation for multiple files
"""

import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict
from datetime import datetime
import logging

from protos.io.paths import get_protos_paths
from protos.io.formats.format_registry import format_registry, ProcessorType


class InputManager:
    """
    Lightweight coordinator for processing input files.
    
    Key responsibilities:
    - Scan input folder and detect file types
    - Group files by processor type
    - Delegate to appropriate loaders
    - Handle dataset creation for multiple files
    """
    
    def __init__(self):
        """Initialize the input manager."""
        self.paths = get_protos_paths()
        self.input_dir = Path(self.paths.data_root) / "input"
        self.input_dir.mkdir(exist_ok=True)
        
        # Set up logging
        self.logger = logging.getLogger(__name__)
        
        # Create README
        self._create_readme()
        
    def _create_readme(self):
        """Create README in input folder."""
        readme_path = self.input_dir / "README.txt"
        if not readme_path.exists():
            readme_content = """PROTOS INPUT FOLDER
===================

Place your data files here for registration.

Usage:
1. Copy files to this folder
2. Run: protos register-input

For single files:
- File will be registered with its filename as the entity name

For multiple files of the same type:
- You'll be prompted for a dataset name
- All files will be registered and added to the dataset

Supported formats:
- Structures: .cif, .pdb
- Sequences: .fasta, .fa
- Properties: .csv, .tsv (you'll specify the type)
- Ligands: .sdf, .mol
- And more...

Files are automatically moved to the correct location after registration.
"""
            readme_path.write_text(readme_content)
    
    def scan_folder(self) -> Dict[ProcessorType, List[Path]]:
        """
        Scan input folder and group files by processor type.
        
        Returns:
            Dictionary mapping processor types to file lists
        """
        files_by_type = defaultdict(list)
        
        for file_path in self.input_dir.iterdir():
            # Skip hidden files, README, and directories
            if (file_path.is_dir() or 
                file_path.name.startswith('.') or 
                file_path.name == 'README.txt'):
                continue
                
            # Detect format
            file_format = format_registry.get_format_by_path(file_path)
            if file_format:
                files_by_type[file_format.processor].append(file_path)
            else:
                self.logger.warning(f"Unknown file format: {file_path.name}")
        
        return dict(files_by_type)
    
    def process_folder(self, 
                      dataset_names: Optional[Dict[ProcessorType, str]] = None,
                      interactive: bool = True) -> Dict[str, Any]:
        """
        Process all files in the input folder.
        
        Args:
            dataset_names: Optional dict mapping processor types to dataset names
            interactive: Whether to prompt for dataset names
            
        Returns:
            Processing summary with results
        """
        # Scan folder
        files_by_type = self.scan_folder()
        
        if not files_by_type:
            self.logger.info("No files found in input folder")
            return {"processed": 0, "errors": []}
        
        # Display what was found
        print("\nFiles found in input folder:")
        for proc_type, files in files_by_type.items():
            print(f"\n{proc_type.value.upper()} files ({len(files)}):")
            for f in files:
                print(f"  - {f.name}")
        
        # Process each type
        results = {
            "processed": 0,
            "failed": 0,
            "datasets_created": [],
            "errors": []
        }
        
        for proc_type, files in files_by_type.items():
            try:
                # Get the appropriate loader
                loader = self._get_loader(proc_type)
                if not loader:
                    self.logger.error(f"No loader available for {proc_type.value}")
                    results["errors"].append(f"No loader for {proc_type.value}")
                    continue
                
                # Handle single vs multiple files
                if len(files) == 1:
                    # Single file - just import it
                    file_path = files[0]
                    entity_name = file_path.stem  # Use filename without extension
                    
                    print(f"\nProcessing {file_path.name} as '{entity_name}'...")
                    success = self._process_single_file(loader, file_path, entity_name)
                    
                    if success:
                        results["processed"] += 1
                    else:
                        results["failed"] += 1
                        
                else:
                    # Multiple files - create dataset
                    dataset_name = None
                    
                    # Get dataset name
                    if dataset_names and proc_type in dataset_names:
                        dataset_name = dataset_names[proc_type]
                    elif interactive:
                        dataset_name = input(
                            f"\nEnter dataset name for {len(files)} {proc_type.value} files: "
                        ).strip()
                    
                    if not dataset_name:
                        dataset_name = f"{proc_type.value}_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                        print(f"Using default dataset name: {dataset_name}")
                    
                    print(f"\nProcessing {len(files)} files into dataset '{dataset_name}'...")
                    processed, failed = self._process_multiple_files(
                        loader, files, dataset_name
                    )
                    
                    results["processed"] += len(processed)
                    results["failed"] += len(failed)
                    if processed:
                        results["datasets_created"].append(dataset_name)
                        
            except Exception as e:
                self.logger.exception(f"Error processing {proc_type.value} files")
                results["errors"].append(f"{proc_type.value}: {str(e)}")
        
        # Summary
        print("\n" + "="*50)
        print("PROCESSING COMPLETE")
        print(f"Files processed: {results['processed']}")
        print(f"Files failed: {results['failed']}")
        if results['datasets_created']:
            print(f"Datasets created: {', '.join(results['datasets_created'])}")
        if results['errors']:
            print("\nErrors:")
            for error in results['errors']:
                print(f"  - {error}")
        
        return results
    
    def _get_loader(self, proc_type: ProcessorType):
        """Get the appropriate loader for a processor type."""
        try:
            if proc_type == ProcessorType.STRUCTURE:
                from protos.io.ingest.structure_loader import StructureLoader
                return StructureLoader()
            elif proc_type == ProcessorType.SEQUENCE:
                from protos.io.ingest.sequence_loader import SequenceLoader
                return SequenceLoader()
            elif proc_type == ProcessorType.MOLECULE:
                from protos.io.ingest.ligand_loader import LigandLoader
                return LigandLoader()
            # Add more loaders as they're implemented
            else:
                return None
        except ImportError:
            return None
    
    def _process_single_file(self, loader, file_path: Path, entity_name: str) -> bool:
        """
        Process a single file using the appropriate loader.
        
        Returns:
            True if successful, False otherwise
        """
        try:
            # Use loader to import and register
            registered = loader.download_and_register(
                identifier=str(file_path),
                name=entity_name,
                source='local'
            )
            
            if registered:
                # Remove from input folder (loader has copied it)
                file_path.unlink()
                self.logger.info(f"Processed {file_path.name} as '{entity_name}'")
                return True
            else:
                self.logger.error(f"Failed to process {file_path.name}")
                return False
                
        except Exception as e:
            self.logger.exception(f"Error processing {file_path}")
            return False
    
    def _process_multiple_files(self, loader, files: List[Path], 
                               dataset_name: str) -> Tuple[List[str], List[str]]:
        """
        Process multiple files and create a dataset.
        
        Returns:
            Tuple of (successful entity names, failed file names)
        """
        successful = []
        failed = []
        
        # Process each file
        for file_path in files:
            entity_name = file_path.stem
            
            try:
                registered = loader.download_and_register(
                    identifier=str(file_path),
                    name=entity_name,
                    source='local'
                )
                
                if registered:
                    successful.append(registered)
                    file_path.unlink()  # Remove from input
                else:
                    failed.append(file_path.name)
                    
            except Exception as e:
                self.logger.exception(f"Error processing {file_path}")
                failed.append(file_path.name)
        
        # Create dataset if we have successes
        if successful:
            try:
                loader.dataset_manager.create_dataset(
                    name=dataset_name,
                    entities=successful,
                    metadata={
                        'source': 'input_folder',
                        'import_date': datetime.now().isoformat(),
                        'file_count': len(successful)
                    }
                )
                self.logger.info(f"Created dataset '{dataset_name}' with {len(successful)} entities")
            except Exception as e:
                self.logger.error(f"Failed to create dataset: {e}")
        
        return successful, [f for f in files if f.name in failed]
    
    def process_ambiguous_files(self, files: List[Path], 
                               target_type: ProcessorType) -> Tuple[List[str], List[str]]:
        """
        Process files with ambiguous types (e.g., CSV could be property or GRN).
        
        Args:
            files: List of ambiguous files
            target_type: The processor type to use
            
        Returns:
            Tuple of (successful, failed)
        """
        loader = self._get_loader(target_type)
        if not loader:
            return [], [str(f) for f in files]
        
        # Determine if we need a dataset
        if len(files) == 1:
            return self._process_single_file(loader, files[0], files[0].stem), []
        else:
            dataset_name = f"{target_type.value}_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            return self._process_multiple_files(loader, files, dataset_name)


# CLI-friendly interface
def register_input(dataset_names: Optional[Dict[str, str]] = None, 
                   interactive: bool = True) -> Dict[str, Any]:
    """
    Process files in the input folder.
    
    This is the main entry point for the CLI command.
    
    Args:
        dataset_names: Optional mapping of processor types to dataset names
        interactive: Whether to prompt for dataset names
        
    Returns:
        Processing results
    """
    manager = InputManager()
    return manager.process_folder(dataset_names, interactive)
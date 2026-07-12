"""
Structure loader for fetching protein structures from various sources.

This loader handles downloading structures from:
- RCSB PDB (Protein Data Bank)
- AlphaFold Database
- Local files

The loader automatically registers downloaded structures with the entity registry,
making them immediately available for processing.
"""

import gzip
import re
import requests
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from Bio.PDB import PDBList

from protos.io.core.base_loader import BaseLoader
from protos.io.formats.cif_utils import load_structure_from_cif
from protos.processing.structure import StructureProcessor


class StructureLoader(BaseLoader):
    """
    Loader for protein structure data from multiple sources.
    
    Supports:
    - RCSB PDB structures in mmCIF format
    - AlphaFold predicted structures
    - Local file imports
    """
    
    loader_type = "structure"
    
    def __init__(self, name: str = "structure_loader", *, processor: Optional[StructureProcessor] = None):
        """Initialize structure loader."""
        super().__init__(name=name)
        self._processor: Optional[StructureProcessor] = processor
        
        # Initialize PDBList for RCSB downloads
        self.pdb_list = PDBList()
        
        # Define subdirectories for different sources
        self.source_dirs = {
            'rcsb': 'mmcif',
            'alphafold': 'mmcif/alphafold',
            'local': 'mmcif/imported'
        }
        
        # Allow multiple user-facing names to map onto the canonical sources
        self.source_aliases = {
            'rcsb': 'rcsb',
            'pdb': 'rcsb',
            'cif': 'rcsb',
            'mmcif': 'rcsb',
            'alphafold': 'alphafold',
            'alpha_fold': 'alphafold',
            'alphafold_db': 'alphafold',
            'af': 'alphafold',
            'local': 'local',
            'file': 'local',
            'filesystem': 'local',
        }
    
    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        """
        Parse structure identifier to determine source and ID.
        
        Args:
            identifier: Structure identifier (PDB ID, UniProt ID, or AlphaFold ID)
            
        Returns:
            Dict with parsed information
            
        Examples:
            - "1ubq" -> {"id": "1ubq", "source": "rcsb", "type": "experimental"}
            - "AF-P00533-F1" -> {"id": "P00533", "source": "alphafold", "version": 1}
            - "P00533" -> {"id": "P00533", "source": "alphafold", "version": "latest"}
        """
        # AlphaFold identifier pattern: AF-{uniprot_id}-F1-model_v{version}
        af_pattern = re.compile(r'^AF-([A-Z0-9]+)-F1(?:-model_v(\d+))?$')
        
        # Standard PDB ID pattern: 4 characters (alphanumeric)
        pdb_pattern = re.compile(r'^[0-9][A-Za-z0-9]{3}$')
        
        # UniProt ID pattern
        uniprot_pattern = re.compile(r'^[A-Z][0-9][A-Z0-9]{3}[0-9](?:-\d+)?$')
        
        # Check for AlphaFold format
        af_match = af_pattern.match(identifier)
        if af_match:
            return {
                'id': af_match.group(1),
                'source': 'alphafold',
                'version': int(af_match.group(2)) if af_match.group(2) else 4,  # Default v4
                'original_id': identifier
            }
        
        # Check for PDB format
        if pdb_pattern.match(identifier.lower()):
            return {
                'id': identifier.lower(),
                'source': 'rcsb',
                'type': 'experimental',
                'original_id': identifier
            }
        
        # Check for UniProt format (assume AlphaFold)
        if uniprot_pattern.match(identifier):
            return {
                'id': identifier,
                'source': 'alphafold', 
                'version': 4,  # Default to latest version
                'original_id': identifier
            }
        
        # Check if it's a file path (for local files)
        if '/' in identifier or '\\' in identifier or identifier.endswith('.cif') or identifier.endswith('.pdb'):
            # Extract filename without path and extension as ID
            from pathlib import Path
            path = Path(identifier)
            return {
                'id': path.stem,
                'source': 'local',
                'type': 'local_file',
                'original_id': identifier,
                'filename': path.name
            }
        
        # Unknown format - might be a custom identifier
        # Return it as-is for local/custom use
        return {
            'id': identifier,
            'source': 'unknown',
            'type': 'custom',
            'original_id': identifier
        }
    
    def _normalize_source(self, source: Optional[str]) -> Optional[str]:
        """Normalize user-provided source names to canonical loader keys."""
        if source is None:
            return None
        normalized = source.strip().lower()
        return self.source_aliases.get(normalized, normalized)
    
    def fetch_entity(self, identifier: str, source: Optional[str] = None, **kwargs) -> Optional[Path]:
        """
        Fetch a structure from the appropriate source.
        
        Args:
            identifier: Structure identifier
            source: Force specific source ('rcsb', 'alphafold', or 'local')
            **kwargs: Additional parameters (version, overwrite, etc.)
            
        Returns:
            Path to downloaded file, or None if failed
        """
        # Parse identifier if source not specified
        if source is None:
            try:
                id_info = self.parse_identifier(identifier)
                source = id_info['source']
            except ValueError:
                # Try RCSB by default for unknown formats
                source = 'rcsb'
        else:
            source = self._normalize_source(source)

        # Route to appropriate download method
        if source == 'rcsb':
            return self._fetch_from_rcsb(identifier, **kwargs)
        elif source == 'alphafold':
            return self._fetch_from_alphafold(identifier, **kwargs)
        elif source == 'local':
            return self._import_local_file(identifier, **kwargs)
        else:
            raise ValueError(f"Unknown source: {source}")
    
    def _fetch_from_rcsb(self, pdb_id: str, overwrite: bool = False, **kwargs) -> Optional[Path]:
        """
        Fetch structure from RCSB PDB.
        
        Args:
            pdb_id: 4-character PDB ID
            overwrite: Whether to overwrite existing files
            
        Returns:
            Path to downloaded file, or None if failed
        """
        pdb_id = pdb_id.lower()
        
        # Determine target path
        target_dir = self.get_download_path(pdb_id, self.source_dirs['rcsb'])
        target_path = target_dir / f"{pdb_id}.cif"
        
        # Check if already exists
        if target_path.exists() and not overwrite:
            self.logger.info(f"Structure already exists: {target_path}")
            return target_path
        
        try:
            # Download using BioPython PDBList
            temp_dir = self.get_download_path(pdb_id, 'temp')
            
            file_path = self.pdb_list.retrieve_pdb_file(
                pdb_id,
                pdir=str(temp_dir),
                file_format="mmCif",
                overwrite=True
            )

            if file_path and Path(file_path).exists():
                retrieved_path = Path(file_path)
                if retrieved_path.suffix == '.gz':
                    with gzip.open(retrieved_path, 'rb') as gz, open(target_path, 'wb') as out:
                        shutil.copyfileobj(gz, out)
                    retrieved_path.unlink(missing_ok=True)
                else:
                    shutil.move(str(retrieved_path), target_path)

                # Clean up temp directory structure
                self._cleanup_pdb_temp_dir(temp_dir, pdb_id)

                self.logger.info(f"Downloaded structure {pdb_id} to {target_path}")
                return target_path
            else:
                self.logger.error(f"Failed to download {pdb_id} from RCSB")
                return None
                
        except Exception as e:
            self.logger.error(f"Error downloading {pdb_id} from RCSB: {e}")
            return None
    
    def _fetch_from_alphafold(self, identifier: str, version: Optional[int] = None, 
                             overwrite: bool = False, **kwargs) -> Optional[Path]:
        """
        Fetch structure from AlphaFold database.
        
        Args:
            identifier: UniProt ID or AlphaFold ID
            version: Model version (default: 4)
            overwrite: Whether to overwrite existing files
            
        Returns:
            Path to downloaded file, or None if failed
        """
        # Parse identifier to get UniProt ID
        id_info = self.parse_identifier(identifier)
        uniprot_id = id_info['id']
        
        if version is None:
            version = id_info.get('version', 4)
        
        # Determine target path
        af_id = f"AF-{uniprot_id}-F1-model_v{version}"
        target_dir = self.get_download_path(af_id, self.source_dirs['alphafold'])
        target_path = target_dir / f"{af_id}.cif"
        
        # Check if already exists
        if target_path.exists() and not overwrite:
            self.logger.info(f"AlphaFold structure already exists: {target_path}")
            return target_path
        
        # Construct URL
        url = f"https://alphafold.ebi.ac.uk/files/{af_id}.cif"
        
        try:
            # Download file
            self.logger.info(f"Downloading AlphaFold structure from {url}")
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # Save to file
                with open(target_path, 'wb') as f:
                    f.write(response.content)
                
                self.logger.info(f"Downloaded AlphaFold structure to {target_path}")
                return target_path
            else:
                self.logger.error(
                    f"Failed to download AlphaFold structure {af_id}. "
                    f"Status: {response.status_code}"
                )
                return None
                
        except Exception as e:
            self.logger.error(f"Error downloading AlphaFold structure {af_id}: {e}")
            return None
    
    def _import_local_file(self, file_path: str, name: Optional[str] = None, **kwargs) -> Optional[Path]:
        """
        Import a local structure file.
        
        Args:
            file_path: Path to local file
            name: Name for the imported structure
            
        Returns:
            Path to imported file, or None if failed
        """
        source_path = Path(file_path)
        
        if not source_path.exists():
            self.logger.error(f"Local file not found: {file_path}")
            return None
        
        # Determine name from file if not provided
        if name is None:
            name = source_path.stem
        
        # Copy to managed location
        target_dir = self.get_download_path(name, self.source_dirs['local'])
        target_path = target_dir / f"{name}.cif"
        
        try:
            # Handle gzipped files
            if source_path.suffix == '.gz':
                with gzip.open(source_path, 'rb') as f_in:
                    with open(target_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                shutil.copy2(source_path, target_path)
            
            self.logger.info(f"Imported local structure to {target_path}")
            return target_path
            
        except Exception as e:
            self.logger.error(f"Error importing local file {file_path}: {e}")
            return None
    
    def _cleanup_pdb_temp_dir(self, temp_dir: Path, pdb_id: str):
        """Clean up PDBList temporary directory structure."""
        try:
            # PDBList creates subdirectories like 'ab' for '1abc'
            if len(pdb_id) >= 3:
                subdir = temp_dir / pdb_id[1:3]
                if subdir.exists() and not any(subdir.iterdir()):
                    subdir.rmdir()
        except Exception:
            # Ignore cleanup errors
            pass
    
    def list_sources(self) -> List[str]:
        """List available structure sources."""
        sources = list(self.source_dirs.keys())
        for alias in sorted(self.source_aliases.keys()):
            if alias not in sources:
                sources.append(alias)
        return sources
    
    def download_and_register(
        self,
        identifier: str,
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> Optional[str]:
        """Fetch, canonicalize, persist, and register a structure atomically.

        Unlike the generic loader path, structures are not usable from their raw
        CIF payload alone.  A successful return therefore means that the registry
        points at a canonical PKL which can immediately be loaded by the structure
        processor.
        """
        registered = name or identifier
        processor = self._get_processor()
        previous_info = self.entity_registry.find_entity(registered, self.loader_type)

        if previous_info and self._is_canonical_registration(previous_info.file_path):
            loaded = processor.load_entity(registered)
            if loaded is not None and not loaded.empty:
                return registered

        raw_path = self._registered_source_path(previous_info)
        fetched = False
        if raw_path is None:
            raw_path = self.fetch_entity(identifier, **kwargs)
            fetched = True
        if raw_path is None:
            self.logger.error("Failed to fetch entity: %s", identifier)
            return None

        raw_path = Path(raw_path)
        pkl_path = processor.path_pkl_dir / f"{registered}.pkl"
        pkl_existed = pkl_path.exists()
        backup_path: Optional[Path] = None
        if pkl_existed:
            with tempfile.NamedTemporaryFile(
                dir=pkl_path.parent,
                prefix=f".{registered}.",
                suffix=".pkl.bak",
                delete=False,
            ) as backup:
                backup_path = Path(backup.name)
            shutil.copy2(pkl_path, backup_path)

        try:
            df = load_structure_from_cif(str(raw_path), structure_id=registered)
            if df.empty:
                raise ValueError("CIF parsing produced an empty structure")

            combined_metadata = self._build_registration_metadata(
                identifier,
                raw_path,
                previous_info,
                metadata,
                fetched=fetched,
            )
            processor.save_entity(registered, df, metadata=combined_metadata)

            entity_info = self.entity_registry.find_entity(registered, self.loader_type)
            if not entity_info or not self._is_canonical_registration(entity_info.file_path):
                raise RuntimeError("structure registry does not point to a canonical PKL")

            registered_path = self._resolve_registered_path(entity_info.file_path)
            if not registered_path.is_file():
                raise RuntimeError(f"canonical PKL was not persisted: {registered_path}")

            # Force the success check through disk rather than the processor's frame cache.
            processor._remove_frame(registered)
            loaded = processor.load_entity(registered)
            if loaded is None or loaded.empty:
                raise RuntimeError("canonical PKL could not be loaded after registration")
        except Exception as exc:
            self.logger.warning(
                "Failed to canonicalize structure '%s': %s",
                registered,
                exc,
            )
            processor._remove_frame(registered)
            self._restore_registration(registered, previous_info)
            if backup_path is not None:
                shutil.move(str(backup_path), pkl_path)
            elif not pkl_existed:
                pkl_path.unlink(missing_ok=True)
            return None
        finally:
            if backup_path is not None:
                backup_path.unlink(missing_ok=True)

        return registered

    @staticmethod
    def _is_canonical_registration(file_path: str) -> bool:
        return Path(file_path).suffix.lower() == ".pkl"

    def _registered_source_path(self, entity_info: Any) -> Optional[Path]:
        """Return a usable CIF from an existing raw or repaired registration."""
        if not entity_info:
            return None

        candidates = [entity_info.file_path]
        source_file = (entity_info.metadata or {}).get("source_file")
        if source_file:
            candidates.append(source_file)

        for candidate in candidates:
            path = self._resolve_registered_path(candidate)
            suffixes = "".join(path.suffixes).lower()
            if path.is_file() and any(
                suffixes.endswith(ext)
                for ext in (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")
            ):
                return path
        return None

    def _build_registration_metadata(
        self,
        identifier: str,
        raw_path: Path,
        previous_info: Any,
        metadata: Optional[Dict[str, Any]],
        *,
        fetched: bool,
    ) -> Dict[str, Any]:
        combined = dict(previous_info.metadata or {}) if previous_info else {}
        combined.update(metadata or {})
        combined.update(self.parse_identifier(identifier))
        combined.update(
            {
                "source_id": identifier,
                "loader": self.__class__.__name__,
                "source_file": str(raw_path),
            }
        )
        if fetched or "download_date" not in combined:
            combined["download_date"] = datetime.now().isoformat()
        return combined

    def _restore_registration(self, name: str, previous_info: Any) -> None:
        """Restore the registry entry that existed before canonicalization."""
        if previous_info is None:
            self.entity_registry.remove_format(name, self.loader_type)
            return
        self.entity_registry.register_entity(
            name=name,
            format_type=self.loader_type,
            file_path=previous_info.file_path,
            metadata=dict(previous_info.metadata or {}),
        )

    def _resolve_registered_path(self, raw_path: str) -> Path:
        """Resolve a registered file path across platforms."""

        candidates: List[Path] = []
        path_obj = Path(raw_path)

        # As-is candidate
        candidates.append(path_obj)

        # Relative to data root
        if not path_obj.is_absolute():
            candidates.append(Path(self.paths.data_root) / path_obj)

        # Convert WSL-style paths (/mnt/<drive>/...) to Windows drive paths
        as_posix = path_obj.as_posix()
        if as_posix.startswith('/mnt/') and len(as_posix) > 6:
            drive_letter = as_posix[5]
            if drive_letter.isalpha():
                converted = Path(f"{drive_letter.upper()}:" + as_posix[6:])
                candidates.append(converted)

        for candidate in candidates:
            try:
                if candidate.exists():
                    return candidate
            except OSError:
                continue

        # Fall back to the first candidate even if it does not exist
        return candidates[-1]

    def download_and_register_alphafold(self, uniprot_id: str, 
                                      name: Optional[str] = None,
                                      version: int = 4,
                                      **kwargs) -> Optional[str]:
        """
        Convenience method for downloading AlphaFold structures.
        
        Args:
            uniprot_id: UniProt identifier
            name: Human-readable name (defaults to UniProt ID)
            version: AlphaFold models version
            
        Returns:
            Registered entity name if successful
            
        Example:
            >>> loader.download_and_register_alphafold("P00533", name="EGFR_HUMAN")
        """
        # Use UniProt ID as identifier
        return self.download_and_register(
            identifier=uniprot_id,
            name=name,
            source='alphafold',
            version=version,
            **kwargs
        )
    
    def import_local_structures(self, file_paths: List[str], 
                              names: Optional[List[str]] = None) -> Tuple[List[str], List[str]]:
        """
        Import multiple local structure files.
        
        Args:
            file_paths: List of paths to structure files
            names: Optional list of names (must match length of file_paths)
            
        Returns:
            Tuple of (successful names, failed paths)
        """
        if names and len(names) != len(file_paths):
            raise ValueError("Number of names must match number of file paths")
        
        successful = []
        failed = []
        
        for i, file_path in enumerate(file_paths):
            name = names[i] if names else None
            registered = self.download_and_register(
                identifier=file_path,
                name=name,
                source='local'
            )
            
            if registered:
                successful.append(registered)
            else:
                failed.append(file_path)
        
        return successful, failed

    def export_structure(
        self,
        structure_id: str,
        output_path: Path,
        *,
        format: str = "cif",
        overwrite: bool = True,
    ) -> Path:
        """Export a structure via the StructureProcessor."""

        processor = self._get_processor()
        return processor.export_entity(
            structure_id,
            Path(output_path),
            format=format,
            overwrite=overwrite,
        )

    def _get_processor(self) -> StructureProcessor:
        if self._processor is None:
            self._processor = StructureProcessor(name=f"{self.name}_processor")
        return self._processor
    
    def register_from_input_folder(self) -> Tuple[List[str], List[str]]:
        """
        Register structure files from the input folder using InputManager.
        
        This provides a safer workflow where users place files in a designated
        input folder for validation and registration.
        
        Returns:
            Tuple of (successful entity names, failed file names)
        """
        try:
            from protos.io.core.input_manager import InputManager
            from protos.io.core.conflict_resolver import ConflictResolutionStrategy
            
            # Create InputManager instance
            input_manager = InputManager()
            
            # Scan for structure files
            input_files = input_manager.scan_input_folder()
            structure_files = [f for f in input_files if f.file_type == 'structure']
            
            if not structure_files:
                self.logger.info("No structure files found in input folder")
                return [], []
            
            self.logger.info(f"Found {len(structure_files)} structure files in input folder")
            
            # Process the files
            report = input_manager.process_input_files(
                conflict_strategy=ConflictResolutionStrategy.SKIP
            )
            
            # Extract successful entity names
            successful = [r.entity_name for r in report.processed 
                         if r.input_file.file_type == 'structure']
            
            # Extract failed file names
            failed = [r.input_file.path.name for r in report.rejected + report.errors
                     if r.input_file.file_type == 'structure']
            
            # Display report
            report.display()
            
            return successful, failed
            
        except ImportError:
            self.logger.error("InputManager not available. Using direct import instead.")
            # Fall back to scanning a default input directory
            input_dir = self.paths.data_root / "input"
            if input_dir.exists():
                files = list(input_dir.glob("*.cif")) + list(input_dir.glob("*.pdb"))
                return self.import_local_structures([str(f) for f in files])
            return [], []

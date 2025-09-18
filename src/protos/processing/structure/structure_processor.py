from __future__ import annotations

import hashlib
import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd
import numpy as np
import requests

from protos.core.base_processor import BaseProcessor
from protos.io.cif_utils import write_cif_file
from protos.processing.structure.struct_utils import load_structure as load_structure_util, STRUCT_COLUMN_DTYPE


class StructureProcessor(BaseProcessor):
    """
    Structure processor with PKL-canonical architecture.
    
    Key principles:
    - PKL is the only canonical format (CIF is transient for ingest/export)
    - Per-entity storage with MultiIndex (structure_id, atom_id)
    - Lazy stacking with frames dict and dirty flag
    - All paths via ProtosPaths, no hardcoding
    """
    
    def __init__(
        self,
        name: str = "structure_processor",
        paths=None,
        processor_type: str = "structure",
        **kwargs
    ):
        super().__init__(
            name=name,
            paths=paths,
            processor_type=processor_type,
            **kwargs
        )
        
        # Storage: dict of DataFrames with lazy stacking
        self.frames: Dict[str, pd.DataFrame] = {}
        self._data: Optional[pd.DataFrame] = None
        self._dirty: bool = False
        
    # ---------- Path Properties ----------
    
    @property
    def path_pkl_dir(self) -> Path:
        """Get PKL directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'cache_dir'))
    
    @property
    def path_cif_dir(self) -> Path:
        """Get CIF directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'structure_dir'))
    
    @property
    def path_dataset_dir(self) -> Path:
        """Get dataset directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'dataset_dir'))
    
    @property
    def path_temp_dir(self) -> Path:
        """Get temp directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'temp_dir'))
    
    # ---------- Core Methods ----------
    
    def _ensure_canonical(self, df: pd.DataFrame, structure_id: str) -> pd.DataFrame:
        """
        Transform any DataFrame to canonical form.
        
        - Remove pdb_id column if exists (map to structure_id first)
        - Ensure structure_id and atom_id columns
        - Set MultiIndex (structure_id, atom_id)
        - Apply STRUCT_COLUMN_DTYPE with Int64 for integers
        - Ensure optional columns have proper defaults
        - Validate coordinate columns are numeric
        """
        df = df.copy()
        
        # Map pdb_id to structure_id if needed
        if 'pdb_id' in df.columns and 'structure_id' not in df.columns:
            df['structure_id'] = df['pdb_id']
        
        # Ensure structure_id column
        if 'structure_id' not in df.columns:
            df['structure_id'] = structure_id
        
        # Drop pdb_id if exists
        if 'pdb_id' in df.columns:
            df = df.drop(columns=['pdb_id'])
        
        # Ensure atom_id column
        if 'atom_id' not in df.columns:
            raise ValueError("DataFrame must contain 'atom_id' column")
        
        # Apply data types
        for col, dtype in STRUCT_COLUMN_DTYPE.items():
            if col in df.columns:
                try:
                    # Handle nullable integers
                    if dtype == 'int' or dtype == int:
                        df[col] = pd.array(df[col], dtype='Int64')
                    else:
                        df[col] = df[col].astype(dtype)
                except Exception:
                    pass
        
        # Ensure numeric coordinates
        for coord in ['x', 'y', 'z']:
            if coord in df.columns:
                df[coord] = pd.to_numeric(df[coord], errors='coerce')
        
        # Ensure optional columns have defaults
        if 'grn' in df.columns:
            df['grn'] = df['grn'].fillna('')
        
        # Set MultiIndex
        df = df.set_index(['structure_id', 'atom_id']).sort_index()
        
        return df
    
    def compute_content_hash(self, df: pd.DataFrame) -> str:
        """
        Hash core structural content for change detection.
        
        - Select core columns: chain, seq_id, residue, atom, x, y, z
        - Sort by index for consistency
        - Generate sha256 hash
        """
        # Core columns for hashing
        core_cols = []
        for col in ['auth_chain_id', 'auth_seq_id', 'res_name', 'atom_name', 'x', 'y', 'z']:
            if col in df.columns:
                core_cols.append(col)
        
        if not core_cols:
            return hashlib.sha256(b'empty').hexdigest()
        
        # Sort by index and select core columns
        df_hash = df.sort_index()[core_cols]
        
        # Convert to bytes and hash
        content = df_hash.to_csv(index=False).encode('utf-8')
        return hashlib.sha256(content).hexdigest()
    
    def load_entity(self, structure_id: str, *, auto_ingest: bool = True) -> Optional[pd.DataFrame]:
        """
        PKL-first loading with auto-ingest fallback.
        
        1. Check if already in memory
        2. Try registry lookup (format_type="structure")
        3. Try PKL in cache dir (auto-register if found)
        4. If auto_ingest and CIF exists: ingest_cif()
        5. Update frames storage
        """
        # Check memory first
        if structure_id in self.frames:
            return self.frames[structure_id]
        
        # Try registry lookup
        entity_info = self.entity_registry.find_entity(structure_id, format_type='structure')
        if entity_info:
            pkl_path = Path(self.paths.data_root) / entity_info.file_path
            if pkl_path.exists():
                try:
                    df = pd.read_pickle(pkl_path)
                    df = self._ensure_canonical(df, structure_id)
                    self._set_frame(structure_id, df)
                    return df
                except Exception as e:
                    self.logger.warning(f"Failed to load PKL from registry for {structure_id}: {e}")
        
        # Try PKL in cache dir
        pkl_path = self.path_pkl_dir / f"{structure_id}.pkl"
        if pkl_path.exists():
            try:
                df = pd.read_pickle(pkl_path)
                df = self._ensure_canonical(df, structure_id)
                # Auto-register
                self._register_entity(structure_id, pkl_path, {'auto_discovered': True})
                self._set_frame(structure_id, df)
                return df
            except Exception as e:
                self.logger.warning(f"Failed to load PKL from cache for {structure_id}: {e}")
        
        # Try auto-ingest from CIF
        if auto_ingest:
            cif_path = self.path_cif_dir / f"{structure_id}.cif"
            if cif_path.exists():
                try:
                    return self.ingest_cif(structure_id, cif_path)
                except Exception as e:
                    self.logger.error(f"Failed to auto-ingest CIF for {structure_id}: {e}")
        
        return None
    
    def save_entity(
        self, 
        structure_id: str, 
        df: pd.DataFrame, 
        format: str = 'pkl', 
        metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Save canonical PKL and register.
        
        - Validate format='pkl' only
        - Canonicalize DataFrame
        - Compute content hash
        - Write PKL to cache dir
        - Register/update with format_type="structure"
        - Include metadata: source, hash, schema version
        """
        if format != 'pkl':
            raise ValueError("Only 'pkl' format is allowed for save_entity. Use export_entity for other formats.")
        
        # Canonicalize
        df = self._ensure_canonical(df, structure_id)
        
        # Compute hash
        content_hash = self.compute_content_hash(df)
        
        # Write PKL
        pkl_path = self.path_pkl_dir / f"{structure_id}.pkl"
        df.to_pickle(pkl_path)
        
        # Prepare metadata
        entity_metadata = {
            'content_hash': content_hash,
            'schema_version': '1.0',
            'saved_at': datetime.utcnow().isoformat(),
            'atom_count': len(df),
            **(metadata or {})
        }
        
        # Register
        self._register_entity(structure_id, pkl_path, entity_metadata)
        
        # Update storage
        self._set_frame(structure_id, df)
    
    def ingest_cif(self, structure_id: str, cif_path: Optional[Path] = None) -> pd.DataFrame:
        """
        One-time CIF to PKL conversion.
        
        - Find CIF (provided path or use load_structure_util with folder)
        - Parse using existing CIF parser
        - Canonicalize DataFrame
        - Save as PKL via save_entity()
        - Update storage
        """
        if cif_path is not None:
            # Use specific path provided
            if not cif_path.exists():
                raise FileNotFoundError(f"CIF file not found: {cif_path}")
            # Extract folder and filename parts for load_structure_util
            folder = str(cif_path.parent)
            # Remove .cif extension for pdb_id
            pdb_id = cif_path.stem
            df = load_structure_util(pdb_id, folder=folder + '/')
        else:
            # Use standard location via ProtosPaths
            folder = str(self.path_cif_dir) + '/'
            df = load_structure_util(structure_id, folder=folder)
        
        # Save as PKL (will canonicalize and register)
        self.save_entity(
            structure_id, 
            df, 
            format='pkl',
            metadata={
                'source': 'cif_ingest',
                'source_folder': folder
            }
        )
        
        return self.frames[structure_id]
    
    def export_entity(
        self, 
        structure_id: str, 
        format: str = 'cif', 
        out_path: Optional[Path] = None
    ) -> Path:
        """
        Export PKL to other formats (no registration).
        
        - Load from PKL (must exist)
        - Reset index for CIF compatibility
        - Write CIF to structure dir
        - Return path (no registry update)
        """
        if format != 'cif':
            raise ValueError("Only 'cif' format is supported for export")
        
        # Load entity
        df = self.load_entity(structure_id, auto_ingest=False)
        if df is None:
            raise ValueError(f"Structure {structure_id} not found in PKL storage")
        
        # Determine output path
        if out_path is None:
            out_path = self.path_cif_dir / f"{structure_id}.cif"
        
        # Reset index for CIF format
        df_cif = df.reset_index()
        
        # Write CIF
        write_cif_file(str(out_path), df_cif)
        
        return out_path
    
    # ---------- Frame Management ----------
    
    def _set_frame(self, structure_id: str, df: pd.DataFrame):
        """Update frame and mark dirty"""
        self.frames[structure_id] = df
        self._dirty = True
    
    def _remove_frame(self, structure_id: str):
        """Remove frame and mark dirty"""
        if structure_id in self.frames:
            del self.frames[structure_id]
            self._dirty = True
    
    @property
    def data(self) -> Optional[pd.DataFrame]:
        """Lazy stacked view with rebuild on access"""
        if self._dirty and self.frames:
            self._data = pd.concat(self.frames.values()).sort_index()
            self._dirty = False
        return self._data
    
    # ---------- Dataset Operations ----------
    
    def load_dataset(self, dataset_name: str, return_format: str = 'stacked') -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
        """
        Load all dataset members from individual PKLs.
        
        - Get member list from dataset
        - Load each via load_entity(auto_ingest=False)
        - Update frames in bulk
        - Return stacked or dict format
        """
        dataset = self.get_dataset(dataset_name)
        if not dataset:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        structure_ids = dataset.get('entities', [])
        
        # Load all structures
        for structure_id in structure_ids:
            self.load_entity(structure_id, auto_ingest=False)
        
        if return_format == 'dict':
            return {sid: self.frames[sid] for sid in structure_ids if sid in self.frames}
        else:  # 'stacked'
            return self.data
    
    def save_dataset(
        self, 
        dataset_name: str, 
        structure_ids: List[str], 
        metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Ensure PKLs exist for all members.
        
        - Check each member has registered PKL
        - Save missing ones from frames if available
        - Create logical dataset entry
        """
        # Ensure all structures are saved as PKL
        for structure_id in structure_ids:
            if structure_id in self.frames:
                # Already in memory, ensure PKL exists
                entity_info = self.entity_registry.find_entity(structure_id, format_type='structure')
                if not entity_info:
                    self.save_entity(structure_id, self.frames[structure_id])
            else:
                # Not in memory, check if registered
                entity_info = self.entity_registry.find_entity(structure_id, format_type='structure')
                if not entity_info:
                    self.logger.warning(f"Structure {structure_id} not found for dataset {dataset_name}")
        
        # Create logical dataset
        self.create_dataset(dataset_name, structure_ids, metadata or {})
    
    # ---------- Utility Methods ----------
    
    def _register_entity(self, structure_id: str, file_path: Path, metadata: Dict[str, Any]) -> None:
        """Register entity with generic 'structure' format type"""
        rel_path = file_path.relative_to(self.paths.data_root)
        self.entity_registry.register_entity(
            name=structure_id,
            format_type='structure',
            file_path=str(rel_path),
            metadata=metadata
        )
    
    @property
    def structure_ids(self) -> List[str]:
        """Get list of loaded structure IDs"""
        return list(self.frames.keys())
    
    # ---------- Download Method ----------
    
    def download_structure(
        self,
        structure_id: str,
        *,
        source: str = "rcsb",
        save_to_cache: bool = True,
        overwrite: bool = False,
        metadata: Optional[Dict[str, Any]] = None
    ) -> Optional[pd.DataFrame]:
        """
        Download structure from remote source.
        
        Following the PKL-canonical approach:
        - Download CIF to structure directory
        - Ingest to PKL if save_to_cache=True
        - Otherwise just parse and return canonicalized
        
        Args:
            structure_id: Structure identifier (PDB ID or UniProt ID)
            source: 'rcsb' or 'alphafold'
            save_to_cache: If True, save as PKL (default)
            overwrite: If True, overwrite existing files
            metadata: Additional metadata to store
            
        Returns:
            DataFrame with structure data or None if download failed
        """
        # Check if already exists and not overwriting
        if not overwrite and save_to_cache:
            existing = self.load_entity(structure_id, auto_ingest=False)
            if existing is not None:
                self.logger.info(f"Structure {structure_id} already exists in cache")
                return existing
        
        # Download to CIF directory
        cif_path = self.path_cif_dir / f"{structure_id}.cif"
        
        # Check if CIF exists and not overwriting
        if cif_path.exists() and not overwrite:
            if save_to_cache:
                # Ingest existing CIF
                return self.ingest_cif(structure_id, cif_path)
            else:
                # Just parse and return
                folder = str(self.path_cif_dir) + '/'
                df = load_structure_util(structure_id, folder=folder)
                return self._ensure_canonical(df, structure_id)
        
        # Download URL mapping
        if source == "rcsb":
            # Use existing download_structures utility
            from protos.loaders.download_structures import download_structures_with_processor
            successful, failed = download_structures_with_processor(
                pdb_ids=[structure_id],
                processor=self,
                overwrite=overwrite
            )
            
            if structure_id.lower() in failed:
                self.logger.error(f"Failed to download {structure_id} from RCSB")
                return None
                
        elif source == "alphafold":
            # Use existing alphafold_utils
            from protos.loaders.alphafold_utils import download_alphafold_with_processor
            download_alphafold_with_processor(
                uid=structure_id,
                processor=self,
                max_models=1  # Just get v4 model
            )
            
            # Check if download succeeded
            expected_file = self.path_cif_dir / f"AF-{structure_id}-F1-model_v4.cif"
            if not expected_file.exists():
                self.logger.error(f"Failed to download {structure_id} from AlphaFold")
                return None
            # Use the AlphaFold filename
            cif_path = expected_file
            
        else:
            raise ValueError(f"Unknown source: {source}")
        
        # Process the downloaded file
        if save_to_cache:
            # Parse and save with download metadata
            folder = str(self.path_cif_dir) + '/'
            df = load_structure_util(structure_id, folder=folder)
            self.save_entity(
                structure_id,
                df,
                metadata={
                    'source': source,
                    'source_file': str(cif_path),
                    'downloaded_at': datetime.utcnow().isoformat(),
                    **(metadata or {})
                }
            )
            return self.frames[structure_id]
        else:
            # Just parse and return canonicalized
            folder = str(self.path_cif_dir) + '/'
            df = load_structure_util(structure_id, folder=folder)
            return self._ensure_canonical(df, structure_id)
    
    # ---------- Deprecation Warnings ----------
    
    @property
    def pdb_ids(self) -> List[str]:
        """Deprecated: use structure_ids"""
        warnings.warn("pdb_ids is deprecated, use structure_ids", DeprecationWarning, stacklevel=2)
        return self.structure_ids
    
    def load_structure(self, identifier: str, **kwargs) -> Optional[pd.DataFrame]:
        """Deprecated: use load_entity"""
        warnings.warn("load_structure is deprecated, use load_entity", DeprecationWarning, stacklevel=2)
        return self.load_entity(identifier)
    
    def save_structure(self, name: str, structure_df: pd.DataFrame, **kwargs) -> None:
        """Deprecated: use save_entity or export_entity"""
        warnings.warn("save_structure is deprecated, use save_entity or export_entity", DeprecationWarning, stacklevel=2)
        # Check if user wants non-PKL format
        if 'format' in kwargs and kwargs['format'] != 'pkl':
            self.export_entity(name, format=kwargs['format'])
        else:
            self.save_entity(name, structure_df, metadata=kwargs.get('metadata'))
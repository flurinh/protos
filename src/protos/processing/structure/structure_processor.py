"""
Structure Processor with BaseProcessor integration.

This module provides a StructureProcessor class that extends BaseProcessor
to provide standardized data management capabilities for structural data.

UPDATED: Follows DATA_MANAGEMENT_UNIFIED.md principles
- NO custom path handling - ProtosPaths manages ALL paths
- Zero configuration required
- Human-readable names for all operations
- Implements abstract methods from BaseProcessor
"""

import os
import time
import json
import warnings
import tempfile
import shutil
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union, Set, Any, Tuple

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from protos.io import cif_utils
from protos.io.data_access import generate_entity_id
from protos.core.base_processor import BaseProcessor
from protos.processing.structure.struct_utils import (
    load_structure as load_structure_util,
    STRUCT_COLUMNS,
    STRUCT_COLUMN_DTYPE,
    SORTED_STRUCT_COLUMNS,
    ALPHA_CARBON
)


class StructureProcessor(BaseProcessor):
    """
    Processor for structural data in mmCIF format.
    
    This class extends BaseProcessor to provide standardized data management
    capabilities for structural data, such as loading and saving structures,
    filtering by various criteria, and extracting specific features.
    
    Key changes from legacy version:
    - NO path parameters in constructor
    - ALL paths managed by ProtosPaths
    - Uses self.paths instead of self.path_resolver
    - Implements abstract methods from BaseProcessor
    """
    
    def __init__(
            self,
            name="cif_processor",
            paths=None,  # ProtosPaths instance
            pdb_ids_file=None,
            limit=None,
            preload=False,
            remove_hetatm=False,
            allow_exception=False
    ):
        """
        Initialize the CIF processor.
        
        Args:
            name: Processor instance name
            paths: ProtosPaths instance (created if not provided)
            pdb_ids_file: File containing PDB IDs to load (relative to data path)
            limit: Maximum number of PDB files to process
            preload: Whether to load pdb_ids on initialization
            remove_hetatm: Whether to remove HETATM records
            allow_exception: Whether to allow exceptions during processing
        """
        # Initialize BaseProcessor with ProtosPaths
        super().__init__(name=name, paths=paths)
        
        # Track PDB ID to entity ID mappings
        self._pdb_entity_map = {}
        
        # Set up parameters
        self.pdb_ids_file = pdb_ids_file
            
        # Configure processing options
        self.limit = limit
        self.remove_hetatm = remove_hetatm
        self.allow_exception = allow_exception
        
        # Set up data containers
        self.pdb_ids = []
        self.structure_filenames = []
        self.ref_seq = None
        self.dfl = []
        self.dfl_list = []
        self.chain_dict = {}
        self.chain_dict_ca_atom_ids = {}
        self.data = None
        self.structure_info = {}
        self.chain_ids = []
        
        # Initialize PDB IDs if requested
        if preload:
            self.initialize_pdb_ids()
    
    # Path properties using ProtosPaths
    
    @property
    def path_structure_dir(self):
        """Get path to structure files directory (mmcif)."""
        return self.get_subdirectory_path('structure_dir')
    
    @property  
    def path_dataset_dir(self):
        """Get path to dataset files directory."""
        return self.get_subdirectory_path('dataset_dir')
    
    @property
    def path_alignment_dir(self):
        """Get path to alignment files directory.""" 
        return self.get_subdirectory_path('alignments_dir')
    
    @property
    def path_cache_dir(self):
        """Get path to cache directory."""
        return self.get_subdirectory_path('cache_dir')
    
    @property
    def temp_dir(self):
        """Get path to temporary files directory."""
        return self.get_subdirectory_path('temp_dir')
    
    @property
    def path_pdb_ids_file(self):
        """Get path to PDB IDs file."""
        if self.pdb_ids_file is None:
            return None
        return self.data_path / self.pdb_ids_file
    
    # Implement abstract methods from BaseProcessor
    
    def load_entity(self, name: str) -> Optional[pd.DataFrame]:
        """
        Load structure entity by PDB ID.
        
        Args:
            name: PDB ID (e.g., '1ubq', '3sn6')
            
        Returns:
            DataFrame with structure data or None if not found
        """
        # First check cache
        cache_file = self.path_cache_dir / f"{name}.pkl"
        if cache_file.exists():
            try:
                df = pd.read_pickle(cache_file)
                # Auto-register if not already registered
                if not self.entity_exists(name):
                    self.entity_registry.register_entity(
                        name=name,
                        format_type=self.processor_type,
                        file_path=str(cache_file.relative_to(self.paths.data_root)),
                        metadata={"cached": True, "auto_discovered": True}
                    )
                return df
            except Exception as e:
                self.logger.warning(f"Failed to load cache for {name}: {e}")
        
        # Try to load from mmcif
        cif_file = self.path_structure_dir / f"{name}.cif"
        if cif_file.exists():
            try:
                # Load from CIF file directly (avoid recursion)
                structure_dir = str(self.path_structure_dir)
                df = load_structure_util(name, folder=structure_dir)
                
                # Auto-register if not already registered
                if not self.entity_exists(name):
                    self.entity_registry.register_entity(
                        name=name,
                        format_type=self.processor_type,
                        file_path=str(cif_file.relative_to(self.paths.data_root)),
                        metadata={"auto_discovered": True}
                    )
                return df
            except Exception as e:
                self.logger.error(f"Failed to load structure {name}: {e}")
                return None
        
        # Check if registered
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            if file_path.exists():
                try:
                    if file_path.suffix == '.pkl':
                        return pd.read_pickle(file_path)
                    else:
                        return self.load_structure(name)
                except Exception as e:
                    self.logger.error(f"Failed to load registered entity {name}: {e}")
        
        return None
    
    def save_entity(self, name: str, data: pd.DataFrame, metadata: Optional[dict] = None):
        """
        Save structure entity.
        
        Args:
            name: PDB ID
            data: Structure DataFrame
            metadata: Optional metadata dict
        """
        # Save to cache
        cache_file = self.path_cache_dir / f"{name}.pkl"
        # ProtosPaths handles directory creation
        
        data.to_pickle(cache_file)
        
        # Register entity
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(cache_file.relative_to(self.paths.data_root)),
            metadata=metadata or {}
        )
    
    # Initialize and setup methods
    def initialize_pdb_ids(self):
        """Set up PDB IDs from file or directory listing."""
        if self.path_pdb_ids_file is not None and self.path_pdb_ids_file.exists():
            self.logger.info(f"Loading PDB ids from {self.path_pdb_ids_file}.")
            with open(self.path_pdb_ids_file) as f:
                self.pdb_ids = sorted([
                    pdb_id.strip() 
                    for pdb_id in f.readlines() 
                    if pdb_id.strip()
                ])
                
        elif self.path_structure_dir.exists():
            # Load from directory
            self.logger.info(f"Loading PDB ids from {self.path_structure_dir}.")
            self.pdb_ids = sorted([
                path.stem
                for path in self.path_structure_dir.glob("*.cif")
            ])
        else:
            self.logger.warning("No PDB IDs file or structure directory found.")
            self.pdb_ids = []
        
        # Apply limit if specified
        if self.limit is not None:
            self.pdb_ids = self.pdb_ids[:self.limit]
            
        self.logger.info(f"Loaded {len(self.pdb_ids)} PDB IDs.")
    
    def get_pdb_ids_from_file(self, path_pdb_ids_file=None):
        """
        Get PDB IDs from a file.
        
        Args:
            path_pdb_ids_file: Path to file containing PDB IDs
            
        Returns:
            List of PDB IDs
        """
        if path_pdb_ids_file is not None:
            self.pdb_ids_file = path_pdb_ids_file
            
        if self.path_pdb_ids_file is None:
            raise ValueError("No PDB IDs file specified")
        
        with open(self.path_pdb_ids_file, 'r') as f:
            pdb_ids = f.readlines()[0].split(',')
            pdb_ids = [pdb_id.strip() for pdb_id in pdb_ids]
        
        return pdb_ids
    
    def get_available_pdb_files(self):
        """
        Get PDB IDs based on available CIF files.
        
        Returns:
            List of PDB IDs from CIF files
        """
        if not self.path_structure_dir.exists():
            return []
            
        # Use pathlib for better path handling
        pdb_ids = [file.stem for file in self.path_structure_dir.glob("*.cif")]
        return sorted(pdb_ids)
    
    def load_structure(self, identifier, apply_dtypes=True, debug=False, use_cache=True, save_processed=True):
        """
        Load a structure from CIF or cached PKL file.
        
        Loading priority:
        1. If use_cache=True, check for processed PKL file
        2. If not found, check for raw CIF file and parse it
        3. If save_processed=True, save parsed structure as PKL for future use
        
        Args:
            identifier: Structure identifier (human-readable name, e.g., "1ubq", "my_protein")
            apply_dtypes: Whether to apply proper data types to the loaded structure
            debug: Whether to print debug information during data type formatting
            use_cache: Whether to check for cached PKL files first (default=True)
            save_processed: Whether to save parsed CIF as PKL for future use (default=True)
            
        Returns:
            DataFrame with structure data or None if loading fails
        """
        # Use human-readable name directly
        name = str(identifier)
        
        # Step 1: Check for cached PKL file if use_cache is True
        if use_cache:
            try:
                # Try to load from entity system first (which handles PKL files)
                structure = self.load_entity(name)
                if structure is not None and isinstance(structure, pd.DataFrame):
                    if debug:
                        self.logger.info(f"Loaded structure '{name}' from cache")
                    
                    # Apply data types if requested
                    if apply_dtypes:
                        # Store in self.data temporarily for format_data_types
                        self.data = structure
                        self.format_data_types(debug=debug)
                        structure = self.data
                    
                    return structure
            except Exception as e:
                if debug:
                    self.logger.debug(f"No cached version found for '{name}': {e}")
        
        # Step 2: Try to load from CIF file
        try:
            # Get the structure subdirectory path from ProtosPaths
            # ProtosPaths manages all path resolution
            structure_dir = str(self.path_structure_dir)
            structure = load_structure_util(name, folder=structure_dir)
            
            # Ensure the structure is properly loaded
            if structure is None or len(structure) == 0:
                self.logger.warning(f"Failed to load structure for {name}")
                return None
            
            # Add the name to our list if not already present
            if name not in self.pdb_ids:
                self.pdb_ids.append(name)
            
            # Filter HETATM records if requested
            if self.remove_hetatm:
                structure = structure[structure["group"] == "ATOM"]
            
            # Apply data types if requested
            if apply_dtypes:
                # Store in self.data temporarily for format_data_types
                original_data = self.data
                self.data = structure
                self.format_data_types(debug=debug)
                structure = self.data
                self.data = original_data  # Restore original
            
            # Step 3: Save as PKL for future use if requested
            if save_processed:
                try:
                    self.save_entity(name, structure, metadata={
                        "source": "cif",
                        "atom_count": len(structure),
                        "chains": structure["auth_chain_id"].unique().tolist() if "auth_chain_id" in structure.columns else []
                    })
                    if debug:
                        self.logger.info(f"Cached processed structure '{name}' as PKL")
                except Exception as e:
                    self.logger.warning(f"Failed to cache structure '{name}': {e}")
            
            return structure
        except Exception as e:
            if self.allow_exception:
                self.logger.error(f"Error loading structure {name}: {str(e)}")
                return None
            else:
                raise
    
    def save_structure(self, name: str, structure_df: pd.DataFrame, 
                      format: str = 'pkl', metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Save a structure with a human-readable name.
        
        Args:
            name: Human-readable name for the structure (e.g., "1ubq", "my_protein")
            structure_df: Structure DataFrame to save
            format: Save format - 'pkl' (default), 'cif', or 'both'
            metadata: Optional metadata for the structure
            
        Returns:
            The human-readable name used
        """
        # Update metadata with structure info
        if metadata is None:
            metadata = {}
        metadata.update({
            "atom_count": len(structure_df),
            "chains": structure_df["auth_chain_id"].unique().tolist() if "auth_chain_id" in structure_df.columns else [],
            "format_saved": format
        })
        
        # Save based on format
        if format in ['pkl', 'both']:
            # Use our overridden save_entity for PKL format
            self.save_entity(name, structure_df, metadata=metadata)
            
        if format in ['cif', 'both']:
            # For CIF format, save to mmcif directory
            cif_path = self.path_structure_dir / f"{name}.cif"
            # ProtosPaths handles directory creation
            
            # Convert DataFrame to CIF format
            try:
                self.save_cif(structure_df, str(cif_path))
                
                # Register the CIF file
                self.entity_registry.register_entity(
                    name=name,
                    format_type=self.processor_type,
                    file_path=str(cif_path.relative_to(self.paths.data_root)),
                    metadata={**metadata, "format": "cif"}
                )
            except Exception as e:
                self.logger.warning(f"Failed to save CIF format: {e}")
        
        return name
    
    def load_structures(self, pdb_ids=None, apply_dtypes=True, debug=False):
        """
        Load multiple structures.
        
        Args:
            pdb_ids: List of PDB IDs to load (defaults to self.pdb_ids)
            apply_dtypes: Whether to apply proper data types to loaded structures
            debug: Whether to print debug information during data type formatting
            
        Returns:
            Concatenated DataFrame of loaded structures
        """
        if pdb_ids is None:
            pdb_ids = self.pdb_ids
        
        # Clear existing data
        self.dfl = []
        
        # Load each structure with type conversion disabled initially
        # We'll do a single type conversion at the end for better performance
        for pdb_id in tqdm(pdb_ids, desc="Loading structures"):
            # Don't apply types for individual structures yet
            structure = self.load_structure(pdb_id, apply_dtypes=False)
            if structure is not None:
                self.dfl.append(structure)
        
        # Update the data structure
        self.concat_data()
        
        # Apply data types once to the entire dataset if requested
        if apply_dtypes and self.data is not None and not self.data.empty:
            self.format_data_types(debug=debug)
        
        return self.data
    
    # Data management methods
    def concat_data(self):
        """Concatenate all DataFrames in dfl into the data attribute."""
        if not self.dfl:
            return
        
        # Filter out None values and empty DataFrames
        valid_dfl = [df for df in self.dfl if df is not None and not df.empty]
        if not valid_dfl:
            return
        
        self.data = pd.concat(valid_dfl, ignore_index=True)
        self.reset_index()
    
    def reset_index(self):
        """Reset and restructure the data index."""
        if self.data is None or len(self.data) == 0:
            return
        
        self.data.reset_index(drop=True, inplace=True)
        pdb_id_mapping = {pdb_id: i for i, pdb_id in enumerate(self.data['pdb_id'].unique(), 1)}
        self.data['idx'] = self.data['pdb_id'].map(pdb_id_mapping)
        self.data.set_index(['idx', self.data.groupby('pdb_id').cumcount()], inplace=True)
        self.data.index.names = [None, None]
    
    def update_pdb_ids(self):
        """Update pdb_ids list from current data."""
        if self.data is not None:
            self.pdb_ids = self.data['pdb_id'].unique().tolist()
    
    def update_chain_data(self):
        """Update chain dictionaries after data changes."""
        self.update_pdb_chain_ids()
    
    def update_pdb_chain_ids(self):
        """Update chain_ids list from chain_dict."""
        # List of all chains in all structures (<pdb_id>_<chain_id>)
        self.chain_ids = list(self.chain_dict.keys())
    
    def format_data_types(self, debug=False):
        """
        Apply correct data types to all columns in the structure data.
        
        Args:
            debug (bool): Whether to print debug information
            
        Returns:
            None (modifies self.data in place)
        """
        if self.data is None or self.data.empty:
            self.logger.warning("No data to format")
            return
            
        if debug:
            self.logger.info(f"Formatting data types for {len(self.pdb_ids)} structures")
            
        # Apply standard column types
        for col in self.data.columns:
            if col in STRUCT_COLUMN_DTYPE:
                try:
                    self.data[col] = self.data[col].astype(STRUCT_COLUMN_DTYPE[col])
                except Exception as e:
                    self.logger.warning(f"Failed to convert column {col}: {e}")
        
        # Ensure numeric coordinates
        for coord in ['x', 'y', 'z']:
            if coord in self.data.columns:
                self.data[coord] = pd.to_numeric(self.data[coord], errors='coerce')
        
        # For integer columns, handle NaN values
        int_columns = ['auth_seq_id', 'gen_seq_id', 'atom_id', 'res_id']
        for col in int_columns:
            if col in self.data.columns:
                try:
                    # Convert to nullable integer type
                    self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
                    mask = ~pd.isna(self.data[col])
                    if mask.any():
                        int_values = pd.Series(index=self.data.index, dtype='Int64')
                        int_values.loc[mask] = self.data.loc[mask, col].astype('int64')
                        self.data[col] = int_values
                except Exception as e:
                    if debug:
                        self.logger.warning(f"Failed to convert {col} to int: {e}")
        
        if debug:
            self.logger.info("Data type formatting complete")
    
    def reset_data(self, preserve_ids: bool = False):
        """
        Reset the processor data.

        Args:
            preserve_ids: If True, keep the list of PDB IDs but reset all other data structures
        """
        # Store current PDB IDs if needed
        current_ids = self.pdb_ids.copy() if preserve_ids else []

        # Reset all data structures
        self.data = None  # Explicitly set to None
        self.dfl = []
        self.chain_dict = {}
        self.pdb_ids = current_ids  # Restore IDs if preserving
    
    # Structure analysis methods
    def get_chains(self, pdb_id):
        """
        Get chain IDs for a structure.
        
        Args:
            pdb_id: PDB ID of the structure
            
        Returns:
            List of chain IDs
        """
        if self.data is None:
            return []
        
        return self.data[self.data['pdb_id'] == pdb_id]['auth_chain_id'].unique().tolist()
    
    def get_ca_coordinates(self, pdb_id, chain_id):
        """
        Get CA (alpha carbon) coordinates for a specific chain.
        
        Args:
            pdb_id: PDB ID of the structure
            chain_id: Chain ID
            
        Returns:
            Array of coordinates
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
        
        df = self.data[(self.data['pdb_id'] == pdb_id) &
                      (self.data['auth_chain_id'] == chain_id) &
                      (self.data['res_atom_name'] == ALPHA_CARBON)]
        
        if df.empty:
            raise ValueError(f"No CA atoms found for {pdb_id} chain {chain_id}")
        
        return df[['x', 'y', 'z']].to_numpy().astype(np.float64)
    
    def get_sequence(self, pdb_id, chain_id='A'):
        """
        Get amino acid sequence for a chain.
        
        Args:
            pdb_id: PDB ID of the structure
            chain_id: Chain ID
            
        Returns:
            Amino acid sequence as string
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
        
        seq = self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id) &
                        (self.data['res_atom_name'] == ALPHA_CARBON)]['res_name1l'].tolist()
        
        return ''.join(seq)
    
    def get_seq(self, pdb_id, chain_id='A'):
        """Alias for get_sequence for backward compatibility."""
        return self.get_sequence(pdb_id, chain_id)
    
    def extract_sequence(self, chain_id='A', pdb_id=None):
        """
        Extract sequence from structure for a specific chain.
        
        This is a convenience method that wraps get_seq for easier use.
        
        Args:
            chain_id: Chain identifier (default 'A')
            pdb_id: PDB ID (if None, uses first loaded structure)
            
        Returns:
            String sequence or None if not found
        """
        if self.data is None or self.data.empty:
            self.logger.warning("No structure data loaded")
            return None
            
        # Use first PDB ID if not specified
        if pdb_id is None:
            if self.pdb_ids:
                pdb_id = self.pdb_ids[0]
            else:
                self.logger.warning("No PDB IDs available")
                return None
                
        try:
            sequence = self.get_seq(pdb_id, chain_id)
            if sequence:
                self.logger.info(f"Extracted sequence from {pdb_id} chain {chain_id} (length: {len(sequence)})")
            else:
                self.logger.warning(f"No sequence found for {pdb_id} chain {chain_id}")
            return sequence
        except Exception as e:
            self.logger.error(f"Error extracting sequence: {e}")
            return None
    
    def get_all_sequences(self) -> Dict[str, str]:
        """
        Get sequences for all loaded chains.
        
        Returns:
            Dictionary mapping chain_id to sequence
        """
        sequences = {}
        
        if self.data is None:
            return sequences
        
        # Get unique combinations of pdb_id and chain_id
        chains = self.data[['pdb_id', 'auth_chain_id']].drop_duplicates()
        
        for _, row in chains.iterrows():
            pdb_id = row['pdb_id']
            chain_id = row['auth_chain_id']
            
            try:
                seq = self.get_sequence(pdb_id, chain_id)
                if seq:
                    sequences[f"{pdb_id}_{chain_id}"] = seq
            except Exception as e:
                self.logger.warning(f"Failed to get sequence for {pdb_id}_{chain_id}: {e}")
        
        return sequences
    
    def get_seq_dict(self, load_file=False, version='v1', chain_type=None):
        """
        Get a dictionary mapping chain IDs to sequences.
        
        Args:
            load_file: Whether to load from a file
            version: Version of the data
            chain_type: Type of chain to filter for
            
        Returns:
            Dictionary mapping chain IDs to sequences
        """
        if load_file:
            return self.load_chain_dict_from_fasta(version)
            
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        chain_dict = {}
        
        # Filter data based on chain_type if specified
        if chain_type is not None and 'chain_type' in self.data.columns:
            backbone = self.data[(self.data['chain_type'] == chain_type) &
                                (self.data['res_atom_name'] == ALPHA_CARBON)]
        else:
            backbone = self.data[(self.data['res_atom_name'] == ALPHA_CARBON) &
                                (self.data['group'] == 'ATOM')]
                                
        for pdb_id in tqdm(self.pdb_ids, desc="Creating sequence dictionary"):
            chains_in_pdb = backbone[backbone['pdb_id'] == pdb_id]['auth_chain_id'].unique()
            for chain_id in chains_in_pdb:
                # Sort by sequence ID to ensure we get residues in order
                chain_data = backbone[(backbone['pdb_id'] == pdb_id) &
                                     (backbone['auth_chain_id'] == chain_id)].sort_values(by='gen_seq_id')
                                     
                seq = chain_data['res_name1l'].tolist()
                ca_atom_ids = chain_data['atom_id'].tolist()
                
                # Store atom IDs for reference
                self.chain_dict_ca_atom_ids[f"{pdb_id}_{chain_id}"] = ca_atom_ids
                
                # Store the sequence
                chain_dict[f"{pdb_id}_{chain_id}"] = ''.join(seq).replace(' ', 'X')
        self.update_pdb_chain_ids()
        self.chain_dict = chain_dict
        return chain_dict
    
    def load_chain_dict_from_fasta(self, version='v1'):
        """
        Load chain sequences from a FASTA file.
        
        Args:
            version: Version string for the file
            
        Returns:
            Dictionary mapping chain IDs to sequences
        """
        filename = self.path_dataset_dir / f'chain_dict_{version}.fas'
        if not filename.exists():
            self.logger.warning(f"Chain dictionary file {filename} not found")
            return {}
            
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Create a dictionary from the lines
        chain_dict = {lines[i][1:].strip(): lines[i + 1].strip().upper()
                     for i in range(0, len(lines), 2)}
                     
        self.chain_dict = chain_dict
        return chain_dict
    
    def save_chain_dict_to_fasta(self, version='v1'):
        """
        Save chain sequences to a FASTA file.
        
        Args:
            version: Version string for the file
        """
        if not self.chain_dict:
            self.logger.warning("No chain dictionary to save")
            return
            
        filename = self.path_dataset_dir / f'chain_dict_{version}.fas'
        # ProtosPaths handles directory creation
        
        with open(filename, 'w') as f:
            for key, value in self.chain_dict.items():
                f.write(f'>{key}\n{value}\n')
                
        self.logger.info(f"Saved chain dictionary to {filename}")
    
    # Filtering methods
    def filter_by_pdb_ids(self, pdb_ids):
        """
        Filter data by PDB IDs.
        
        Args:
            pdb_ids: List of PDB IDs to keep
        """
        if self.data is None:
            return
        
        # Convert input to a set for faster operations
        pdb_ids_set = set(pdb_ids)
        
        # Check if all requested pdb_ids are present
        present_pdb_ids = set(self.data['pdb_id'].unique())
        missing_pdb_ids = pdb_ids_set - present_pdb_ids
        
        if missing_pdb_ids:
            self.logger.warning(f"The following PDB IDs are missing: {missing_pdb_ids}")
            
        # Filter pdb_ids
        self.pdb_ids = list(present_pdb_ids & pdb_ids_set)
        
        # Filter data
        self.data = self.data[self.data['pdb_id'].isin(self.pdb_ids)]
        
        # Update chains
        self.update_pdb_chain_ids()
    
    def filter_by_ids(self, pdb_ids: List[str]):
        """
        Filter data by PDB IDs (alias for filter_by_pdb_ids).
        
        Args:
            pdb_ids: List of PDB IDs to keep
        """
        self.filter_by_pdb_ids(pdb_ids)
    
    def filter_by_chain(self, chain_id):
        """
        Filter data to keep only specific chain.
        
        Args:
            chain_id: Chain ID to keep
        """
        if self.data is None:
            return
        
        self.data = self.data[self.data['auth_chain_id'] == chain_id].copy()
        self.reset_index()
    
    def filter_by_residue_range(self, start_res, end_res):
        """
        Filter data to keep only residues in range.
        
        Args:
            start_res: Start residue number
            end_res: End residue number
        """
        if self.data is None:
            return
        
        self.data = self.data[
            (self.data['auth_seq_id'] >= start_res) & 
            (self.data['auth_seq_id'] <= end_res)
        ].copy()
        self.reset_index()
    
    def remove_hetatm_records(self):
        """Remove HETATM records from data."""
        if self.data is None:
            return
        
        self.data = self.data[self.data['group'] == 'ATOM'].copy()
        self.reset_index()
    
    def filter_data_flexibly(self, filter_dict, inplace=True):
        """
        Filter data based on flexible criteria.
        
        Args:
            filter_dict: Dictionary of column: value pairs to filter on
            inplace: If True, modify self.data. If False, return filtered copy
            
        Returns:
            None if inplace=True, filtered DataFrame if inplace=False
        """
        if self.data is None or self.data.empty:
            if inplace:
                return
            else:
                return pd.DataFrame()
        
        # Start with a copy of data if not inplace
        if inplace:
            result = self.data
        else:
            result = self.data.copy()
            
        for key, value in filter_dict.items():
            # Handle operator suffixes like __eq, __ne, __in, etc.
            if '__' in key:
                column, operator = key.rsplit('__', 1)
            else:
                column = key
                operator = 'eq'
                
            if column not in result.columns:
                raise ValueError(f"Filter column '{column}' not found in data")
                
            # Apply filter based on operator
            if operator == 'eq':
                if isinstance(value, list):
                    result = result[result[column].isin(value)]
                else:
                    result = result[result[column] == value]
            elif operator == 'ne':
                result = result[result[column] != value]
            elif operator in ['in', 'is_in']:
                result = result[result[column].isin(value)]
            elif operator == 'nin':
                result = result[~result[column].isin(value)]
            elif operator == 'gt':
                result = result[result[column] > value]
            elif operator == 'gte':
                result = result[result[column] >= value]
            elif operator == 'lt':
                result = result[result[column] < value]
            elif operator == 'lte' or operator == 'le':
                result = result[result[column] <= value]
            elif operator == 'contains':
                result = result[result[column].str.contains(value, na=False)]
            elif operator == 'startswith':
                result = result[result[column].str.startswith(value, na=False)]
            elif operator == 'endswith':
                result = result[result[column].str.endswith(value, na=False)]
            elif operator == 'isna':
                result = result[result[column].isna()]
            elif operator == 'notna':
                result = result[result[column].notna()]
            else:
                raise ValueError(f"Invalid filter operator: __{operator}")
                
        if inplace:
            self.data = result
            # Update pdb_ids based on filtered data
            if 'pdb_id' in self.data.columns:
                self.pdb_ids = list(self.data['pdb_id'].unique())
            return None
        else:
            return result
    
    def get_chain(self, pdb_chain_id):
        """
        Get data for a specific chain.
        
        Args:
            pdb_chain_id: Combined PDB and chain ID (e.g., '1U19_A')
            
        Returns:
            DataFrame with chain data
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if '_' not in pdb_chain_id:
            raise ValueError(f"Invalid pdb_chain_id format: {pdb_chain_id}. Expected format: 'pdbid_chainid'")
            
        parts = pdb_chain_id.split('_')
        chain_id = parts[-1]
        pdb_id = '_'.join(parts[:-1])
        
        return self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id)]
    
    def get_backbone(self, pdb_chain_id):
        """
        Get backbone atoms for a chain.
        
        Args:
            pdb_chain_id: Combined PDB and chain ID
            
        Returns:
            DataFrame with backbone atoms
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if '_' not in pdb_chain_id:
            raise ValueError(f"Invalid pdb_chain_id format: {pdb_chain_id}. Expected format: 'pdbid_chainid'")
            
        parts = pdb_chain_id.split('_')
        chain_id = parts[-1]
        pdb_id = '_'.join(parts[:-1])
        
        return self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id) &
                        (self.data['res_atom_name'] == ALPHA_CARBON)]
    
    def apply_grn_interval(self, grn_interval):
        """
        Filter data to include only specified GRNs.
        
        Args:
            grn_interval: List of GRN values to include
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if 'grn' not in self.data.columns:
            self.logger.warning("Application of a GRN interval selection requires the 'grn' annotation")
            return
            
        self.data = self.data[self.data['grn'].isin(grn_interval)]
        self.update_pdb_ids()
    
    # GRN methods
    def assign_grns(self, grn_table, chain_ids=None, grn_chain_dict=None):
        """
        Assign GRN annotations to structure data.
        
        Args:
            grn_table: DataFrame with GRN mappings
            chain_ids: List of chain IDs to process
            grn_chain_dict: Dictionary mapping chain IDs to GRN chain IDs
        """
        if self.data is None:
            raise ValueError("No data loaded")
            
        # Implementation would go here
        self.logger.warning("GRN assignment not yet implemented in merged version")
    
    def get_grn_dict(self):
        """Get dictionary of GRN values."""
        if self.data is None or 'grn' not in self.data.columns:
            return {}
            
        grn_dict = {}
        for pdb_id in self.pdb_ids:
            chains = self.get_chains(pdb_id)
            for chain_id in chains:
                chain_data = self.data[(self.data['pdb_id'] == pdb_id) & 
                                      (self.data['auth_chain_id'] == chain_id)]
                if 'grn' in chain_data.columns:
                    grn_dict[f"{pdb_id}_{chain_id}"] = chain_data['grn'].unique().tolist()
        
        return grn_dict
    
    # Download methods
    def download_structure(self, pdb_id: str, save_to_cache: bool = True) -> Optional[pd.DataFrame]:
        """
        Download a structure from RCSB PDB.
        
        Args:
            pdb_id: PDB identifier
            save_to_cache: Whether to save to cache after download
            
        Returns:
            Structure DataFrame or None if download fails
        """
        # Download URL
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        
        try:
            # Download to temp file
            with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp_file:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                tmp_file.write(response.content)
                tmp_path = tmp_file.name
            
            # Parse the structure
            df = load_structure_util(
                tmp_path,
                remove_hetatm=self.remove_hetatm,
                allow_exception=self.allow_exception
            )
            
            # Clean up temp file
            os.unlink(tmp_path)
            
            # Save to mmcif directory
            cif_path = self.path_structure_dir / f"{pdb_id}.cif"
            # ProtosPaths handles directory creation
            
            with open(cif_path, 'wb') as f:
                f.write(response.content)
            
            # Save to cache if requested
            if save_to_cache and df is not None:
                self.save_entity(pdb_id, df)
            
            return df
            
        except Exception as e:
            self.logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    @staticmethod
    def download_cif(pdb_id, save_dir=None):
        """
        Download a CIF file from the PDB (static method).
        
        Args:
            pdb_id: PDB ID to download
            save_dir: Directory to save the file
            
        Returns:
            Path to saved file or None if failed
        """
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            if save_dir is None:
                save_dir = Path.cwd()
            else:
                save_dir = Path(save_dir)
                
            file_path = save_dir / f"{pdb_id}.cif"
            with open(file_path, 'wb') as f:
                f.write(response.content)
                
            return str(file_path)
        except Exception as e:
            logging.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    def check_and_download_missing_cifs(self, pdb_ids=None):
        """
        Check for missing CIF files and download them.
        
        Args:
            pdb_ids: List of PDB IDs to check (defaults to self.pdb_ids)
        """
        if pdb_ids is None:
            pdb_ids = self.pdb_ids
            
        missing_pdbs = []
        for pdb_id in pdb_ids:
            cif_path = self.path_structure_dir / f"{pdb_id}.cif"
            if not cif_path.exists():
                missing_pdbs.append(pdb_id)
        
        if missing_pdbs:
            self.logger.info(f"Downloading {len(missing_pdbs)} missing CIF files...")
            for pdb_id in tqdm(missing_pdbs, desc="Downloading CIF files"):
                self.download_structure(pdb_id)
    
    # Extract binding pocket
    def extract_binding_pocket(self, pdb_id, ligand='RET', distance=5.0):
        """
        Extract atoms in the binding pocket around a ligand.
        
        Args:
            pdb_id: PDB ID of the structure
            ligand: Ligand residue name
            distance: Maximum distance to include
            
        Returns:
            DataFrame with atoms in the binding pocket
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
        
        # Get ligand data
        ligand_data = self.data[(self.data['pdb_id'] == pdb_id) &
                               (self.data['res_name3l'] == ligand)]
        
        if ligand_data.empty:
            self.logger.warning(f"No ligand '{ligand}' found in {pdb_id}")
            return pd.DataFrame()
        
        # Get chain with ligand
        chain_id = ligand_data['auth_chain_id'].iloc[0]
        
        # Get ligand coordinates
        ligand_coords = ligand_data[['x', 'y', 'z']].values.astype(float)
        
        # Get protein coordinates
        protein_data = self.data[(self.data['pdb_id'] == pdb_id) &
                                (self.data['auth_chain_id'] == chain_id) &
                                (self.data['group'] == 'ATOM')]
        
        protein_coords = protein_data[['x', 'y', 'z']].values.astype(float)
        
        # Find atoms within distance of ligand
        binding_pocket_indices = []
        
        for ligand_coord in ligand_coords:
            # Calculate distances
            distances = np.sqrt(np.sum((protein_coords - ligand_coord) ** 2, axis=1))
            # Find atoms within distance
            close_indices = np.where(distances <= distance)[0]
            binding_pocket_indices.extend(close_indices)
        
        # Remove duplicates
        binding_pocket_indices = list(set(binding_pocket_indices))
        
        # Return the binding pocket atoms
        if binding_pocket_indices:
            return protein_data.iloc[binding_pocket_indices].copy()
        else:
            return pd.DataFrame()
    
    # Ligand methods
    def add_ligand(self, ligand_df):
        """
        Add ligand data to the structure.
        
        Args:
            ligand_df: DataFrame with ligand data
        """
        if self.data is None:
            self.data = ligand_df
        else:
            self.data = pd.concat([self.data, ligand_df], ignore_index=True)
        self.reset_index()
    
    # Export methods
    def save_cif(self, df, output_path):
        """
        Save DataFrame as CIF file.
        
        Args:
            df: Structure DataFrame
            output_path: Path to save CIF file
        """
        # Convert DataFrame to CIF format
        cif_utils.dataframe_to_cif(df, output_path)
    
    def to_cif(self, output_path=None):
        """
        Convert current data to CIF format.
        
        Args:
            output_path: Path to save CIF file (optional)
            
        Returns:
            CIF string if no output path specified
        """
        if self.data is None:
            raise ValueError("No data loaded")
            
        if output_path:
            self.save_cif(self.data, output_path)
        else:
            return cif_utils.dataframe_to_cif_string(self.data)
    
    def save_temp_cif(self, df, name="temp"):
        """
        Save DataFrame as temporary CIF file.
        
        Args:
            df: Structure DataFrame
            name: Name for temp file
            
        Returns:
            Path to temp file
        """
        temp_path = self.temp_dir / f"{name}_{int(time.time())}.cif"
        # ProtosPaths handles directory creation
        self.save_cif(df, str(temp_path))
        return str(temp_path)
    
    def export_structures(self, output_dir, format='cif'):
        """
        Export all loaded structures.
        
        Args:
            output_dir: Directory to export to
            format: Export format ('cif' or 'pkl')
        """
        output_dir = Path(output_dir)
        # User is responsible for creating output directory
        if not output_dir.exists():
            raise ValueError(f"Output directory does not exist: {output_dir}")
        
        for pdb_id in tqdm(self.pdb_ids, desc="Exporting structures"):
            df = self.data[self.data['pdb_id'] == pdb_id]
            
            if format == 'cif':
                output_path = output_dir / f"{pdb_id}.cif"
                self.save_cif(df, str(output_path))
            elif format == 'pkl':
                output_path = output_dir / f"{pdb_id}.pkl"
                df.to_pickle(output_path)
    
    def cleanup_temp_files(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            for file in self.temp_dir.glob("temp_*.cif"):
                try:
                    file.unlink()
                except:
                    pass
    
    # Dataset management methods (legacy compatibility)
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all available structure entities.
        
        Args:
            dataset: Optional dataset ID to filter by
            
        Returns:
            List of PDB IDs (not hash IDs!)
        """
        if dataset:
            # If dataset specified, return structures in that dataset only
            dataset_info = self.get_dataset_info(dataset)
            if dataset_info and "entities" in dataset_info:
                return dataset_info["entities"]
            return []
        
        # Otherwise, list all available structure files
        return self.get_available_pdb_files()
    
    def list_structures(self, dataset: Optional[str] = None) -> List[str]:
        """List all available structures (backward compatibility)."""
        return self.list_entities(dataset=dataset)
    
    def list_datasets_detailed(self):
        """List available structure datasets with detailed information."""
        datasets = []
        
        # First check legacy datasets.json
        datasets_json_path = self.path_dataset_dir / 'datasets.json'
        if datasets_json_path.exists():
            try:
                with open(datasets_json_path, 'r') as f:
                    datasets_data = json.load(f)
                    for dataset_id, pdb_ids in datasets_data.items():
                        datasets.append({
                            'id': dataset_id,
                            'type': 'structure_dataset',
                            'structure_count': len(pdb_ids),
                            'source': 'datasets.json'
                        })
            except:
                pass
        
        # Then check for individual dataset JSON files
        if self.path_dataset_dir.exists():
            for file in self.path_dataset_dir.glob("*.json"):
                if file.name != 'datasets.json':
                    dataset_id = file.stem
                    try:
                        with open(file, 'r') as f:
                            dataset_data = json.load(f)
                            datasets.append({
                                'id': dataset_id,
                                'type': 'structure_dataset',
                                'name': dataset_data.get('name', dataset_id),
                                'description': dataset_data.get('description', ''),
                                'structure_count': len(dataset_data.get('pdb_ids', [])),
                                'source': file.name
                            })
                    except:
                        pass
        
        return datasets
    
    # Use BaseProcessor's list_datasets which returns just names
    # def list_datasets(self):
    #     """List available structure datasets."""
    #     return super().list_datasets()
    
    def get_dataset(self, dataset_name):
        """
        Get dataset information without loading structures.
        
        Args:
            dataset_name: Name of the dataset
            
        Returns:
            Dataset object with metadata and PDB ID list
        """
        # Use the standardized dataset API
        if self.dataset_manager is not None:
            dataset = self.dataset_manager.load_dataset(dataset_name)
            if dataset is not None:
                # Store PDB IDs but don't load structures
                self.pdb_ids = dataset.get('entities', dataset.get('pdb_ids', []))
                return dataset
        
        # Legacy fallback
        datasets_json_path = self.path_dataset_dir / 'datasets.json'
        if datasets_json_path.exists():
            try:
                with open(datasets_json_path, 'r') as f:
                    datasets_data = json.load(f)
                    if dataset_name in datasets_data:
                        self.pdb_ids = datasets_data[dataset_name]
                        return {
                            'id': dataset_name,
                            'pdb_ids': self.pdb_ids,
                            'source': 'datasets.json'
                        }
            except:
                pass
        
        return None
    
    def load_dataset(self, dataset_id, apply_dtypes=True, debug=False):
        """
        Load a dataset and all its structures.
        
        Args:
            dataset_id: Dataset identifier
            apply_dtypes: Whether to apply proper data types
            debug: Whether to print debug info
            
        Returns:
            List of loaded PDB IDs
        """
        # Get dataset info
        dataset = self.get_dataset(dataset_id)
        if dataset:
            # Handle both 'entities' (new format) and 'pdb_ids' (legacy format)
            pdb_ids = dataset.get('entities', dataset.get('pdb_ids', dataset.get('content', [])))
            if pdb_ids:
                self.load_structures(pdb_ids, apply_dtypes=apply_dtypes, debug=debug)
                self.update_pdb_ids()
                return self.pdb_ids
        
        self.logger.error(f"Could not load dataset '{dataset_id}'")
        return []
    
    def create_dataset(self, dataset_id: str, pdb_ids: List[str] = None, 
                      metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Create a new dataset.
        
        Args:
            dataset_id: Unique identifier for the dataset
            pdb_ids: List of PDB IDs (defaults to current pdb_ids)
            metadata: Additional metadata
            
        Returns:
            Dataset ID
        """
        # Use current PDB IDs if not provided
        if pdb_ids is None:
            pdb_ids = self.pdb_ids
            
        # Ensure we have PDB IDs
        if not pdb_ids:
            self.logger.error("No PDB IDs provided or loaded to create dataset")
            return None
        
        # Create metadata
        if metadata is None:
            metadata = {}
            
        # Use the standardized dataset API if available
        if self.dataset_manager is not None:
            name = metadata.get('name', dataset_id)
            description = metadata.get('description', f'Structure dataset with {len(pdb_ids)} structures')
            
            # Call parent class method to avoid recursion
            # BaseProcessor expects (dataset_name, entity_names, metadata)
            return super().create_dataset(
                dataset_id,  # dataset_name
                pdb_ids,     # entity_names
                metadata
            )
        else:
            # Legacy approach - save to datasets.json
            datasets_json_path = self.path_dataset_dir / 'datasets.json'
            datasets = {}
            
            if datasets_json_path.exists():
                try:
                    with open(datasets_json_path, 'r') as f:
                        datasets = json.load(f)
                except:
                    pass
            
            datasets[dataset_id] = pdb_ids
            
            # ProtosPaths handles directory creation
            with open(datasets_json_path, 'w') as f:
                json.dump(datasets, f, indent=4)
            
            return dataset_id
    
    def save_dataset(self, dataset_id: str, pdb_ids: List[str] = None):
        """Save dataset (alias for create_dataset)."""
        return self.create_dataset(dataset_id, pdb_ids)
    
    def delete_dataset(self, dataset_id: str) -> bool:
        """
        Delete a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if deletion was successful
        """
        # Use the standardized dataset API if available
        if self.dataset_manager is not None:
            try:
                self.dataset_manager.delete_dataset(dataset_id)
                return True
            except Exception as e:
                self.logger.error(f"Failed to delete dataset '{dataset_id}': {e}")
                return False
            
        # Legacy approach
        datasets_json_path = self.path_dataset_dir / 'datasets.json'
        if datasets_json_path.exists():
            try:
                with open(datasets_json_path, 'r') as f:
                    datasets = json.load(f)
                
                if dataset_id in datasets:
                    del datasets[dataset_id]
                    
                    with open(datasets_json_path, 'w') as f:
                        json.dump(datasets, f, indent=4)
                    
                    self.logger.info(f"Dataset '{dataset_id}' deleted successfully")
                    return True
            except Exception as e:
                self.logger.error(f"Failed to delete dataset '{dataset_id}': {e}")
                
        return False
    
    def list_dataset_structures(self, dataset_name: str) -> List[str]:
        """List structures in a specific dataset."""
        dataset = self.get_dataset(dataset_name)
        if dataset:
            return dataset.get('pdb_ids', dataset.get('content', []))
        return []
    
    # Persistence methods using BaseProcessor
    def save_data(self, dataset_id, data=None, **kwargs):
        """
        Save structure data to a dataset PKL in structure_dataset/.
        
        Args:
            dataset_id: Dataset identifier
            data: Data to save (uses self.data if None)
            **kwargs: Additional parameters for saving
            
        Returns:
            Path to saved file
        """
        if data is None:
            if self.data is None:
                raise ValueError("No data to save")
            data = self.data
        
        # Set default file format
        file_format = kwargs.pop("file_format", "pkl")
        
        # Save to structure_dataset directory
        dataset_path = self.path_dataset_dir / f"{dataset_id}.{file_format}"
        
        # ProtosPaths ensures directory exists
        if file_format == 'pkl':
            data.to_pickle(dataset_path)
        elif file_format == 'csv':
            data.to_csv(dataset_path, index=False)
        else:
            raise ValueError(f"Unsupported format: {file_format}")
        
        self.logger.info(f"Saved dataset '{dataset_id}' to {dataset_path}")
        return str(dataset_path)
        
    def load_data(self, dataset_id, apply_dtypes=True, debug=False, **kwargs):
        """
        Load structure data from a dataset.
        
        Priority:
        1. Check for preprocessed PKL in structure_dataset/
        2. Fall back to loading individual structures
        
        Args:
            dataset_id: Dataset identifier
            apply_dtypes: Whether to apply proper data types
            debug: Whether to print debug info
            **kwargs: Additional parameters
            
        Returns:
            Loaded data
        """
        # First, try to load preprocessed dataset from structure_dataset/
        file_format = kwargs.get("file_format", "pkl")
        
        # Validate format early
        if file_format not in ['pkl', 'csv']:
            raise ValueError(f"Unsupported format: {file_format}")
        
        dataset_path = self.path_dataset_dir / f"{dataset_id}.{file_format}"
        
        if dataset_path.exists():
            self.logger.info(f"Loading preprocessed dataset from {dataset_path}")
            try:
                if file_format == 'pkl':
                    self.data = pd.read_pickle(dataset_path)
                elif file_format == 'csv':
                    self.data = pd.read_csv(dataset_path)
                else:
                    raise ValueError(f"Unsupported format: {file_format}")
                
                # Apply data types if requested
                if apply_dtypes and self.data is not None and not self.data.empty:
                    self.format_data_types(debug=debug)
                
                # Update state
                if self.data is not None:
                    self.update_pdb_ids()
                    self.dfl = [self.data[self.data['pdb_id'] == pdb_id] for pdb_id in self.pdb_ids]
                
                return self.data
            except Exception as e:
                self.logger.warning(f"Failed to load preprocessed dataset: {e}")
        
        # Fall back to loading individual structures
        self.logger.info(f"No preprocessed dataset found, loading individual structures for '{dataset_id}'")
        
        # Get dataset info to find PDB IDs
        dataset_info = self.get_dataset(dataset_id)
        if dataset_info:
            pdb_ids = dataset_info.get('pdb_ids', dataset_info.get('content', dataset_info.get('entities', [])))
            if pdb_ids:
                # Load individual structures
                self.load_structures(pdb_ids, apply_dtypes=apply_dtypes, debug=debug)
                return self.data
        
        self.logger.error(f"Could not load dataset '{dataset_id}'")
        return None
    
    # Entity management methods
    
    
    def load_structure_by_entity(self, entity_id: str, apply_dtypes: bool = True, debug: bool = False):
        """Load a structure by its entity ID (legacy compatibility)."""
        # In the new system, we use human-readable names directly
        return self.load_structure(entity_id, apply_dtypes=apply_dtypes, debug=debug)
    
    def get_entity_id_for_pdb(self, pdb_id: str) -> Optional[str]:
        """Get the entity ID for a PDB ID (legacy compatibility)."""
        # Generate hash-based entity ID for compatibility with tests
        return generate_entity_id(pdb_id)
    
    def add_ligand(self, target_pdb_id: str, ligand_df: pd.DataFrame, 
                   ligand_chain_id: str = 'Z', ligand_res_seq_id: int = 9001) -> None:
        """
        Add a ligand to an existing structure.
        
        Args:
            target_pdb_id: PDB ID of the structure to add ligand to
            ligand_df: DataFrame containing ligand atoms with columns:
                      [atom_name, res_name, x, y, z, element (optional)]
            ligand_chain_id: Chain ID to assign to the ligand (default: 'Z')
            ligand_res_seq_id: Residue sequence ID for the ligand (default: 9001)
            
        Raises:
            ValueError: If target PDB not loaded, ligand_df is invalid, or missing columns
        """
        # Validate target PDB is loaded
        if self.data is None or self.data.empty:
            raise ValueError("No structures loaded. Load structures first.")
            
        if target_pdb_id not in self.data['pdb_id'].unique():
            raise ValueError(f"Target PDB ID '{target_pdb_id}' not found in loaded structures")
            
        # Validate ligand DataFrame
        if ligand_df is None or ligand_df.empty:
            raise ValueError("Ligand DataFrame cannot be empty")
            
        # Check required columns
        required_cols = ['atom_name', 'res_name', 'x', 'y', 'z']
        missing_cols = [col for col in required_cols if col not in ligand_df.columns]
        if missing_cols:
            raise ValueError(f"Ligand DataFrame is missing required columns: {missing_cols}")
            
        # Get the target structure data
        target_data = self.data[self.data['pdb_id'] == target_pdb_id]
        max_atom_id = target_data['atom_id'].max()
        
        # Prepare ligand data
        ligand_data = ligand_df.copy()
        
        # Add missing columns with appropriate values
        ligand_data['pdb_id'] = target_pdb_id
        ligand_data['auth_chain_id'] = ligand_chain_id
        ligand_data['auth_seq_id'] = ligand_res_seq_id
        ligand_data['group'] = 'HETATM'
        ligand_data['model_num'] = 0
        
        # Assign new atom IDs
        num_ligand_atoms = len(ligand_data)
        ligand_data['atom_id'] = range(max_atom_id + 1, max_atom_id + num_ligand_atoms + 1)
        
        # Add element column if not present
        if 'element' not in ligand_data.columns:
            # Infer element from atom name (first letter usually)
            ligand_data['element'] = ligand_data['atom_name'].str[0]
            
        # Add other columns that might be expected
        if 'b_factor' not in ligand_data.columns:
            ligand_data['b_factor'] = 50.0  # Default B-factor
            
        # Ensure all expected columns are present (copy from first row of target)
        first_row = target_data.iloc[0]
        for col in target_data.columns:
            if col not in ligand_data.columns:
                # Use appropriate default values
                if col in ['label_atom_id', 'label_comp_id']:
                    ligand_data[col] = ligand_data.get('atom_name', '')
                elif col == 'res_name3l':
                    ligand_data[col] = ligand_data.get('res_name', 'LIG')
                elif col in ['entity_id', 'label_entity_id']:
                    ligand_data[col] = 999  # Ligand entity
                elif col == 'occupancy':
                    ligand_data[col] = 1.0
                elif col in ['alt_id', 'label_alt_id']:
                    ligand_data[col] = '.'
                else:
                    # For other columns, use NaN or empty string
                    if first_row[col] is not None and isinstance(first_row[col], str):
                        ligand_data[col] = ''
                    else:
                        ligand_data[col] = np.nan
        
        # Append ligand data to main data
        self.data = pd.concat([self.data, ligand_data], ignore_index=True)
        
        # Sort by PDB ID and atom ID for consistency
        self.data = self.data.sort_values(['pdb_id', 'atom_id']).reset_index(drop=True)
    
    def save_structure_as_entity(self, structure_df: pd.DataFrame, pdb_id: str, 
                                datasets: Optional[List[str]] = None,
                                metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Save a structure DataFrame as an entity.
        
        Args:
            structure_df: DataFrame containing structure data
            pdb_id: PDB ID for the structure
            datasets: List of dataset IDs to associate with this entity
            metadata: Additional metadata for the entity
            
        Returns:
            Entity ID of the saved structure
        """
        # Use ProtosPaths to construct the proper path
        # Following the principle: filepath = f(ProtosPaths, human readable name)
        cif_path = self.path_structure_dir / f"{pdb_id}.cif"
        
        # Save the CIF file using cif_utils
        from protos.io import cif_utils
        cif_utils.write_cif_file(str(cif_path), structure_df, force_overwrite=True)
        
        # Generate entity ID
        entity_id = generate_entity_id(pdb_id)
        
        # Register the entity
        self._register_structure_entity(pdb_id, structure_df, datasets=datasets, metadata=metadata)
        
        return entity_id
    
    def _register_structure_entity(self, pdb_id: str, structure_df: pd.DataFrame,
                                  datasets: Optional[List[str]] = None,
                                  metadata: Optional[Dict[str, Any]] = None) -> str:
        """Register a structure as an entity in the registry."""
        entity_id = generate_entity_id(pdb_id)
        
        # Prepare metadata
        if metadata is None:
            metadata = {}
            
        # Add structure-specific metadata
        chains = structure_df['auth_chain_id'].unique().tolist() if 'auth_chain_id' in structure_df.columns else []
        metadata.update({
            'atom_count': len(structure_df),
            'chains': chains,
            'pdb_id': pdb_id
        })
        
        # Register with entity registry
        # Following the principle: filepath = f(ProtosPaths, human readable name)
        if hasattr(self, 'entity_registry') and self.entity_registry is not None:
            # Construct the path using ProtosPaths
            cif_path = self.path_structure_dir / f"{pdb_id}.cif"
            # Convert to relative path for registry storage
            relative_path = cif_path.relative_to(self.paths.data_root) if hasattr(self.paths, 'data_root') else cif_path
            
            # Add datasets to metadata
            if datasets:
                metadata['datasets'] = datasets
                
            self.entity_registry.register_entity(
                pdb_id,  # Use human-readable name, not hash
                "structure",
                str(relative_path),
                metadata
            )
        
        return pdb_id  # Return human-readable name, not hash
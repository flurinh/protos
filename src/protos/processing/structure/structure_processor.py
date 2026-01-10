from __future__ import annotations

import hashlib
import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple, TYPE_CHECKING, Iterable, Set

import pandas as pd
import numpy as np

from protos.io.core.base_processor import BaseProcessor
from protos.io.formats.structure_schema import STRUCT_COLUMN_DTYPE, SORTED_STRUCT_COLUMNS
from protos.analysis.structure_water_networks import analyze_water_networks

if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from .alignment_engine import StructureAlignmentEngine, AlignmentResult


class StructureProcessor(BaseProcessor):
    """
    Structure processor with PKL-canonical architecture.

    Key principles:
    - PKL is the canonical on-disk format handled by the processor
    - Per-entity storage with MultiIndex (structure_id, atom_id)
    - Lazy stacking with frames dict and dirty flag
    - All paths via ProtosPaths, no hardcoding
    """
    
    # Define processor type
    processor_type = "structure"
    
    def __init__(self, name: str = "structure_processor"):
        """
        Initialize StructureProcessor - NO PATH PARAMETERS!
        
        Args:
            name: Processor instance name
        """
        super().__init__(name=name)
        
        # Storage: dict of DataFrames with lazy stacking
        self.frames: Dict[str, pd.DataFrame] = {}
        self._data: Optional[pd.DataFrame] = None
        self._dirty: bool = False
        self._alignment_engine: Optional["StructureAlignmentEngine"] = None
        self._sequence_loader = None
        # Legacy caches for compatibility with pre-unified workflows
        self.chain_dict: Dict[str, str] = {}

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
    def path_pdb_dir(self) -> Path:
        """Get PDB directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'pdb_dir'))

    @property
    def path_sdf_dir(self) -> Path:
        """Get SDF directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'sdf_dir'))

    @property
    def path_dataset_dir(self) -> Path:
        """Get dataset directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'dataset_dir'))
    
    @property
    def path_temp_dir(self) -> Path:
        """Get temp directory from ProtosPaths"""
        return Path(self.paths.get_subdir_path('structure', 'temp_dir'))

    def _get_sequence_loader(self):
        if self._sequence_loader is None:
            from protos.io.ingest.sequence_loader import SequenceLoader

            self._sequence_loader = SequenceLoader(name=f"{self.name}_chain_loader")
        return self._sequence_loader

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
        
        # If DataFrame already stored in canonical form (MultiIndex), expose index as columns
        if isinstance(df.index, pd.MultiIndex):
            index_names = set(name or '' for name in df.index.names)
            if {'structure_id', 'atom_id'}.issubset(index_names):
                df = df.reset_index()

        # Map pdb_id to structure_id if needed
        if 'pdb_id' in df.columns and 'structure_id' not in df.columns:
            df['structure_id'] = df['pdb_id']

        # Ensure structure_id column
        if 'structure_id' not in df.columns:
            df['structure_id'] = structure_id

        df['structure_id'] = df['structure_id'].fillna(structure_id).astype(str)

        # Drop legacy columns we no longer track in canonical form
        legacy_drop = {'pdb_id'}
        existing_drop = [col for col in legacy_drop if col in df.columns]
        if existing_drop:
            df = df.drop(columns=existing_drop)

        # Ensure atom_id column exists and is numeric
        if 'atom_id' not in df.columns:
            raise ValueError("DataFrame must contain 'atom_id' column")

        atom_ids = pd.to_numeric(df['atom_id'], errors='coerce')
        if atom_ids.isna().any():
            raise ValueError("'atom_id' column contains non-numeric values")
        df['atom_id'] = atom_ids.astype('Int64')

        # Ensure all expected columns exist with sensible defaults
        for column, dtype in STRUCT_COLUMN_DTYPE.items():
            if column in ('structure_id', 'atom_id'):
                continue
            if column not in df.columns:
                if dtype in (int, 'int'):
                    df[column] = pd.Series(pd.NA, index=df.index, dtype='Int64')
                elif dtype is float:
                    df[column] = np.nan
                else:
                    df[column] = ''

        # Apply dtypes consistently
        for column, dtype in STRUCT_COLUMN_DTYPE.items():
            if column not in df.columns:
                continue
            try:
                if dtype in (int, 'int'):
                    df[column] = pd.to_numeric(df[column], errors='coerce').astype('Int64')
                elif dtype is float:
                    df[column] = pd.to_numeric(df[column], errors='coerce')
                else:
                    df[column] = df[column].fillna('').astype(str)
            except Exception:
                # Leave column as-is if conversion fails; data will be handled upstream
                pass

        # Default values for optional textual columns
        if 'grn' in df.columns:
            df['grn'] = df['grn'].fillna('')

        # Reorder columns to canonical order before setting index
        base_cols = ['structure_id', 'atom_id']
        ordered_cols = [col for col in SORTED_STRUCT_COLUMNS if col not in base_cols and col in df.columns]
        remaining_cols = [col for col in df.columns if col not in base_cols and col not in ordered_cols]
        column_order = [col for col in base_cols if col in df.columns] + ordered_cols + remaining_cols
        df = df[column_order]

        # Set MultiIndex and ensure deterministic ordering
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
    
    def load_entity(self, structure_id: str) -> Optional[pd.DataFrame]:
        """
        Load a structure exclusively from registered PKL storage.

        1. Check if already in memory
        2. Try registry lookup (format_type="structure")
        3. Try PKL in cache dir (auto-register if found)
        4. Return None if nothing found
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
            raise ValueError("Only 'pkl' format is allowed for save_entity.")
        
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

    def export_entity(
        self,
        name: str,
        out_path: Optional[Path] = None,
        format: Optional[str] = None,
        overwrite: bool = False,
        **kwargs,
    ) -> Path:
        exporter = self._get_exporter()
        return exporter.export_entity(
            name,
            out_path,
            format=format,
            overwrite=overwrite,
            **kwargs,
        )

    def export_dataset(
        self,
        dataset_name: str,
        output_dir: Optional[Path] = None,
        format: Optional[str] = None,
        overwrite: bool = False,
        name_pattern: Optional[str] = None,
        **kwargs,
    ) -> Dict[str, Path]:
        exporter = self._get_exporter()
        return exporter.export_dataset(
            dataset_name,
            output_dir,
            format=format,
            overwrite=overwrite,
            name_pattern=name_pattern,
            **kwargs,
        )

    def summarize_ligands(
        self,
        structure_id: str,
        *,
        group_by: Optional[Union[str, List[str]]] = None,
        include_chains: Optional[List[str]] = None,
        include_comp_ids: Optional[List[str]] = None,
        min_atoms: int = 1,
    ) -> Dict[str, Any]:
        """Return a nested dict summary of ligands for a structure.

        - Filters HETATM only
        - Groups by 'res_id' when available (default), otherwise by provided columns
        - Returns counts by chain and per ligand
        """
        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        if 'group' in df.columns:
            df = df[df['group'].str.upper() == 'HETATM']
        if df.empty:
            return {"entity": structure_id, "total_groups": 0, "chains": {}}

        if include_chains and 'auth_chain_id' in df.columns:
            df = df[df['auth_chain_id'].isin(include_chains)]
        if include_comp_ids and 'auth_comp_id' in df.columns:
            df = df[df['auth_comp_id'].isin(include_comp_ids)]

        # Determine grouping
        if group_by is None:
            if 'res_id' in df.columns:
                group_cols = ['res_id']
            else:
                group_cols = [c for c in ['auth_chain_id', 'auth_seq_id', 'auth_comp_id', 'insertion'] if c in df.columns]
                if not group_cols:
                    group_cols = ['auth_comp_id'] if 'auth_comp_id' in df.columns else []
        else:
            group_cols = [group_by] if isinstance(group_by, str) else list(group_by)

        if group_cols:
            grouped = df.groupby(group_cols, dropna=False)
        else:
            grouped = [(structure_id, df)]

        chains: Dict[str, Dict[str, Any]] = {}
        total = 0
        for key, sub in grouped:
            if len(sub) < min_atoms:
                continue
            total += 1
            # Resolve chain/comp/seq
            chain = str(sub.get('auth_chain_id', pd.Series(['?'])).iloc[0]) if 'auth_chain_id' in sub.columns else '?'
            comp = str(sub.get('auth_comp_id', pd.Series(['?'])).iloc[0]) if 'auth_comp_id' in sub.columns else '?'
            seq = sub.get('auth_seq_id', pd.Series([None])).iloc[0] if 'auth_seq_id' in sub.columns else None
            ins = sub.get('insertion', pd.Series([''])).iloc[0] if 'insertion' in sub.columns else ''
            res_id = sub.get('res_id', pd.Series([None])).iloc[0] if 'res_id' in sub.columns else None
            if not res_id:
                res_id = f"{comp}_{seq if seq is not None else ''}{ins}".strip('_')

            chains.setdefault(chain, {})[res_id] = {
                'comp_id': comp,
                'seq_id': int(seq) if seq is not None and pd.notna(seq) else None,
                'insertion': ins if ins and ins == ins else '',
                'atom_count': int(len(sub)),
            }

        return {
            'entity': structure_id,
            'grouping': group_cols,
            'total_groups': total,
            'chains': chains,
        }

    # ---------- Ligand & Contact Analysis ----------

    # Common hetero codes to exclude when listing ligands
    _COMMON_HETERO = {'HOH', 'WAT', 'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'SO4', 'PO4'}
    _WATER_CODES = {'HOH', 'WAT'}
    _ION_CODES = {'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'FE', 'CU', 'MN', 'CO', 'NI'}

    def list_ligands(
        self,
        structure_id: str,
        *,
        exclude_common: bool = True,
        min_atoms: int = 5,
    ) -> List[Dict[str, Any]]:
        """
        List all ligands in a structure.

        Args:
            structure_id: Structure identifier
            exclude_common: Exclude water, common ions, and buffer molecules
            min_atoms: Minimum number of atoms for a molecule to be considered a ligand

        Returns:
            List of ligand dictionaries with keys:
            - res_name: 3-letter ligand code
            - chain_id: Chain containing the ligand
            - res_id: Residue sequence ID
            - num_atoms: Number of atoms in the ligand
            - centroid: (x, y, z) coordinates of ligand center
        """
        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        # Filter to HETATM only
        if 'group' not in df.columns:
            return []
        hetatoms = df[df['group'].str.upper() == 'HETATM']
        if hetatoms.empty:
            return []

        ligands = []
        # Group by residue
        group_cols = ['res_name3l', 'auth_chain_id', 'auth_seq_id']
        group_cols = [c for c in group_cols if c in hetatoms.columns]
        if not group_cols:
            return []

        for keys, atoms in hetatoms.groupby(group_cols, dropna=False):
            res_name = keys[0] if len(keys) > 0 else 'UNK'
            chain_id = keys[1] if len(keys) > 1 else '?'
            res_id = keys[2] if len(keys) > 2 else None

            # Skip common molecules if requested
            if exclude_common and res_name in self._COMMON_HETERO:
                continue

            # Skip small molecules
            if len(atoms) < min_atoms:
                continue

            # Calculate centroid
            coords = atoms[['x', 'y', 'z']].values
            centroid = tuple(coords.mean(axis=0))

            ligands.append({
                'res_name': res_name,
                'chain_id': chain_id,
                'res_id': int(res_id) if pd.notna(res_id) else None,
                'num_atoms': len(atoms),
                'centroid': centroid,
            })

        return ligands

    def get_ligand_interactions(
        self,
        structure_id: str,
        ligand_id: Optional[str] = None,
        *,
        chain_id: Optional[str] = None,
        cutoff: float = 4.5,
        include_water_bridges: bool = False,
    ) -> Dict[str, Any]:
        """
        Get protein residues that interact with a ligand.

        Args:
            structure_id: Structure identifier
            ligand_id: 3-letter ligand code (if None, uses largest ligand)
            chain_id: Specific chain for the ligand (optional)
            cutoff: Distance cutoff in Angstroms for binding residues
            include_water_bridges: Include water-mediated contacts

        Returns:
            Dictionary with:
            - ligand: ligand info dict
            - binding_residues: list of residue dicts with chain_id, res_id, res_name, min_distance, grn (if available)
            - summary: counts of different interaction types
        """
        from .ligand_interactions import LigandInteractionAnalyzer

        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        # Initialize analyzer
        analyzer = LigandInteractionAnalyzer(df)

        # Find ligand(s)
        ligands = analyzer.extract_ligands(exclude_common=True)
        if not ligands:
            return {'ligand': None, 'binding_residues': [], 'summary': {}}

        # Select specific ligand or largest
        if ligand_id:
            matching = [l for l in ligands if l['res_name3l'] == ligand_id]
            if chain_id:
                matching = [l for l in matching if l['chain_id'] == chain_id]
            if not matching:
                return {'ligand': None, 'binding_residues': [], 'summary': {}}
            ligand = matching[0]
        else:
            ligand = max(ligands, key=lambda x: x.get('num_atoms', 0))

        ligand_atoms = ligand['atoms']

        # Get binding residues
        binding_df = analyzer.get_binding_site_residues(ligand_atoms, cutoff=cutoff)

        # Build result
        binding_residues = []
        if not binding_df.empty:
            # Check if structure has GRN annotations
            has_grn = 'grn' in df.columns

            for _, row in binding_df.iterrows():
                res_info = {
                    'chain_id': row.get('chain_id'),
                    'res_id': int(row['res_id']) if pd.notna(row.get('res_id')) else None,
                    'res_name': row.get('res_name'),
                    'min_distance': float(row.get('min_distance', 0)),
                    'num_contacts': int(row.get('num_contacts', 1)),
                }

                # Add GRN if available
                if has_grn and res_info['chain_id'] and res_info['res_id']:
                    grn_mask = (
                        (df['auth_chain_id'] == res_info['chain_id']) &
                        (df['auth_seq_id'] == res_info['res_id'])
                    )
                    grn_vals = df.loc[grn_mask, 'grn'].dropna().unique()
                    if len(grn_vals) > 0 and grn_vals[0]:
                        res_info['grn'] = str(grn_vals[0])

                binding_residues.append(res_info)

        # Get interaction summary
        summary = {
            'num_binding_residues': len(binding_residues),
        }

        # Optionally include water bridges
        if include_water_bridges:
            water_bridges = analyzer.get_water_mediated_contacts(ligand_atoms)
            summary['num_water_bridges'] = len(water_bridges)

        return {
            'ligand': {
                'res_name': ligand['res_name3l'],
                'chain_id': ligand['chain_id'],
                'res_id': ligand['res_id'],
                'num_atoms': ligand['num_atoms'],
            },
            'binding_residues': binding_residues,
            'summary': summary,
        }

    def get_water_contacts(
        self,
        structure_id: str,
        *,
        cutoff: float = 3.5,
        protein_chain: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Get protein residues in contact with water molecules.

        Args:
            structure_id: Structure identifier
            cutoff: Distance cutoff in Angstroms
            protein_chain: Filter to specific protein chain

        Returns:
            List of contact dicts with water_id, protein residue info, and distance
        """
        from scipy.spatial import cKDTree

        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        # Get water molecules
        waters = df[
            (df['group'].str.upper() == 'HETATM') &
            (df['res_name3l'].isin(self._WATER_CODES))
        ]
        if waters.empty:
            return []

        # Get protein atoms
        protein = df[df['group'].str.upper() == 'ATOM']
        if protein_chain:
            protein = protein[protein['auth_chain_id'] == protein_chain]
        if protein.empty:
            return []

        # Build KD-tree for protein atoms
        protein_coords = protein[['x', 'y', 'z']].values
        protein_tree = cKDTree(protein_coords)

        contacts = []
        has_grn = 'grn' in df.columns

        # Group waters by residue
        water_groups = waters.groupby(['auth_chain_id', 'auth_seq_id'], dropna=False)

        for (water_chain, water_seq), water_atoms in water_groups:
            water_coords = water_atoms[['x', 'y', 'z']].values

            # Find nearby protein atoms
            nearby_indices = set()
            for coord in water_coords:
                indices = protein_tree.query_ball_point(coord, cutoff)
                nearby_indices.update(indices)

            if not nearby_indices:
                continue

            # Get unique residues
            nearby_protein = protein.iloc[list(nearby_indices)]
            for (chain, res_id, res_name), res_atoms in nearby_protein.groupby(
                ['auth_chain_id', 'auth_seq_id', 'res_name3l'], dropna=False
            ):
                # Calculate minimum distance
                res_coords = res_atoms[['x', 'y', 'z']].values
                from scipy.spatial import distance_matrix
                dists = distance_matrix(water_coords, res_coords)
                min_dist = float(dists.min())

                contact = {
                    'water_chain': water_chain,
                    'water_seq': int(water_seq) if pd.notna(water_seq) else None,
                    'protein_chain': chain,
                    'protein_res_id': int(res_id) if pd.notna(res_id) else None,
                    'protein_res_name': res_name,
                    'min_distance': min_dist,
                }

                # Add GRN if available
                if has_grn:
                    grn_mask = (
                        (df['auth_chain_id'] == chain) &
                        (df['auth_seq_id'] == res_id)
                    )
                    grn_vals = df.loc[grn_mask, 'grn'].dropna().unique()
                    if len(grn_vals) > 0 and grn_vals[0]:
                        contact['grn'] = str(grn_vals[0])

                contacts.append(contact)

        return contacts

    def get_ion_contacts(
        self,
        structure_id: str,
        *,
        cutoff: float = 3.5,
        protein_chain: Optional[str] = None,
        ion_types: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Get protein residues in contact with ions.

        Args:
            structure_id: Structure identifier
            cutoff: Distance cutoff in Angstroms
            protein_chain: Filter to specific protein chain
            ion_types: List of ion codes to include (default: common ions)

        Returns:
            List of contact dicts with ion info, protein residue info, and distance
        """
        from scipy.spatial import cKDTree

        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        # Get ions
        ion_codes = set(ion_types) if ion_types else self._ION_CODES
        ions = df[
            (df['group'].str.upper() == 'HETATM') &
            (df['res_name3l'].isin(ion_codes))
        ]
        if ions.empty:
            return []

        # Get protein atoms
        protein = df[df['group'].str.upper() == 'ATOM']
        if protein_chain:
            protein = protein[protein['auth_chain_id'] == protein_chain]
        if protein.empty:
            return []

        # Build KD-tree for protein atoms
        protein_coords = protein[['x', 'y', 'z']].values
        protein_tree = cKDTree(protein_coords)

        contacts = []
        has_grn = 'grn' in df.columns

        for _, ion in ions.iterrows():
            ion_coord = ion[['x', 'y', 'z']].values.reshape(1, -1)

            # Find nearby protein atoms
            indices = protein_tree.query_ball_point(ion_coord[0], cutoff)
            if not indices:
                continue

            # Get unique residues
            nearby_protein = protein.iloc[indices]
            for (chain, res_id, res_name), res_atoms in nearby_protein.groupby(
                ['auth_chain_id', 'auth_seq_id', 'res_name3l'], dropna=False
            ):
                # Calculate minimum distance
                res_coords = res_atoms[['x', 'y', 'z']].values
                from scipy.spatial import distance_matrix
                dists = distance_matrix(ion_coord, res_coords)
                min_dist = float(dists.min())

                contact = {
                    'ion_type': ion['res_name3l'],
                    'ion_chain': ion['auth_chain_id'],
                    'ion_seq': int(ion['auth_seq_id']) if pd.notna(ion['auth_seq_id']) else None,
                    'protein_chain': chain,
                    'protein_res_id': int(res_id) if pd.notna(res_id) else None,
                    'protein_res_name': res_name,
                    'min_distance': min_dist,
                }

                # Add GRN if available
                if has_grn:
                    grn_mask = (
                        (df['auth_chain_id'] == chain) &
                        (df['auth_seq_id'] == res_id)
                    )
                    grn_vals = df.loc[grn_mask, 'grn'].dropna().unique()
                    if len(grn_vals) > 0 and grn_vals[0]:
                        contact['grn'] = str(grn_vals[0])

                contacts.append(contact)

        return contacts

    def annotate_with_grn(
        self,
        structure_id: str,
        *,
        reference_table: str = "gpcrdb_ref",
        protein_family: str = "gpcr_a",
        chains: Optional[List[str]] = None,
        save: bool = True,
    ) -> pd.DataFrame:
        """
        Directly annotate a structure with Generic Residue Numbers (GRN).

        This method handles the full workflow internally:
        1. Extracts chain sequences from the structure
        2. Aligns sequences to GRN reference
        3. Maps GRN positions back to structure residues

        Args:
            structure_id: Structure identifier
            reference_table: GRN reference table name
            protein_family: Protein family for GRN reference
            chains: Specific chains to annotate (default: all protein chains)
            save: Whether to save the annotated structure

        Returns:
            Structure DataFrame with 'grn' column populated
        """
        from protos.processing.grn import GRNProcessor

        df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' not found")
        df = df.reset_index()

        # Initialize GRN column if needed
        if 'grn' not in df.columns:
            df['grn'] = ''

        # Get protein atoms only
        protein = df[df['group'].str.upper() == 'ATOM']
        if protein.empty:
            self.logger.warning(f"No protein atoms in {structure_id}")
            return self._ensure_canonical(df, structure_id)

        # Filter chains if specified
        available_chains = protein['auth_chain_id'].unique().tolist()
        target_chains = chains if chains else available_chains

        grn_proc = GRNProcessor()

        for chain_id in target_chains:
            if chain_id not in available_chains:
                continue

            # Extract sequence for this chain
            chain_atoms = protein[protein['auth_chain_id'] == chain_id]

            # Get unique residues in order
            residue_data = []
            for (res_id, res_name), atoms in chain_atoms.groupby(
                ['auth_seq_id', 'res_name3l'], sort=True, dropna=False
            ):
                if pd.isna(res_id):
                    continue
                residue_data.append((int(res_id), res_name))

            if not residue_data:
                continue

            # Convert to sequence (standard amino acid mapping)
            THREE_TO_ONE = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            }
            sequence = ''
            res_id_list = []
            for res_id, res_name in sorted(residue_data, key=lambda x: x[0]):
                aa = THREE_TO_ONE.get(res_name, 'X')
                sequence += aa
                res_id_list.append(res_id)

            if not sequence or len(sequence) < 10:
                continue

            # Annotate sequence with GRN
            seq_name = f"{structure_id}_{chain_id}"
            try:
                annotations_df, summary = grn_proc.annotate_sequences(
                    {seq_name: sequence},
                    reference_table=reference_table,
                    protein_family=protein_family,
                )

                if seq_name not in annotations_df.index:
                    continue

                # Get the annotation row for this sequence
                # Values are like 'M1', 'E2' where the number is 1-indexed position in sequence
                row = annotations_df.loc[seq_name]

                # Map GRN annotations back to structure residues
                import re
                for grn_position, value in row.items():
                    if value == '-' or not value:
                        continue
                    # Parse position from value like 'M1', 'E2', etc.
                    match = re.match(r'([A-Z])(\d+)', str(value))
                    if not match:
                        continue
                    seq_pos = int(match.group(2))  # 1-indexed position in sequence
                    if seq_pos < 1 or seq_pos > len(res_id_list):
                        continue
                    res_id = res_id_list[seq_pos - 1]  # Convert to 0-indexed

                    mask = (
                        (df['auth_chain_id'] == chain_id) &
                        (df['auth_seq_id'] == res_id)
                    )
                    df.loc[mask, 'grn'] = grn_position

            except Exception as e:
                self.logger.warning(f"GRN annotation failed for {structure_id} chain {chain_id}: {e}")
                continue

        # Canonicalize and optionally save
        df_canonical = self._ensure_canonical(df, structure_id)
        self._set_frame(structure_id, df_canonical)

        if save:
            self.save_entity(structure_id, df_canonical)

        return df_canonical

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

    def get_dataset(self, dataset_name: str) -> Dict[str, Any]:
        """Return the raw dataset definition stored on disk."""
        return self.dataset_manager.load_dataset(dataset_name)

    def get_dataset_entities(self, dataset_name: str) -> List[str]:
        """Expose dataset entity resolution with current human-readable names."""
        return self.dataset_manager.get_dataset_entities(dataset_name)

    def load_dataset(self, dataset_name: str, return_format: str = 'stacked') -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
        """
        Load all dataset members from individual PKLs.

        - Get member list from dataset
        - Load each via load_entity()
        - Update frames in bulk
        - Return stacked or dict format
        """
        structure_ids = self.get_dataset_entities(dataset_name)

        # Load all structures
        for structure_id in structure_ids:
            self.load_entity(structure_id)

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
        
        # Create logical dataset definition via DatasetManager
        self.create_dataset(dataset_name, structure_ids, metadata or {})

    def delete_entity(self, name: str) -> bool:
        """Delete entity and clear any cached frame."""
        removed = super().delete_entity(name)
        if removed:
            self._remove_frame(name)
        return removed

    # ---------- Utility Methods ----------

    def _get_exporter(self):
        from protos.io.export.structure_exporter import StructureExporter

        return StructureExporter(self)
    
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

    @property
    def alignment_engine(self) -> "StructureAlignmentEngine":
        """Lazily instantiate the structure alignment engine."""

        if self._alignment_engine is None:
            from .alignment_engine import StructureAlignmentEngine

            self._alignment_engine = StructureAlignmentEngine(self)
        return self._alignment_engine

    # ---------- Data Modification Methods ----------
    
    def delete_atoms(self, structure_id: str, atom_ids: List[int]) -> pd.DataFrame:
        """
        Delete specific atoms from structure.
        
        Args:
            structure_id: Structure identifier
            atom_ids: List of atom IDs to delete
            
        Returns:
            Modified structure DataFrame
        """
        # Get structure from frames
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Delete atoms by filtering out the specified atom_ids
        df = df[~df.index.get_level_values('atom_id').isin(atom_ids)]
        
        # Update frames
        self._set_frame(structure_id, df)
        
        return df
    
    def delete_residues(self, structure_id: str, chain_id: str, residue_ids: List[int]) -> pd.DataFrame:
        """
        Delete specific residues from structure.
        
        Args:
            structure_id: Structure identifier
            chain_id: Chain identifier
            residue_ids: List of residue IDs (auth_seq_id) to delete
            
        Returns:
            Modified structure DataFrame
        """
        # Get structure
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Reset index to access columns
        df_reset = df.reset_index()
        
        # Filter out specified residues
        mask = ~((df_reset['auth_chain_id'] == chain_id) & 
                 (df_reset['auth_seq_id'].isin(residue_ids)))
        df_filtered = df_reset[mask]
        
        # Re-canonicalize
        df_canonical = self._ensure_canonical(df_filtered, structure_id)
        
        # Update frames
        self._set_frame(structure_id, df_canonical)
        
        return df_canonical
    
    def delete_chain(self, structure_id: str, chain_id: str) -> pd.DataFrame:
        """
        Delete entire chain from structure.
        
        Args:
            structure_id: Structure identifier
            chain_id: Chain ID to delete
            
        Returns:
            Modified structure DataFrame
        """
        # Get structure
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Reset index to access chain column
        df_reset = df.reset_index()
        
        # Filter out the specified chain
        df_filtered = df_reset[df_reset['auth_chain_id'] != chain_id]
        
        # Re-canonicalize
        df_canonical = self._ensure_canonical(df_filtered, structure_id)
        
        # Update frames
        self._set_frame(structure_id, df_canonical)
        
        return df_canonical
    
    def filter_by_chain(
        self,
        structure_id: str,
        chain_ids: List[str],
        new_id: Optional[str] = None,
        *,
        register: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> pd.DataFrame:
        """
        Filter structure to keep only specified chains.
        Creates a new structure entity.
        
        Args:
            structure_id: Source structure identifier
            chain_ids: List of chain IDs to keep
            new_id: ID for filtered structure (auto-generated if not provided)
            
        Returns:
            Filtered structure DataFrame
        """
        # Get structure
        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        
        # Reset index to filter
        df_reset = df.reset_index()
        
        # Filter by chains
        df_filtered = df_reset[df_reset['auth_chain_id'].isin(chain_ids)]
        
        if df_filtered.empty:
            raise ValueError(f"No atoms found for chains {chain_ids}")
        
        # Create new structure ID if not provided
        if new_id is None:
            new_id = f"{structure_id}_chains_{''.join(sorted(chain_ids))}"
        
        # Ensure canonical form with new ID
        df_canonical = self._ensure_canonical(df_filtered, new_id)
        
        # Save as new entity
        self._set_frame(new_id, df_canonical)

        if register:
            self.save_entity(new_id, df_canonical, metadata=metadata)

        return df_canonical

    def filter_by_residue_range(
        self,
        structure_id: str,
        chain_id: str,
        start: int,
        end: int,
        new_id: Optional[str] = None,
        *,
        register: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> pd.DataFrame:
        """
        Filter structure to residue range.
        Creates a new structure entity.
        
        Args:
            structure_id: Source structure identifier
            chain_id: Chain to filter
            start: Start residue number (auth_seq_id)
            end: End residue number (auth_seq_id)
            new_id: ID for filtered structure
            
        Returns:
            Filtered structure DataFrame
        """
        # Get structure
        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        
        # Reset index to filter
        df_reset = df.reset_index()
        
        # Filter by residue range
        mask = ((df_reset['auth_chain_id'] == chain_id) & 
                (df_reset['auth_seq_id'] >= start) & 
                (df_reset['auth_seq_id'] <= end))
        df_filtered = df_reset[mask]
        
        if df_filtered.empty:
            raise ValueError(f"No atoms found in range {start}-{end} for chain {chain_id}")
        
        # Create new ID if not provided
        if new_id is None:
            new_id = f"{structure_id}_{chain_id}_{start}_{end}"
        
        # Ensure canonical form
        df_canonical = self._ensure_canonical(df_filtered, new_id)
        
        # Save as new entity
        self._set_frame(new_id, df_canonical)

        if register:
            self.save_entity(new_id, df_canonical, metadata=metadata)

        return df_canonical

    def filter_structure(
        self,
        structure_id: str,
        *,
        filters: Optional[Dict[str, Union[Any, Iterable[Any]]]] = None,
        query: Optional[str] = None,
        new_id: Optional[str] = None,
        register: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> pd.DataFrame:
        """Generic column-based filtering helper."""

        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")

        df_reset = df.reset_index()

        if filters:
            for column, value in filters.items():
                if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
                    df_reset = df_reset[df_reset[column].isin(list(value))]
                else:
                    df_reset = df_reset[df_reset[column] == value]

        if query:
            df_reset = df_reset.query(query)

        if df_reset.empty:
            raise ValueError("Filtering removed all atoms; adjust filters or query")

        target_id = new_id or structure_id
        df_canonical = self._ensure_canonical(df_reset, target_id)
        self._set_frame(target_id, df_canonical)

        if register:
            self.save_entity(target_id, df_canonical, metadata=metadata)

        return df_canonical

    def add_atoms(
        self,
        structure_id: str,
        atoms: List[Dict[str, Any]],
        *,
        assign_new_ids: bool = True,
        new_id: Optional[str] = None,
    ) -> pd.DataFrame:
        """Append atom records to the structure."""

        if not atoms:
            raise ValueError("No atoms provided")

        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")

        df_reset = df.reset_index()

        target_id = new_id or structure_id
        if target_id != structure_id:
            df_reset = df_reset.copy()
            df_reset['structure_id'] = target_id

        next_atom_id = int(df_reset['atom_id'].max()) + 1 if not df_reset.empty else 1
        prepared_rows = []
        for atom in atoms:
            row = dict(atom)
            row['structure_id'] = target_id
            if assign_new_ids or 'atom_id' not in row or pd.isna(row['atom_id']):
                row['atom_id'] = next_atom_id
                next_atom_id += 1
            prepared_rows.append(row)

        new_rows = pd.DataFrame(prepared_rows)
        combined = pd.concat([df_reset, new_rows], ignore_index=True, sort=False)

        df_canonical = self._ensure_canonical(combined, target_id)
        self._set_frame(target_id, df_canonical)

        return df_canonical

    def annotate_structure(
        self,
        structure_id: str,
        annotations: Dict[str, Any],
        *,
        new_id: Optional[str] = None,
        register: bool = False,
    ) -> pd.DataFrame:
        """Apply structured annotations across structure/chain/residue/atom scopes.

        The ``annotations`` payload supports the following optional keys:

        * ``structure`` → {column: value} applied to all atoms.
        * ``chains`` → {chain_id: {column: value}} scoped to chain identifiers.
        * ``residues`` → {(chain_id, auth_seq_id) | "chain:res": {column: value}}.
        * ``atoms`` → {atom_id: {column: value}} targeting specific atoms.

        Unknown targets are ignored, allowing callers to pass best-effort data
        without pre-validating against the structure.
        """

        if not annotations:
            raise ValueError("No annotations provided")

        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")

        target_id = new_id or structure_id
        df_reset = df.reset_index()
        if target_id != structure_id:
            df_reset = df_reset.copy()
            df_reset['structure_id'] = target_id

        def ensure_column(column: str) -> None:
            if column not in df_reset.columns:
                df_reset[column] = pd.NA

        structure_payload = annotations.get('structure') or {}
        for column, value in structure_payload.items():
            ensure_column(column)
            df_reset[column] = value

        chain_payload = annotations.get('chains') or annotations.get('chain') or {}
        for chain_id, column_map in chain_payload.items():
            mask = df_reset['auth_chain_id'] == chain_id
            if not mask.any():
                continue
            for column, value in column_map.items():
                ensure_column(column)
                df_reset.loc[mask, column] = value

        residue_payload = annotations.get('residues') or annotations.get('residue') or {}
        for key, column_map in residue_payload.items():
            if isinstance(key, tuple):
                if len(key) != 2:
                    raise ValueError(f"Residue key tuple must have length 2: {key}")
                chain_id, resid = key
            else:
                if ':' not in str(key):
                    raise ValueError(f"Residue key '{key}' must be tuple or 'chain:resid'")
                chain_id, resid = str(key).split(':', 1)
            try:
                resid_int = int(resid)
            except ValueError as exc:
                raise ValueError(f"Residue identifier must be integer: {resid}") from exc

            mask = (
                (df_reset['auth_chain_id'] == chain_id)
                & (df_reset['auth_seq_id'] == resid_int)
            )
            if not mask.any():
                continue
            for column, value in column_map.items():
                ensure_column(column)
                df_reset.loc[mask, column] = value

        atom_payload = annotations.get('atoms') or annotations.get('atom') or {}
        for atom_id, column_map in atom_payload.items():
            try:
                atom_id_int = int(atom_id)
            except ValueError as exc:
                raise ValueError(f"Atom identifier must be integer: {atom_id}") from exc

            mask = df_reset['atom_id'] == atom_id_int
            if not mask.any():
                continue
            for column, value in column_map.items():
                ensure_column(column)
                df_reset.loc[mask, column] = value

        df_canonical = self._ensure_canonical(df_reset, target_id)
        self._set_frame(target_id, df_canonical)

        if register:
            self.save_entity(target_id, df_canonical)

        return df_canonical

    def add_ligand(
        self,
        structure_id: str,
        ligand_code: str,
        atoms: List[Dict[str, Any]],
        *,
        chain_id: str = 'L',
        residue_id: Optional[int] = None,
        insertion_code: str = '',
    ) -> pd.DataFrame:
        """Convenience wrapper for adding ligand atoms."""

        if residue_id is None:
            residue_id = self._next_residue_id(structure_id, chain_id)

        ligand_atoms: List[Dict[str, Any]] = []
        for atom in atoms:
            row = dict(atom)
            row.setdefault('auth_chain_id', chain_id)
            row.setdefault('label_chain_id', chain_id)
            row.setdefault('auth_seq_id', residue_id)
            row.setdefault('label_seq_id', residue_id)
            row.setdefault('res_name', ligand_code)
            row.setdefault('res_name3l', ligand_code)
            row.setdefault('group', 'HETATM')
            row.setdefault('pdb_ins_code', insertion_code)
            ligand_atoms.append(row)

        return self.add_atoms(structure_id, ligand_atoms)

    def _next_residue_id(self, structure_id: str, chain_id: str) -> int:
        df = self.frames.get(structure_id)
        if df is None:
            df = self.load_entity(structure_id)
            if df is None:
                return 1
        df_reset = df.reset_index()
        mask = df_reset['auth_chain_id'] == chain_id
        return int(df_reset[mask]['auth_seq_id'].max()) + 1 if mask.any() else 1

    def reindex_atom_ids(
        self,
        structure_id: str,
        *,
        start: int = 1,
        new_id: Optional[str] = None,
    ) -> pd.DataFrame:
        """Renumber atom_id sequentially starting from ``start``."""

        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")

        df_reset = df.reset_index()
        df_reset['atom_id'] = range(start, start + len(df_reset))

        target_id = new_id or structure_id
        df_canonical = self._ensure_canonical(df_reset, target_id)
        self._set_frame(target_id, df_canonical)

        return df_canonical

    def annotate_chain(
        self,
        structure_id: str,
        chain_id: str,
        column: str,
        value: Any,
        *,
        default_value: Optional[Any] = None,
        new_id: Optional[str] = None,
    ) -> pd.DataFrame:
        """Add or update a column for all atoms of ``chain_id``."""

        annotations = {
            'chains': {
                chain_id: {column: value},
            }
        }

        if default_value is not None:
            annotations.setdefault('structure', {})[column] = default_value

        return self.annotate_structure(
            structure_id,
            annotations,
            new_id=new_id,
            register=False,
        )
    
    def remove_hetatm(self, structure_id: str) -> pd.DataFrame:
        """
        Remove HETATM records, keeping only ATOM records.
        Modifies structure in place.
        
        Args:
            structure_id: Structure identifier
            
        Returns:
            Modified structure DataFrame
        """
        # Get structure
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Reset index to filter
        df_reset = df.reset_index()
        
        # Keep only ATOM records
        df_filtered = df_reset[df_reset['group'] == 'ATOM']
        
        # Re-canonicalize
        df_canonical = self._ensure_canonical(df_filtered, structure_id)
        
        # Update frames
        self._set_frame(structure_id, df_canonical)
        
        return df_canonical
    
    def apply_transformation(
        self, 
        structure_id: str, 
        rotation: Optional[np.ndarray] = None,
        translation: Optional[np.ndarray] = None
    ) -> pd.DataFrame:
        """
        Apply rotation and/or translation to structure coordinates.
        
        Args:
            structure_id: Structure identifier
            rotation: 3x3 rotation matrix (optional)
            translation: 3D translation vector (optional)
            
        Returns:
            Transformed structure DataFrame
        """
        # Get structure
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Get coordinates
        coords = df[['x', 'y', 'z']].values
        
        # Apply rotation if provided
        if rotation is not None:
            if rotation.shape != (3, 3):
                raise ValueError("Rotation must be a 3x3 matrix")
            coords = coords @ rotation.T
        
        # Apply translation if provided
        if translation is not None:
            if translation.shape != (3,):
                raise ValueError("Translation must be a 3D vector")
            coords += translation
        
        # Update coordinates
        df.loc[:, 'x'] = coords[:, 0]
        df.loc[:, 'y'] = coords[:, 1]
        df.loc[:, 'z'] = coords[:, 2]
        
        # Update frames
        self._set_frame(structure_id, df)
        
        return df
    
    def assign_grns(self, structure_id: str, grn_mapping: Dict[Tuple[str, int], str]) -> pd.DataFrame:
        """
        Assign Generic Residue Numbers to structure residues.
        
        Args:
            structure_id: Structure identifier
            grn_mapping: Dict mapping (chain_id, auth_seq_id) to GRN string
            
        Returns:
            Structure DataFrame with GRN assignments
        """
        # Get structure
        if structure_id not in self.frames:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        else:
            df = self.frames[structure_id].copy()
        
        # Reset index to access residue info
        df_reset = df.reset_index()
        
        # Initialize grn column if not exists
        if 'grn' not in df_reset.columns:
            df_reset['grn'] = ''
        
        # Apply GRN mappings
        for (chain_id, auth_seq_id), grn in grn_mapping.items():
            mask = ((df_reset['auth_chain_id'] == chain_id) & 
                    (df_reset['auth_seq_id'] == auth_seq_id))
            df_reset.loc[mask, 'grn'] = grn
        
        # Re-canonicalize
        df_canonical = self._ensure_canonical(df_reset, structure_id)
        
        # Update frames
        self._set_frame(structure_id, df_canonical)
        
        return df_canonical
    
    def merge_structures(
        self,
        structure_ids: List[str],
        new_id: str,
        chain_mapping: Optional[Dict[str, str]] = None
    ) -> pd.DataFrame:
        """
        Combine multiple structures into single multi-chain structure.
        
        Args:
            structure_ids: List of structure IDs to merge
            new_id: ID for the merged structure
            chain_mapping: Optional mapping of original chain IDs to new chain IDs
                          Format: {f"{structure_id}:{chain_id}": "new_chain_id"}
            
        Returns:
            Merged structure DataFrame
        """
        if not structure_ids:
            raise ValueError("No structures provided to merge")
        
        merged_dfs = []
        used_chains = set()
        
        for struct_id in structure_ids:
            # Load structure
            if struct_id in self.frames:
                df = self.frames[struct_id].copy()
            else:
                df = self.load_entity(struct_id)
                if df is None:
                    raise ValueError(f"Structure {struct_id} not found")
            
            # Reset index to manipulate
            df_reset = df.reset_index()
            
            # Apply chain mapping if provided
            if chain_mapping:
                for (chain_id, auth_seq_id), group in df_reset.groupby(['auth_chain_id', 'auth_seq_id']):
                    mapping_key = f"{struct_id}:{chain_id}"
                    if mapping_key in chain_mapping:
                        new_chain = chain_mapping[mapping_key]
                        df_reset.loc[group.index, 'auth_chain_id'] = new_chain
                        df_reset.loc[group.index, 'label_chain_id'] = new_chain
            else:
                # Auto-rename chains to avoid conflicts
                chains = df_reset['auth_chain_id'].unique()
                for chain in chains:
                    if chain in used_chains:
                        # Find new chain name
                        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                            if c not in used_chains:
                                mask = df_reset['auth_chain_id'] == chain
                                df_reset.loc[mask, 'auth_chain_id'] = c
                                df_reset.loc[mask, 'label_chain_id'] = c
                                used_chains.add(c)
                                break
                    else:
                        used_chains.add(chain)
            
            merged_dfs.append(df_reset)
        
        # Concatenate all structures
        merged_df = pd.concat(merged_dfs, ignore_index=True)
        
        # Ensure canonical form with new structure ID
        df_canonical = self._ensure_canonical(merged_df, new_id)
        
        # Save merged structure
        self._set_frame(new_id, df_canonical)

        return df_canonical

    # ---------- Sequence Extraction & Registration ----------

    def collect_chain_sequences(
        self,
        structure_ids: Union[str, Iterable[str]],
        *,
        chain_filter: Optional[Union[Iterable[str], Dict[str, Iterable[str]]]] = None,
        one_letter: bool = True,
        min_length: int = 1,
    ) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Collect per-chain sequences for the provided structures without registering them."""

        from protos.analysis.structure.sequence import extract_all_sequences

        if isinstance(structure_ids, str):
            structure_list = [structure_ids]
        else:
            structure_list = list(structure_ids)

        filter_map: Dict[str, Optional[Set[str]]] = {}
        if chain_filter is None:
            filter_map = {sid: None for sid in structure_list}
        elif isinstance(chain_filter, dict):
            for sid, chains in chain_filter.items():
                if isinstance(chains, str):
                    filter_map[sid] = {chains}
                elif chains is None:
                    filter_map[sid] = None
                else:
                    filter_map[sid] = set(chains)
        else:
            if isinstance(chain_filter, str):
                chain_set = {chain_filter}
            else:
                chain_set = set(chain_filter)
            filter_map = {sid: chain_set for sid in structure_list}

        collected: Dict[str, Dict[str, Dict[str, Any]]] = {}

        for structure_id in structure_list:
            df = self.frames.get(structure_id)
            if df is None:
                df = self.load_entity(structure_id)
            if df is None:
                self.logger.warning(f"Structure {structure_id} not found while collecting chain sequences")
                continue

            sequences = extract_all_sequences(df, one_letter=one_letter, min_length=min_length)

            allowed_chains = filter_map.get(structure_id)
            if allowed_chains is not None:
                sequences = {cid: seq for cid, seq in sequences.items() if cid in allowed_chains}

            df_reset = df.reset_index()
            chain_payloads: Dict[str, Dict[str, Any]] = {}

            for chain_id, sequence in sequences.items():
                chain_df = df_reset[df_reset['auth_chain_id'] == chain_id]
                residue_series = chain_df['auth_seq_id'].dropna()
                start_res = int(residue_series.min()) if not residue_series.empty else None
                end_res = int(residue_series.max()) if not residue_series.empty else None

                entity_name = self._make_chain_entity_name(structure_id, chain_id)
                residue_span = [start_res, end_res] if start_res is not None and end_res is not None else None

                chain_payloads[chain_id] = {
                    'entity_name': entity_name,
                    'sequence': sequence,
                    'length': len(sequence),
                    'residue_span': residue_span,
                    'metadata': {
                        'structure_id': structure_id,
                        'chain_id': chain_id,
                        'sequence_length': len(sequence),
                        'residue_span': residue_span,
                    },
                }

            collected[structure_id] = chain_payloads

        return collected

    def register_chain_sequences(
        self,
        structure_ids: Union[str, Iterable[str]],
        *,
        chain_filter: Optional[Union[Iterable[str], Dict[str, Iterable[str]]]] = None,
        one_letter: bool = True,
        min_length: int = 1,
        semantic_role: str = 'chain_sequence',
        dataset_prefix: Optional[str] = None,
        create_dataset: bool = True,
        overwrite: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        """Register chain sequences for the selected structures and emit relationships."""

        collected = self.collect_chain_sequences(
            structure_ids,
            chain_filter=chain_filter,
            one_letter=one_letter,
            min_length=min_length,
        )

        structure_list = [structure_ids] if isinstance(structure_ids, str) else list(structure_ids)

        loader = self._get_sequence_loader()
        summary: Dict[str, Dict[str, Any]] = {}

        for structure_id in structure_list:
            chains = collected.get(structure_id, {})
            if not chains:
                summary[structure_id] = {
                    'chains': {},
                    'registered_entities': [],
                    'dataset': None,
                }
                continue

            records = []
            for chain_id, payload in chains.items():
                metadata = dict(payload['metadata'])
                metadata['semantic_role'] = semantic_role
                metadata['source_processor'] = self.processor_type

                records.append({
                    'name': payload['entity_name'],
                    'sequence': payload['sequence'],
                    'metadata': metadata,
                })

            dataset_name = None
            if create_dataset:
                if dataset_prefix:
                    dataset_name = f"{dataset_prefix}_{structure_id}"
                else:
                    dataset_name = f"{structure_id}_chains"

            dataset_metadata = {
                'structure_id': structure_id,
                'semantic_role': semantic_role,
                'chain_count': len(records),
            }

            registration = loader.register_sequence_records(
                records,
                dataset_name=dataset_name,
                dataset_metadata=dataset_metadata,
                overwrite=overwrite,
            )

            for chain_id, payload in chains.items():
                seq_name = payload['entity_name']
                rel_metadata = {
                    'chain_id': chain_id,
                    'structure_id': structure_id,
                    'sequence_length': payload['length'],
                    'residue_span': payload['residue_span'],
                    'semantic_role': semantic_role,
                }
                self.entity_registry.add_relationship(
                    seq_name,
                    structure_id,
                    'derived_from',
                    metadata=rel_metadata,
                )

            summary[structure_id] = {
                'chains': chains,
                'registered_entities': registration.get('entities', []),
                'dataset': registration.get('dataset'),
            }

        return summary

    def _make_chain_entity_name(self, structure_id: str, chain_id: str) -> str:
        sanitized_chain = chain_id.replace(' ', '_')
        return f"{structure_id}_chain_{sanitized_chain}"

    def compute_water_networks(
        self,
        structure_ids: Optional[Union[str, Iterable[str]]] = None,
        *,
        residue_cutoff: float = 3.4,
        water_water_cutoff: float = 3.4,
        hydrogen_bond_cutoff: float = 3.2,
        property_table_name: Optional[str] = None,
        property_metadata: Optional[Dict[str, Any]] = None,
        allow_create_property_table: bool = False,
    ) -> Dict[str, Any]:
        """Analyze water-mediated networks for the provided structures."""

        if structure_ids is None:
            structure_list = self.list_entities()
        elif isinstance(structure_ids, str):
            structure_list = [structure_ids]
        else:
            structure_list = list(structure_ids)

        structure_list = list(dict.fromkeys(structure_list))

        analysis_results: Dict[str, Any] = {}
        errors: Dict[str, str] = {}
        property_rows: List[Dict[str, Any]] = []

        for structure_id in structure_list:
            frame = self.load_entity(structure_id)
            if frame is None:
                errors[structure_id] = "structure not found"
                continue

            try:
                analysis = analyze_water_networks(
                    frame,
                    residue_cutoff=residue_cutoff,
                    water_water_cutoff=water_water_cutoff,
                    hydrogen_bond_cutoff=hydrogen_bond_cutoff,
                )
            except Exception as exc:  # pragma: no cover - defensive logging
                self.logger.exception(
                    "Failed to compute water networks for %s", structure_id
                )
                errors[structure_id] = str(exc)
                continue

            analysis_results[structure_id] = analysis

            summary = analysis.get('summary') or {}
            if property_table_name and summary:
                property_rows.append(
                    {
                        'scope': [{'format': 'structure', 'name': structure_id}],
                        'entity_name': structure_id,
                        'network_count': summary.get('network_count'),
                        'residue_count': summary.get('residue_count'),
                        'water_count': summary.get('water_count'),
                        'bridging_water_count': summary.get('bridging_water_count'),
                        'max_residue_path_length': summary.get('max_residue_path_length'),
                        'residues_with_grn': summary.get('residues_with_grn'),
                        'parameters': {
                            'residue_cutoff': residue_cutoff,
                            'water_water_cutoff': water_water_cutoff,
                            'hydrogen_bond_cutoff': hydrogen_bond_cutoff,
                        },
                    }
                )

        property_table_recorded: Optional[str] = None
        if property_rows and property_table_name:
            from protos.processing.property import PropertyProcessor

            prop_proc = PropertyProcessor()
            metadata = property_metadata.copy() if property_metadata else {}
            metadata.update(
                {
                    'residue_cutoff': residue_cutoff,
                    'water_water_cutoff': water_water_cutoff,
                    'hydrogen_bond_cutoff': hydrogen_bond_cutoff,
                }
            )

            prop_proc.record_properties(
                property_table_name,
                property_rows,
                metadata=metadata,
                allow_create=allow_create_property_table,
                materialize_entries=False,
            )
            property_table_recorded = property_table_name

        return {
            'structures': analysis_results,
            'errors': errors,
            'property_table': property_table_recorded,
        }

    def list_related_sequences(
        self,
        structure_ids: Optional[Union[str, Iterable[str]]] = None,
        *,
        rel_type: str = 'derived_from',
        direction: str = 'incoming',
        include_unloaded: bool = False,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Expose related sequence entities for the given structures."""

        if structure_ids is None:
            if include_unloaded:
                structure_ids = self.list_entities()
            else:
                structure_ids = self.structure_ids

        return self.resolve_related_entities(
            structure_ids,
            rel_type=rel_type,
            direction=direction,
            format_type='sequence',
        )

    def list_dataset_related_sequences(
        self,
        dataset_name: str,
        *,
        rel_type: str = 'derived_from',
        direction: str = 'incoming',
    ) -> Dict[str, List[Dict[str, Any]]]:
        """List related sequence entities for every structure in a dataset."""

        return self.resolve_related_entities_for_dataset(
            dataset_name,
            rel_type=rel_type,
            direction=direction,
            format_type='sequence',
        )

    def align_structures(
        self,
        structure_ids: List[str],
        reference_id: str,
        method: str = 'cealign',
        atom_selection: str = 'CA',
        apply_transform: bool = True,
        chain_id: Optional[str] = None,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
    ) -> Tuple[Dict[str, Dict[str, float]], Dict[str, 'AlignmentResult']]:
        """Align ``structure_ids`` against ``reference_id``.

        Returns a nested dictionary mapping target→reference→RMSD along with
        the detailed :class:`AlignmentResult` objects.
        """

        mobile_ids = [sid for sid in structure_ids if sid != reference_id]
        rmsd_map, results = self.align_all_to_reference(
            reference_id,
            mobile_ids,
            method=method,
            atom_selection=atom_selection,
            apply_transform=apply_transform,
            chain_id=chain_id,
            cealign_window=cealign_window,
            cealign_max_gap=cealign_max_gap,
        )
        return rmsd_map, results


    def align_and_record(
        self,
        structure_ids: Iterable[str],
        reference_id: str,
        *,
        method: str = 'cealign',
        atom_selection: str = 'CA',
        chain_id: Optional[str] = None,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
        apply_transform: bool = True,
        save_aligned: bool = False,
        summary_name: Optional[str] = None,
        summary_metadata: Optional[Dict[str, Any]] = None,
        aligned_dataset_name: Optional[str] = None,
        aligned_dataset_include_reference: bool = False,
        property_table_name: Optional[str] = None,
    ) -> Tuple[Dict[str, Any], Dict[str, 'AlignmentResult']]:
        """Align structures, persist optional artifacts, and return a summary."""

        structure_list = list(dict.fromkeys(structure_ids))
        structure_list = [sid for sid in structure_list if sid != reference_id]

        rmsd_map, results = self.align_structures(
            structure_list,
            reference_id,
            method=method,
            atom_selection=atom_selection,
            apply_transform=apply_transform,
            chain_id=chain_id,
            cealign_window=cealign_window,
            cealign_max_gap=cealign_max_gap,
        )

        pairwise_rows: List[Dict[str, Any]] = []
        rmsd_values: List[float] = []
        errors: Dict[str, str] = {}
        aligned_entities: List[str] = []

        for target_id, ref_map in rmsd_map.items():
            for ref_id, rmsd in ref_map.items():
                if rmsd is not None and not (isinstance(rmsd, float) and np.isnan(rmsd)):
                    rmsd_values.append(float(rmsd))
                result = results.get(target_id)
                error = result.error if result else None
                if error:
                    errors[target_id] = error
                pairwise_rows.append(
                    {
                        'target_id': target_id,
                        'reference_id': ref_id,
                        'rmsd': None if rmsd is None or (isinstance(rmsd, float) and np.isnan(rmsd)) else float(rmsd),
                        'algorithm': result.algorithm if result else method,
                        'error': error,
                    }
                )

        timestamp = datetime.utcnow().isoformat()
        global_stats: Dict[str, Any]
        if rmsd_values:
            global_stats = {
                'count': len(rmsd_values),
                'min': float(np.min(rmsd_values)),
                'max': float(np.max(rmsd_values)),
                'mean': float(np.mean(rmsd_values)),
                'median': float(np.median(rmsd_values)),
            }
        else:
            global_stats = {
                'count': 0,
                'min': None,
                'max': None,
                'mean': None,
                'median': None,
            }

        aligned_result_map = {
            sid: {
                'aligned_id': res.aligned_id,
                'rmsd': None if res.rmsd is None or np.isnan(res.rmsd) else float(res.rmsd),
                'success': res.error is None,
                'algorithm': res.algorithm,
                'error': res.error,
            }
            for sid, res in results.items()
            if res is not None
        }

        summary_payload: Dict[str, Any] = {
            'reference_id': reference_id,
            'structure_ids': [reference_id] + structure_list,
            'method': method,
            'apply_transform': apply_transform,
            'atom_selection': atom_selection,
            'chain_id': chain_id,
            'parameters': {
                'cealign_window': cealign_window,
                'cealign_max_gap': cealign_max_gap,
            },
            'timestamp': timestamp,
            'rmsd': {
                'global': global_stats,
                'pairwise': pairwise_rows,
            },
            'results': aligned_result_map,
            'errors': errors,
        }

        if summary_metadata:
            extra = {k: v for k, v in summary_metadata.items() if v is not None}
            summary_payload.setdefault('metadata', {}).update(extra)

        if save_aligned:
            for struct_id, res in results.items():
                if res is None or res.error:
                    continue
                if struct_id == reference_id and res.aligned_id == reference_id:
                    continue
                frame = self.frames.get(res.aligned_id)
                if frame is None:
                    frame = self.load_entity(res.aligned_id)
                if frame is None:
                    continue
                aligned_meta = {
                    'aligned_to': reference_id,
                    'alignment_method': res.algorithm,
                    'rmsd': None if res.rmsd is None or np.isnan(res.rmsd) else float(res.rmsd),
                    'aligned_at': timestamp,
                }
                if summary_metadata:
                    aligned_meta.update({k: v for k, v in summary_metadata.items() if v is not None})
                self.save_entity(res.aligned_id, frame, metadata=aligned_meta)
                aligned_entities.append(res.aligned_id)

            aligned_entities = list(dict.fromkeys(aligned_entities))

        aligned_entities = list(dict.fromkeys(aligned_entities))
        summary_payload['aligned_entities'] = aligned_entities
        summary_payload['aligned_entities'] = aligned_entities

        summary_basename = summary_name or f"{reference_id}_alignment"
        safe_summary_name = self._sanitize_filename(summary_basename)
        alignment_dir = Path(self.paths.get_subdir_path('structure', 'alignments_dir'))
        summary_path = alignment_dir / f"{safe_summary_name}.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary_payload, indent=2), encoding='utf-8')
        summary_rel_path = summary_path.relative_to(self.paths.data_root)
        summary_payload['summary_file'] = str(summary_rel_path)
        summary_payload['summary_name'] = summary_basename

        dataset_entities = [reference_id] + structure_list
        dataset_entities = list(dict.fromkeys(dataset_entities))
        dataset_metadata = {
            'summary_file': str(summary_rel_path),
            'aligned_entities': aligned_entities,
            'timestamp': timestamp,
            'method': method,
            'apply_transform': apply_transform,
            'atom_selection': atom_selection,
            'chain_id': chain_id,
            'global_stats': global_stats,
        }
        if summary_metadata:
            dataset_metadata.update({k: v for k, v in summary_metadata.items() if v is not None})

        summary_dataset_name = safe_summary_name

        if property_table_name:
            from protos.processing.property import PropertyProcessor

            prop_proc = PropertyProcessor()
            property_rows = []
            for row in pairwise_rows:
                target = row['target_id']
                reference = row['reference_id']
                if target == reference:
                    continue
                scope = [{'format': 'structure', 'name': target}]
                if reference:
                    scope.append({'format': 'structure', 'name': reference})
                status = 'ok' if not row['error'] else 'error'
                property_row = {
                    'scope': scope,
                    'entity_name': target,
                    'reference': reference,
                    'rmsd': row['rmsd'],
                    'algorithm': row['algorithm'],
                    'status': status,
                    'summary_dataset': summary_dataset_name,
                }
                if row['error']:
                    property_row['error_message'] = row['error']
                property_rows.append(property_row)

            if property_rows:
                prop_metadata = {
                    'reference_id': reference_id,
                    'method': method,
                    'atom_selection': atom_selection,
                    'chain_id': chain_id,
                    'apply_transform': apply_transform,
                    'summary_dataset': summary_dataset_name,
                    'timestamp': timestamp,
                }
                if summary_metadata:
                    prop_metadata.update({k: v for k, v in summary_metadata.items() if v is not None})

                prop_proc.record_properties(
                    property_table_name,
                    property_rows,
                    metadata=prop_metadata,
                    allow_create=True,
                    materialize_entries=False,
                )
                summary_payload['property_table'] = property_table_name
                dataset_metadata['property_table'] = property_table_name

        summary_dataset = self.create_dataset(summary_dataset_name, dataset_entities, dataset_metadata)
        summary_payload['summary_dataset'] = summary_dataset

        if save_aligned and aligned_dataset_name:
            safe_aligned_name = self._sanitize_filename(aligned_dataset_name)
            aligned_meta = {
                'reference_id': reference_id,
                'source_dataset': summary_dataset,
                'method': method,
                'atom_selection': atom_selection,
                'chain_id': chain_id,
                'apply_transform': apply_transform,
                'timestamp': timestamp,
            }
            if summary_metadata:
                aligned_meta.update({k: v for k, v in summary_metadata.items() if v is not None})

            aligned_dataset_entities = list(dict.fromkeys(aligned_entities))
            if aligned_dataset_include_reference and reference_id not in aligned_dataset_entities:
                aligned_dataset_entities.insert(0, reference_id)

            if aligned_dataset_entities:
                self.create_dataset(safe_aligned_name, aligned_dataset_entities, aligned_meta)
                summary_payload['aligned_dataset'] = safe_aligned_name
            else:
                summary_payload['aligned_dataset'] = None
        elif save_aligned and aligned_entities:
            summary_payload['aligned_dataset'] = None

        return summary_payload, results

    def export_aligned_structures(
        self,
        structure_ids: Optional[Iterable[str]] = None,
        *,
        output_dir: Optional[Union[str, Path]] = None,
        overwrite: bool = False,
        dataset_name: Optional[str] = None,
        export_format: str = 'cif',
    ) -> Dict[str, Path]:
        """Export aligned structures to CIF files using the exporter."""

        if dataset_name:
            structure_ids = self.dataset_manager.get_dataset_entities(dataset_name)
        if structure_ids is None:
            raise ValueError('structure_ids or dataset_name must be provided')

        ids = list(dict.fromkeys(structure_ids))
        export_root = Path(output_dir) if output_dir else Path(self.paths.get_subdir_path('structure', 'alignments_dir')) / 'exports'
        export_root.mkdir(parents=True, exist_ok=True)

        exported: Dict[str, Path] = {}
        fmt = export_format.lstrip('.')
        suffix = export_format if export_format.startswith('.') else f'.{fmt}'

        for struct_id in ids:
            out_path = export_root / f"{self._sanitize_filename(struct_id)}{suffix}"
            exported[struct_id] = self.export_entity(struct_id, out_path, format=fmt, overwrite=overwrite)
        return exported

    def align_pair(
        self,
        reference_id: str,
        mobile_id: str,
        *,
        method: str = 'cealign',
        atom_selection: str = 'CA',
        apply_transform: bool = True,
        chain_id: Optional[str] = None,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
    ) -> Optional['AlignmentResult']:
        """Align a single mobile structure against a reference."""

        results = self.alignment_engine.align(
            [mobile_id],
            reference_id,
            method=method,
            atom_selection=atom_selection,
            chain_id=chain_id,
            apply_transform=apply_transform,
            cealign_window=cealign_window,
            cealign_max_gap=cealign_max_gap,
        )
        return results.get(mobile_id)

    def align_all_to_reference(
        self,
        reference_id: str,
        structure_ids: List[str],
        *,
        method: str = 'cealign',
        atom_selection: str = 'CA',
        apply_transform: bool = True,
        chain_id: Optional[str] = None,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
    ) -> Tuple[Dict[str, Dict[str, float]], Dict[str, 'AlignmentResult']]:
        """Align every structure in ``structure_ids`` to the reference."""

        targets = [sid for sid in structure_ids if sid != reference_id]
        raw_results = self.alignment_engine.align(
            targets,
            reference_id,
            method=method,
            atom_selection=atom_selection,
            chain_id=chain_id,
            apply_transform=apply_transform,
            cealign_window=cealign_window,
            cealign_max_gap=cealign_max_gap,
        )

        rmsd_map: Dict[str, Dict[str, float]] = {}
        for struct_id, result in raw_results.items():
            rmsd = result.rmsd if result.error is None else float('nan')
            rmsd_map.setdefault(struct_id, {})[reference_id] = rmsd

        return rmsd_map, raw_results

    def align_one_vs_all(
        self,
        reference_id: str,
        candidates: List[str],
        **kwargs: Any,
    ) -> Tuple[Optional[str], Optional['AlignmentResult'], Dict[str, 'AlignmentResult']]:
        """Align reference against candidates and return the best match."""

        rmsd_map, results = self.align_all_to_reference(
            reference_id,
            candidates,
            **kwargs,
        )

        best_id: Optional[str] = None
        best_result: Optional['AlignmentResult'] = None
        for struct_id, result in results.items():
            if result.error:
                continue
            if best_result is None or result.rmsd < best_result.rmsd:
                best_id = struct_id
                best_result = result

        return best_id, best_result, results

    def align_all_vs_reference(
        self,
        reference_id: str,
        structure_ids: List[str],
        **kwargs: Any,
    ) -> Tuple[Dict[str, Dict[str, float]], Dict[str, 'AlignmentResult']]:
        """Alias for :meth:`align_all_to_reference` for readability."""

        return self.align_all_to_reference(reference_id, structure_ids, **kwargs)

    def align_all_vs_all(
        self,
        structure_ids: List[str],
        *,
        method: str = 'cealign',
        atom_selection: str = 'CA',
        apply_transform: bool = False,
        chain_id: Optional[str] = None,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
    ) -> Tuple[Dict[str, Dict[str, float]], Dict[Tuple[str, str], 'AlignmentResult']]:
        """Perform all-vs-all alignments and return an RMSD mapping."""

        order = list(structure_ids)
        rmsd_map: Dict[str, Dict[str, float]] = {sid: {} for sid in order}
        pair_results: Dict[Tuple[str, str], 'AlignmentResult'] = {}

        for i, reference_id in enumerate(order):
            rmsd_map[reference_id][reference_id] = 0.0
            for j in range(i + 1, len(order)):
                mobile_id = order[j]
                result = self.align_pair(
                    reference_id,
                    mobile_id,
                    method=method,
                    atom_selection=atom_selection,
                    apply_transform=apply_transform,
                    chain_id=chain_id,
                    cealign_window=cealign_window,
                    cealign_max_gap=cealign_max_gap,
                )
                pair_results[(mobile_id, reference_id)] = result
                pair_results[(reference_id, mobile_id)] = result
                rmsd = float('nan')
                if result is not None and result.error is None:
                    rmsd = result.rmsd
                rmsd_map[mobile_id][reference_id] = rmsd
                rmsd_map[reference_id][mobile_id] = rmsd

        return rmsd_map, pair_results

    def annotate_structures_with_grn(
        self,
        structure_ids: Iterable[str],
        *,
        reference_table: str,
        protein_family: str,
        output_table: Optional[str] = None,
        dataset_prefix: str = "grn_chain_sequences",
        chain_filter: Optional[Union[str, Iterable[str], Dict[str, Iterable[str]]]] = None,
        allow_create: bool = False,
        materialize_entries: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
        return_summary: bool = False,
    ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, Any]]]:
        """Annotate structure chains with GRNs using registered references."""

        structure_list = list(structure_ids)
        if not structure_list:
            raise ValueError("No structure IDs provided for GRN annotation")

        self.register_chain_sequences(
            structure_list,
            chain_filter=chain_filter,
            dataset_prefix=dataset_prefix,
            create_dataset=True,
            overwrite=False,
        )

        from protos.processing.sequence import SequenceProcessor
        from protos.processing.grn import GRNProcessor

        seq_proc = SequenceProcessor()
        grn_proc = GRNProcessor()

        sequence_map: Dict[str, str] = {}
        sequence_structure_map: Dict[str, List[Dict[str, Any]]] = {}

        related = self.list_related_sequences(structure_list, include_unloaded=True)
        for structure_id, relations in related.items():
            for payload in relations:
                seq_name = payload.get("name")
                if not seq_name:
                    continue
                chain_id = payload.get("metadata", {}).get("chain_id")
                sequence_data = seq_proc.load_entity(seq_name)
                if isinstance(sequence_data, str) and sequence_data:
                    sequence_map[seq_name] = sequence_data
                    sequence_structure_map.setdefault(seq_name, []).append(
                        {"name": structure_id, "chain_id": chain_id}
                    )

        if not sequence_map:
            raise RuntimeError("No chain sequences available for GRN annotation")

        annotations, summary = grn_proc.annotate_sequences(
            sequence_map,
            reference_table=reference_table,
            protein_family=protein_family,
        )

        if output_table:
            table_metadata: Dict[str, Any] = {
                "reference_table": reference_table,
                "protein_family": protein_family,
                "structure_ids": structure_list,
                "sequence_count": len(annotations),
                "materialize_entries": materialize_entries,
                **summary.get("global", {}),
            }
            if metadata:
                table_metadata.update(metadata)

            per_entity_metadata: Dict[str, Dict[str, Any]] = {}
            for seq_name, info in summary.get("per_sequence", {}).items():
                entry_metadata = dict(info)
                if seq_name in sequence_structure_map:
                    entry_metadata["related_structures"] = sequence_structure_map[seq_name]
                per_entity_metadata[seq_name] = entry_metadata

            grn_proc.record_table(
                output_table,
                annotations,
                metadata=table_metadata,
                per_entity_metadata=per_entity_metadata,
                allow_create=allow_create,
                link_entities=True,
            )

        if return_summary:
            return annotations, summary
        return annotations
    
    def orient_structure(
        self,
        structure_id: str,
        method: str = 'principal_axes',
        axis_order: Optional[List[int]] = None
    ) -> pd.DataFrame:
        """
        Orient structure using standard methods.
        
        Args:
            structure_id: Structure identifier
            method: Orientation method:
                   - 'principal_axes': Align to principal axes of inertia
                   - 'membrane_normal': Orient with membrane normal as Z axis
            axis_order: For principal_axes, order of axes mapping [0,1,2]
            
        Returns:
            Oriented structure DataFrame
        """
        # Get structure
        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        
        if method == 'principal_axes':
            # Calculate center of mass
            from protos.analysis.structure.geometry import calculate_center_of_mass
            com = calculate_center_of_mass(df)
            
            # Center structure
            coords = df[['x', 'y', 'z']].values
            centered_coords = coords - com
            
            # Calculate inertia tensor
            I = np.zeros((3, 3))
            for coord in centered_coords:
                I[0, 0] += coord[1]**2 + coord[2]**2
                I[1, 1] += coord[0]**2 + coord[2]**2
                I[2, 2] += coord[0]**2 + coord[1]**2
                I[0, 1] -= coord[0] * coord[1]
                I[0, 2] -= coord[0] * coord[2]
                I[1, 2] -= coord[1] * coord[2]
            
            I[1, 0] = I[0, 1]
            I[2, 0] = I[0, 2]
            I[2, 1] = I[1, 2]
            
            # Get principal axes (eigenvectors)
            eigenvalues, eigenvectors = np.linalg.eig(I)
            
            # Sort by eigenvalue (largest to smallest)
            idx = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            
            # Apply axis ordering if specified
            if axis_order:
                eigenvectors = eigenvectors[:, axis_order]
            
            # Create rotation matrix to align principal axes with coordinate axes
            rotation = eigenvectors.T
            
            # Apply transformation (translate to origin, rotate, translate back)
            df_oriented = self.apply_transformation(
                structure_id, rotation=rotation, translation=-com
            )
            
            # Translate back to original COM position
            df_oriented = self.apply_transformation(
                structure_id, translation=com
            )
            
        elif method == 'membrane_normal':
            from protos.analysis.structure.membrane import calculate_membrane_normal
            
            # Calculate membrane normal
            normal = calculate_membrane_normal(df)
            
            # Create rotation to align normal with Z axis
            z_axis = np.array([0, 0, 1])
            from protos.analysis.structure.geometry import calculate_rotation_matrix
            rotation = calculate_rotation_matrix(normal, z_axis)
            
            # Apply rotation
            df_oriented = self.apply_transformation(structure_id, rotation=rotation)
            
        else:
            raise ValueError(f"Unknown orientation method: {method}")
        
        # Update frames
        self._set_frame(structure_id, df_oriented)
        
        return df_oriented
    
    def renumber_residues(
        self,
        structure_id: str,
        start: int = 1,
        by_chain: bool = True,
        keep_mapping: bool = True
    ) -> pd.DataFrame:
        """
        Renumber residues sequentially.
        
        Args:
            structure_id: Structure identifier
            start: Starting residue number
            by_chain: If True, restart numbering for each chain
            keep_mapping: If True, stores old->new mapping in metadata
            
        Returns:
            Structure with renumbered residues
        """
        # Get structure
        if structure_id in self.frames:
            df = self.frames[structure_id]
        else:
            df = self.load_entity(structure_id)
            if df is None:
                raise ValueError(f"Structure {structure_id} not found")
        
        # Reset index to access all data
        df_reset = df.reset_index()
        
        # Create mapping
        residue_mapping = {}
        
        if by_chain:
            # Renumber each chain separately
            for chain_id in df_reset['auth_chain_id'].unique():
                chain_mask = df_reset['auth_chain_id'] == chain_id
                chain_residues = df_reset[chain_mask]['auth_seq_id'].unique()
                chain_residues = sorted(chain_residues)
                
                # Create mapping for this chain
                for i, old_resid in enumerate(chain_residues):
                    new_resid = start + i
                    residue_mapping[(chain_id, old_resid)] = new_resid
                    
                    # Apply to dataframe
                    mask = chain_mask & (df_reset['auth_seq_id'] == old_resid)
                    df_reset.loc[mask, 'auth_seq_id'] = new_resid
                    df_reset.loc[mask, 'label_seq_id'] = new_resid
        else:
            # Global renumbering
            all_residues = df_reset.groupby(['auth_chain_id', 'auth_seq_id']).size()
            
            for i, (chain_id, old_resid) in enumerate(all_residues.index):
                new_resid = start + i
                residue_mapping[(chain_id, old_resid)] = new_resid
                
                # Apply to dataframe
                mask = ((df_reset['auth_chain_id'] == chain_id) & 
                        (df_reset['auth_seq_id'] == old_resid))
                df_reset.loc[mask, 'auth_seq_id'] = new_resid
                df_reset.loc[mask, 'label_seq_id'] = new_resid
        
        # Re-canonicalize
        df_canonical = self._ensure_canonical(df_reset, structure_id)
        
        # Update frames
        self._set_frame(structure_id, df_canonical)
        
        # Store mapping in metadata if requested
        if keep_mapping:
            self.logger.info(f"Residue renumbering mapping: {residue_mapping}")
        
        return df_canonical

    # ---------- Legacy Compatibility Helpers ----------

    def load_structures(
        self,
        structure_ids: Optional[Union[str, Iterable[str]]] = None,
        **_,
    ) -> List[str]:
        """Compatibility wrapper for legacy load_structures API."""

        warnings.warn(
            "StructureProcessor.load_structures is deprecated; use load_entity instead",
            DeprecationWarning,
            stacklevel=2,
        )

        if structure_ids is None:
            structure_list = list(self.list_entities())
        elif isinstance(structure_ids, str):
            structure_list = [structure_ids]
        else:
            structure_list = list(structure_ids)

        loaded: List[str] = []
        for structure_id in structure_list:
            try:
                df = self.load_entity(structure_id)
                if df is not None:
                    loaded.append(structure_id)
            except Exception as exc:  # noqa: BLE001
                self.logger.warning(
                    "Could not load structure %s via compatibility wrapper: %s",
                    structure_id,
                    exc,
                )
        return loaded

    def get_chains(self, structure_id: str) -> List[str]:
        """Return chain identifiers for a structure (legacy helper)."""

        df = self.frames.get(structure_id)
        if df is None:
            df = self.load_entity(structure_id)
        if df is None:
            return []

        df_reset = df.reset_index()
        chains_series = (
            df_reset['auth_chain_id']
            if 'auth_chain_id' in df_reset
            else df_reset.get('chain_id')
        )
        if chains_series is None:
            return []

        chain_values = sorted(
            {
                str(c).strip()
                for c in chains_series.dropna().unique()
                if str(c).strip()
            }
        )
        return chain_values

    def get_sequence(self, structure_id: str, chain_id: str = 'A') -> Optional[str]:
        """Return a single chain sequence using the unified extractor."""

        collected = self.collect_chain_sequences(
            [structure_id],
            chain_filter={structure_id: [chain_id]},
        )
        chains = collected.get(structure_id, {})
        payload = chains.get(chain_id)
        if not payload:
            # Attempt relaxed match for legacy callers with stripped chain IDs
            for key, value in chains.items():
                if key.strip() == chain_id:
                    payload = value
                    break
        if not payload:
            return None
        return payload.get('sequence')

    def get_seq(self, structure_id: str, chain_id: str = 'A') -> Optional[str]:
        """Alias maintained for backward compatibility."""

        return self.get_sequence(structure_id, chain_id)

    def get_all_sequences(
        self,
        structure_ids: Optional[Union[str, Iterable[str]]] = None,
    ) -> Dict[str, str]:
        """Aggregate chain sequences for one or more structures."""

        if structure_ids is None:
            if self.frames:
                structure_list = list(self.frames.keys())
            else:
                structure_list = list(self.list_entities())
        elif isinstance(structure_ids, str):
            structure_list = [structure_ids]
        else:
            structure_list = list(structure_ids)

        collected = self.collect_chain_sequences(structure_list)

        sequences: Dict[str, str] = {}
        for struct_id, chains in collected.items():
            for chain_id, payload in chains.items():
                seq_name = payload.get('entity_name') or f"{struct_id}_chain_{chain_id}"
                sequence = payload.get('sequence')
                if sequence:
                    sequences[seq_name] = sequence

        # Cache for compatibility with legacy workflows
        self.chain_dict = sequences.copy()
        return sequences

    def get_seq_dict(
        self,
        load_file: bool = False,
        *_args,
        **_kwargs,
    ) -> Dict[str, str]:
        """Compatibility wrapper mimicking legacy get_seq_dict."""

        if load_file:
            raise NotImplementedError(
                "Loading legacy chain FASTA files is no longer supported"
            )
        return self.get_all_sequences()

    def get_chain_dict(self) -> Dict[str, str]:
        """Alias for legacy get_chain_dict semantics."""

        return self.get_all_sequences()

    def get_ca_coordinates(self, structure_id: str, chain_id: str) -> np.ndarray:
        """Return CA coordinates for a specific chain (legacy helper)."""

        df = self.frames.get(structure_id)
        if df is None:
            df = self.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure {structure_id} not found")

        df_reset = df.reset_index()
        if 'auth_chain_id' not in df_reset and 'chain_id' not in df_reset:
            raise ValueError("Structure data missing 'auth_chain_id' column")

        chain_series = (
            df_reset['auth_chain_id']
            if 'auth_chain_id' in df_reset
            else df_reset['chain_id']
        ).astype(str).str.strip()
        mask = chain_series == chain_id

        atom_col = 'atom_name' if 'atom_name' in df_reset else 'label_atom_id'
        atom_series = df_reset[atom_col].astype(str).str.upper()
        chain_df = df_reset[mask & (atom_series == 'CA')]

        if chain_df.empty:
            raise ValueError(f"No CA atoms found for {structure_id} chain {chain_id}")

        return chain_df[['x', 'y', 'z']].to_numpy(dtype=float)

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
        # Always save as PKL (no other formats supported)
        self.save_entity(name, structure_df, metadata=kwargs.get('metadata'))

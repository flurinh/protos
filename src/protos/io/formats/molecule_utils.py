"""Utilities for converting small molecules to canonical structure DataFrames."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:  # Optional RDKit dependency
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:  # pragma: no cover - optional
    Chem = None

from .structure_schema import STRUCT_COLUMN_DTYPE, SORTED_STRUCT_COLUMNS


def require_rdkit() -> None:
    if Chem is None:
        raise ImportError(
            "RDKit is required for ligand handling. Install with 'pip install rdkit-pypi'."
        )


def mol_to_canonical_dataframe(
    mol: 'Chem.Mol',
    *,
    structure_id: str,
    residue_name: Optional[str] = None,
    chain_id: str = 'L',
) -> pd.DataFrame:
    """Convert an RDKit molecule into the canonical structure DataFrame."""
    require_rdkit()

    if mol.GetNumConformers() == 0:
        raise ValueError("RDKit molecule must contain at least one conformer with coordinates")

    conformer = mol.GetConformer()
    if not conformer.Is3D():
        # Coordinates may still be 3D even if flag unset; proceed regardless
        pass

    name = mol.GetProp('_Name') if mol.HasProp('_Name') else structure_id
    res_name = residue_name or (name[:3].upper() if name else 'LIG')
    if len(res_name) < 3:
        res_name = res_name.ljust(3, '_')

    rows: List[dict] = []
    for atom_idx, atom in enumerate(mol.GetAtoms(), start=1):
        pos = conformer.GetAtomPosition(atom.GetIdx())
        element = atom.GetSymbol()
        atom_name = f"{element}{atom_idx}"[:4]
        row = {
            'structure_id': structure_id,
            'group': 'HETATM',
            'atom_id': atom_idx,
            'atom_name': atom_name,
            'element': element,
            'alt_id': '',
            'res_name': res_name,
            'res_name3l': res_name,
            'res_name1l': '',
            'auth_seq_id': 1,
            'label_seq_id': 1,
            'gen_seq_id': 1,
            'insertion': '',
            'grn': '',
            'auth_chain_id': chain_id,
            'label_chain_id': chain_id,
            'entity_id': 1,
            'x': float(pos.x),
            'y': float(pos.y),
            'z': float(pos.z),
            'occupancy': 1.0,
            'b_factor': 0.0,
            'charge': str(atom.GetFormalCharge()),
            'phi': np.nan,
            'psi': np.nan,
            'omega': np.nan,
            'model_num': 1,
            'res_id': f"{res_name}_{1}",
            'auth_comp_id': res_name,
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # Ensure all canonical columns exist with defaults
    for column, dtype in STRUCT_COLUMN_DTYPE.items():
        if column not in df.columns:
            if dtype in (int, 'int'):
                df[column] = pd.Series(pd.NA, index=df.index, dtype='Int64')
            elif dtype is float:
                df[column] = np.nan
            else:
                df[column] = ''

    df = df[[col for col in SORTED_STRUCT_COLUMNS if col in df.columns] +
            [col for col in df.columns if col not in SORTED_STRUCT_COLUMNS]]
    return df


def canonical_dataframe_to_mol(
    df: pd.DataFrame,
    *,
    metadata: Optional[Dict[str, Any]] = None,
) -> 'Chem.Mol':
    """Convert a canonical structure DataFrame back into an RDKit molecule."""

    require_rdkit()
    meta = metadata or {}

    mol_block = meta.get('mol_block')
    if mol_block:
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
        if mol is not None:
            try:
                conf = mol.GetConformer()
                coords = df.reset_index()[['x', 'y', 'z']].to_numpy()
                for idx, (x_val, y_val, z_val) in enumerate(coords):
                    conf.SetAtomPosition(idx, Point3D(float(x_val), float(y_val), float(z_val)))
            except Exception:
                pass
            return mol

    working = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df.copy()
    mol_edit = Chem.RWMol()
    conformer = Chem.Conformer(len(working))

    for idx, row in working.iterrows():
        symbol = (row.get('element') or '').strip()
        if not symbol:
            atom_name = str(row.get('atom_name', ''))
            symbol = atom_name[:1] if atom_name else 'C'

        atom = Chem.Atom(symbol)
        charge = row.get('charge')
        if charge not in (None, '', 'nan'):
            try:
                atom.SetFormalCharge(int(charge))
            except (ValueError, TypeError):
                pass
        atom_idx = mol_edit.AddAtom(atom)
        conformer.SetAtomPosition(atom_idx, Point3D(float(row['x']), float(row['y']), float(row['z'])))

    mol = mol_edit.GetMol()
    conformer.SetId(0)
    mol.AddConformer(conformer, assignId=True)
    return mol


def sdf_to_molecules(path: Path) -> Iterable['Chem.Mol']:
    """Load molecules from an SDF file using RDKit."""
    require_rdkit()
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for mol in supplier:
        if mol is None:
            continue
        yield mol


__all__ = [
    'mol_to_canonical_dataframe',
    'canonical_dataframe_to_mol',
    'sdf_to_molecules',
]

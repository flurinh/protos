"""
SDF file utilities for Protos.

This module provides utilities for working with Structure Data Format (SDF) files,
commonly used for storing molecular structures and associated data.
"""

import logging
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some SDF functionality will be limited.")
    HAS_RDKIT = False


def read_sdf_file(file_path: str, include_properties: bool = True,
                  sanitize: bool = True) -> List[Dict[str, Any]]:
    """
    Read molecules from an SDF file.
    
    Args:
        file_path: Path to the SDF file
        include_properties: Include molecular properties
        sanitize: Sanitize molecules (RDKit validation)
        
    Returns:
        List of dictionaries with molecule data
        
    Raises:
        ValueError: If RDKit is not available or file is invalid
    """
    if not HAS_RDKIT:
        raise ValueError("RDKit is required to read SDF files")
    
    molecules = []
    
    try:
        suppl = Chem.SDMolSupplier(str(file_path), sanitize=sanitize)
        
        for idx, mol in enumerate(suppl):
            if mol is None:
                logger.warning(f"Failed to read molecule at index {idx}")
                continue
            
            mol_data = {
                'index': idx,
                'smiles': Chem.MolToSmiles(mol),
                'mol': mol  # Store RDKit molecule object
            }
            
            # Extract properties from SDF
            if include_properties:
                props = {}
                for prop_name in mol.GetPropNames():
                    props[prop_name] = mol.GetProp(prop_name)
                mol_data['sdf_properties'] = props
            
            molecules.append(mol_data)
            
    except Exception as e:
        raise ValueError(f"Failed to read SDF file: {e}")
    
    return molecules


def write_sdf_file(file_path: str, molecules: List[Dict[str, Any]],
                   properties_to_write: Optional[List[str]] = None) -> None:
    """
    Write molecules to an SDF file.
    
    Args:
        file_path: Output file path
        molecules: List of molecule dictionaries
        properties_to_write: Specific properties to include (None = all)
        
    Raises:
        ValueError: If RDKit is not available or data is invalid
    """
    if not HAS_RDKIT:
        raise ValueError("RDKit is required to write SDF files")
    
    writer = Chem.SDWriter(str(file_path))
    
    try:
        for mol_data in molecules:
            # Get molecule object
            if 'mol' in mol_data and isinstance(mol_data['mol'], Chem.Mol):
                mol = mol_data['mol']
            elif 'smiles' in mol_data:
                mol = Chem.MolFromSmiles(mol_data['smiles'])
                if mol is None:
                    logger.warning(f"Failed to parse SMILES: {mol_data['smiles']}")
                    continue
            else:
                logger.warning("No molecule data found in entry")
                continue
            
            # Set properties
            if 'properties' in mol_data:
                for key, value in mol_data['properties'].items():
                    if properties_to_write is None or key in properties_to_write:
                        mol.SetProp(key, str(value))
            
            # Set SDF properties
            if 'sdf_properties' in mol_data:
                for key, value in mol_data['sdf_properties'].items():
                    if properties_to_write is None or key in properties_to_write:
                        mol.SetProp(key, str(value))
            
            # Add standard identifiers if available
            if 'chembl_id' in mol_data:
                mol.SetProp('CHEMBL_ID', mol_data['chembl_id'])
            if 'name' in mol_data:
                mol.SetProp('Name', mol_data['name'])
            
            writer.write(mol)
            
    finally:
        writer.close()


def sdf_to_dataframe(file_path: str, include_mol: bool = False) -> pd.DataFrame:
    """
    Convert SDF file to a pandas DataFrame.
    
    Args:
        file_path: Path to SDF file
        include_mol: Include RDKit molecule objects in DataFrame
        
    Returns:
        DataFrame with molecule data and properties
    """
    molecules = read_sdf_file(file_path)
    
    if not molecules:
        return pd.DataFrame()
    
    # Flatten the data for DataFrame
    rows = []
    for mol_data in molecules:
        row = {
            'smiles': mol_data['smiles'],
            'index': mol_data['index']
        }
        
        if include_mol:
            row['mol'] = mol_data.get('mol')
        
        # Add SDF properties
        if 'sdf_properties' in mol_data:
            row.update(mol_data['sdf_properties'])
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def dataframe_to_sdf(df: pd.DataFrame, file_path: str, 
                     smiles_col: str = 'smiles',
                     property_cols: Optional[List[str]] = None) -> None:
    """
    Write a DataFrame to an SDF file.
    
    Args:
        df: DataFrame with molecule data
        file_path: Output file path
        smiles_col: Column containing SMILES strings
        property_cols: Columns to include as properties (None = all except SMILES)
    """
    if smiles_col not in df.columns:
        raise ValueError(f"SMILES column '{smiles_col}' not found in DataFrame")
    
    molecules = []
    
    # Determine which columns to use as properties
    if property_cols is None:
        property_cols = [col for col in df.columns if col != smiles_col]
    
    for _, row in df.iterrows():
        mol_data = {
            'smiles': row[smiles_col],
            'properties': {}
        }
        
        # Add properties
        for col in property_cols:
            if col in row and pd.notna(row[col]):
                mol_data['properties'][col] = row[col]
        
        molecules.append(mol_data)
    
    write_sdf_file(file_path, molecules)


def merge_sdf_files(file_paths: List[str], output_path: str,
                    remove_duplicates: bool = True) -> int:
    """
    Merge multiple SDF files into one.
    
    Args:
        file_paths: List of input SDF files
        output_path: Output file path
        remove_duplicates: Remove duplicate molecules by SMILES
        
    Returns:
        Number of molecules written
    """
    all_molecules = []
    seen_smiles = set()
    
    for file_path in file_paths:
        try:
            molecules = read_sdf_file(file_path)
            
            for mol_data in molecules:
                smiles = mol_data['smiles']
                
                if remove_duplicates:
                    if smiles in seen_smiles:
                        continue
                    seen_smiles.add(smiles)
                
                all_molecules.append(mol_data)
                
        except Exception as e:
            logger.warning(f"Failed to read {file_path}: {e}")
    
    write_sdf_file(output_path, all_molecules)
    return len(all_molecules)


def validate_sdf_file(file_path: str) -> Tuple[bool, List[str]]:
    """
    Validate an SDF file.
    
    Args:
        file_path: Path to SDF file
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    if not HAS_RDKIT:
        return False, ["RDKit is required to validate SDF files"]
    
    errors = []
    
    try:
        suppl = Chem.SDMolSupplier(str(file_path), sanitize=False)
        
        mol_count = 0
        failed_count = 0
        
        for idx, mol in enumerate(suppl):
            if mol is None:
                failed_count += 1
                errors.append(f"Failed to read molecule at index {idx}")
            else:
                mol_count += 1
                
                # Try to sanitize
                try:
                    Chem.SanitizeMol(mol)
                except Exception as e:
                    errors.append(f"Molecule {idx} failed sanitization: {e}")
        
        if mol_count == 0:
            errors.append("No valid molecules found in file")
        
        if failed_count > 0:
            errors.append(f"{failed_count} molecules failed to load")
            
    except Exception as e:
        errors.append(f"Failed to open SDF file: {e}")
    
    return len(errors) == 0, errors


def extract_unique_properties(file_path: str) -> List[str]:
    """
    Extract all unique property names from an SDF file.
    
    Args:
        file_path: Path to SDF file
        
    Returns:
        List of unique property names
    """
    molecules = read_sdf_file(file_path, include_properties=True)
    
    all_props = set()
    for mol_data in molecules:
        if 'sdf_properties' in mol_data:
            all_props.update(mol_data['sdf_properties'].keys())
    
    return sorted(list(all_props))


def filter_sdf_by_property(file_path: str, output_path: str,
                          property_name: str, min_value: Optional[float] = None,
                          max_value: Optional[float] = None) -> int:
    """
    Filter SDF file by property value range.
    
    Args:
        file_path: Input SDF file
        output_path: Output SDF file
        property_name: Property to filter by
        min_value: Minimum value (inclusive)
        max_value: Maximum value (inclusive)
        
    Returns:
        Number of molecules written
    """
    molecules = read_sdf_file(file_path)
    filtered = []
    
    for mol_data in molecules:
        if 'sdf_properties' in mol_data and property_name in mol_data['sdf_properties']:
            try:
                value = float(mol_data['sdf_properties'][property_name])
                
                if min_value is not None and value < min_value:
                    continue
                if max_value is not None and value > max_value:
                    continue
                    
                filtered.append(mol_data)
            except ValueError:
                logger.warning(f"Non-numeric value for {property_name}")
    
    write_sdf_file(output_path, filtered)
    return len(filtered)
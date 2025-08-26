"""
Utility functions for ligand data loaders.

This module provides common utilities for handling ligand data, including:
- SMILES string validation and sanitization
- Identifier mapping between databases
- Molecular property calculations
"""

import re
import logging
from typing import Optional, Dict, List, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

# Try to import RDKit for molecular operations
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    from rdkit.Chem.inchi import MolToInchi, MolToInchiKey
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some ligand functionality will be limited.")
    HAS_RDKIT = False


def sanitize_smiles_filename(smiles: str) -> str:
    """
    Convert SMILES string to a safe filename.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Sanitized string safe for use as filename
        
    Example:
        >>> sanitize_smiles_filename("CC(=O)OC1=CC=CC=C1C(=O)O")
        'CC_eq_O_OC1_eq_CC_eq_CC_eq_C1C_eq_O_O'
    """
    # Define replacements for special characters
    replacements = {
        '=': '_eq_',
        '#': '_triple_',
        '(': '_lp_',
        ')': '_rp_',
        '[': '_lb_',
        ']': '_rb_',
        '+': '_plus_',
        '-': '_minus_',
        '@': '_at_',
        '/': '_fs_',
        '\\': '_bs_',
        '.': '_dot_',
        '%': '_pct_',
        '*': '_star_',
        '&': '_and_',
        '$': '_dollar_',
        '!': '_excl_',
        '?': '_quest_',
        ':': '_colon_',
        ';': '_semi_',
        ',': '_comma_',
        '<': '_lt_',
        '>': '_gt_',
        '{': '_lcb_',
        '}': '_rcb_',
        '|': '_pipe_',
        '~': '_tilde_',
        '^': '_caret_',
        ' ': '_',
        '\t': '_tab_',
        '\n': '_nl_',
        '\r': '_cr_'
    }
    
    # Apply replacements
    safe_name = smiles
    for char, replacement in replacements.items():
        safe_name = safe_name.replace(char, replacement)
    
    # Remove any remaining non-alphanumeric characters
    safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', safe_name)
    
    # Ensure it doesn't start with a number
    if safe_name and safe_name[0].isdigit():
        safe_name = 'mol_' + safe_name
    
    # Limit length to avoid filesystem issues
    if len(safe_name) > 200:
        # Use InChI key if available and SMILES is too long
        if HAS_RDKIT:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    safe_name = MolToInchiKey(mol)
            except:
                # Fallback to truncation
                safe_name = safe_name[:200]
        else:
            safe_name = safe_name[:200]
    
    return safe_name


def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a SMILES string.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        Tuple of (is_valid, canonical_smiles)
    """
    if not HAS_RDKIT:
        logger.warning("RDKit not available, cannot validate SMILES")
        return True, smiles  # Assume valid if we can't check
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, None
        
        # Get canonical SMILES
        canonical = Chem.MolToSmiles(mol, canonical=True)
        return True, canonical
    except Exception as e:
        logger.debug(f"SMILES validation failed for '{smiles}': {e}")
        return False, None


def smiles_to_inchi(smiles: str) -> Optional[Dict[str, str]]:
    """
    Convert SMILES to InChI and InChI Key.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary with 'inchi' and 'inchi_key', or None if conversion fails
    """
    if not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        inchi = MolToInchi(mol)
        inchi_key = MolToInchiKey(mol)
        
        return {
            'inchi': inchi,
            'inchi_key': inchi_key
        }
    except Exception as e:
        logger.debug(f"InChI conversion failed for '{smiles}': {e}")
        return None


def calculate_molecular_properties(smiles: str) -> Optional[Dict[str, float]]:
    """
    Calculate basic molecular properties for filtering.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of molecular properties, or None if calculation fails
    """
    if not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        properties = {
            'mw': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'tpsa': Descriptors.TPSA(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'heavy_atoms': mol.GetNumHeavyAtoms(),
            'formal_charge': Chem.rdmolops.GetFormalCharge(mol)
        }
        
        # Add Lipinski properties
        properties['lipinski_violations'] = sum([
            properties['mw'] > 500,
            properties['logp'] > 5,
            properties['hba'] > 10,
            properties['hbd'] > 5
        ])
        
        return properties
    except Exception as e:
        logger.debug(f"Property calculation failed for '{smiles}': {e}")
        return None


def is_drug_like(smiles: str, strict: bool = False) -> bool:
    """
    Check if a molecule satisfies drug-like criteria.
    
    Args:
        smiles: SMILES string
        strict: If True, apply stricter criteria
        
    Returns:
        True if drug-like, False otherwise
    """
    props = calculate_molecular_properties(smiles)
    if not props:
        return False
    
    if strict:
        # Stricter criteria
        return all([
            150 <= props['mw'] <= 450,
            -0.4 <= props['logp'] <= 4.5,
            props['hba'] <= 8,
            props['hbd'] <= 4,
            props['tpsa'] <= 120,
            props['rotatable_bonds'] <= 8,
            props['lipinski_violations'] == 0
        ])
    else:
        # Standard Lipinski Rule of Five
        return props['lipinski_violations'] <= 1


def extract_protein_mapping(identifier: str) -> Dict[str, str]:
    """
    Extract and map various protein identifier formats.
    
    Args:
        identifier: Protein identifier (UniProt, PDB, gene name, etc.)
        
    Returns:
        Dictionary with detected identifier types and values
    """
    mapping = {'original': identifier}
    
    # Clean the identifier
    identifier = identifier.strip().upper()
    
    # PDB ID pattern (4 characters, alphanumeric)
    if re.match(r'^[0-9][A-Z0-9]{3}$', identifier):
        mapping['pdb_id'] = identifier
        mapping['type'] = 'pdb'
    
    # UniProt ID patterns
    elif re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', identifier):
        mapping['uniprot_id'] = identifier
        mapping['type'] = 'uniprot_acc'
    
    # UniProt mnemonic (e.g., EGFR_HUMAN)
    elif re.match(r'^[A-Z0-9]+_[A-Z]+$', identifier):
        mapping['uniprot_mnemonic'] = identifier
        mapping['type'] = 'uniprot_name'
        # Extract gene name
        mapping['gene_name'] = identifier.split('_')[0]
    
    # ChEMBL target ID
    elif identifier.startswith('CHEMBL'):
        mapping['chembl_id'] = identifier
        mapping['type'] = 'chembl'
    
    # Gene name (simple alphanumeric)
    elif re.match(r'^[A-Z][A-Z0-9]+$', identifier):
        mapping['gene_name'] = identifier
        mapping['type'] = 'gene'
    
    else:
        mapping['type'] = 'unknown'
    
    return mapping


def parse_activity_value(value: str, units: str) -> Optional[float]:
    """
    Parse activity values to nanomolar (nM) units.
    
    Args:
        value: Activity value as string
        units: Units of the activity
        
    Returns:
        Activity value in nM, or None if parsing fails
    """
    try:
        # Extract numeric value
        value_str = re.sub(r'[<>~=]', '', str(value)).strip()
        numeric_value = float(value_str)
        
        # Convert to nM based on units
        units_lower = units.lower()
        
        if 'nm' in units_lower:
            return numeric_value
        elif 'um' in units_lower or 'μm' in units_lower:
            return numeric_value * 1000
        elif 'mm' in units_lower:
            return numeric_value * 1000000
        elif 'm' in units_lower and 'nm' not in units_lower and 'um' not in units_lower:
            # Molar
            return numeric_value * 1e9
        else:
            logger.warning(f"Unknown units: {units}")
            return None
            
    except (ValueError, AttributeError) as e:
        logger.debug(f"Failed to parse activity value '{value}' with units '{units}': {e}")
        return None


def batch_download_helper(download_func, items: List, max_workers: int = 4) -> Tuple[List, List]:
    """
    Helper function for parallel downloads with error handling.
    
    Args:
        download_func: Function to call for each item
        items: List of items to download
        max_workers: Number of parallel workers
        
    Returns:
        Tuple of (successful_items, failed_items)
    """
    import concurrent.futures
    
    successful = []
    failed = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_item = {executor.submit(download_func, item): item for item in items}
        
        # Process completed tasks
        for future in concurrent.futures.as_completed(future_to_item):
            item = future_to_item[future]
            try:
                result = future.result()
                if result:
                    successful.append(item)
                else:
                    failed.append(item)
            except Exception as e:
                logger.error(f"Download failed for {item}: {e}")
                failed.append(item)
    
    return successful, failed


def create_sdf_from_smiles(smiles: str, properties: Optional[Dict] = None) -> Optional[str]:
    """
    Create SDF format string from SMILES.
    
    Args:
        smiles: SMILES string
        properties: Optional dictionary of properties to include
        
    Returns:
        SDF format string, or None if conversion fails
    """
    if not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Add properties if provided
        if properties:
            for key, value in properties.items():
                mol.SetProp(str(key), str(value))
        
        # Generate SDF
        sdf_string = Chem.MolToMolBlock(mol)
        return sdf_string + "$$$$\n"
        
    except Exception as e:
        logger.debug(f"SDF creation failed for '{smiles}': {e}")
        return None


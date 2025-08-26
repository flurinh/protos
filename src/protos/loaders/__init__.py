"""
Protos data loaders.

This module provides loaders for various biological data sources.
"""

from protos.loaders.uniprot_loader import UniprotDL
# from protos.loaders.gpcrdb_loader import GPCRdbDL  # Temporarily disabled due to import issues
from protos.loaders.download_structures import (
    download_protein_structures,
    download_structures_with_processor
)
from protos.loaders import chembl_loader
from protos.loaders.ligand_utils import (
    sanitize_smiles_filename,
    validate_smiles,
    smiles_to_inchi,
    calculate_molecular_properties,
    is_drug_like,
    parse_activity_value,
    create_sdf_from_smiles,
    extract_protein_mapping
)

# Import new ligand database loaders
from protos.loaders import qm9_loader
from protos.loaders import ccd_loader
from protos.loaders import enamine_loader

__all__ = [
    'UniprotDL',
    # 'GPCRdbDL',  # Temporarily disabled
    'download_protein_structures',
    'download_structures_with_processor',
    'chembl_loader',
    'sanitize_smiles_filename',
    'validate_smiles',
    'smiles_to_inchi',
    'calculate_molecular_properties',
    'is_drug_like',
    'parse_activity_value',
    'create_sdf_from_smiles',
    'extract_protein_mapping',
    # New loaders
    'qm9_loader',
    'ccd_loader', 
    'enamine_loader'
]
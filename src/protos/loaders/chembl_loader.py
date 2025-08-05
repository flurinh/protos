"""
ChEMBL data loader for Protos - REFACTORED to follow proper patterns.

This loader provides utility functions for downloading ligand and bioactivity data 
from ChEMBL database. It follows Protos principles:
- NO path management
- NO directory creation
- NO entity registry management
- Pure utility functions only

The loader should be used by processors, not directly by users.
"""

import json
import logging
from typing import List, Dict, Optional, Union, Tuple
from datetime import datetime

from protos.loaders.ligand_utils import (
    extract_protein_mapping, parse_activity_value
)

logger = logging.getLogger(__name__)

# ChEMBL client will be imported on demand
HAS_CHEMBL = None
chembl_client = None

def _get_chembl_client():
    """Lazy import of ChEMBL client."""
    global HAS_CHEMBL, chembl_client
    
    if HAS_CHEMBL is None:
        try:
            from chembl_webresource_client.new_client import new_client as _chembl_client
            chembl_client = _chembl_client
            HAS_CHEMBL = True
        except (ImportError, Exception) as e:
            logger.warning(f"chembl_webresource_client not available: {e}")
            HAS_CHEMBL = False
            chembl_client = None
    
    return chembl_client


# Common gene name aliases that map directly to ChEMBL IDs
GENE_ALIASES = {
    'EGFR': 'CHEMBL203',      # Epidermal growth factor receptor
    'HER2': 'CHEMBL1824',     # ERBB2
    'ERBB2': 'CHEMBL1824',    # HER2
    'COX2': 'CHEMBL230',      # PTGS2
    'PTGS2': 'CHEMBL230',     # COX2
    'ABL': 'CHEMBL1862',      # ABL1
    'ABL1': 'CHEMBL1862',     # BCR-ABL
    'SRC': 'CHEMBL267',       # SRC kinase
    'VEGFR2': 'CHEMBL279',    # KDR
    'KDR': 'CHEMBL279',       # VEGFR2
    'MET': 'CHEMBL3717',      # c-Met
    'ALK': 'CHEMBL4247',      # Anaplastic lymphoma kinase
    'BTK': 'CHEMBL5251',      # Bruton's tyrosine kinase
    'JAK2': 'CHEMBL6133',     # Janus kinase 2
    'CDK4': 'CHEMBL331',      # Cyclin-dependent kinase 4
    'CDK6': 'CHEMBL2508',     # Cyclin-dependent kinase 6
    'BRAF': 'CHEMBL5145',     # B-Raf
    'PI3K': 'CHEMBL3145',     # PIK3CA
    'PIK3CA': 'CHEMBL3145',   # PI3K alpha
    'MTOR': 'CHEMBL2185',     # mTOR
    'HDAC': 'CHEMBL325',      # HDAC1
    'PD1': 'CHEMBL4683',      # PDCD1
    'PDL1': 'CHEMBL5917',     # CD274
    'CTLA4': 'CHEMBL4301',    # CTLA-4
}

# Common PDB to ChEMBL target mappings
PDB_ALIASES = {
    # EGFR structures
    '1M17': 'CHEMBL203',      # EGFR kinase domain
    '2ITY': 'CHEMBL203',      # EGFR with erlotinib
    '4HJO': 'CHEMBL203',      # EGFR T790M with inhibitor
    
    # Other kinases
    '1T46': 'CHEMBL267',      # SRC kinase
    '3LXK': 'CHEMBL1862',     # ABL1 with dasatinib
    '1M52': 'CHEMBL331',      # CDK4
    
    # GPCRs
    '3EML': 'CHEMBL213',      # Adenosine A2A receptor
    '4EJ4': 'CHEMBL210',      # β2 adrenergic receptor
    '5C1M': 'CHEMBL1889',     # μ-opioid receptor
    
    # Other targets
    '1CX2': 'CHEMBL230',      # COX-2
    '4B0Q': 'CHEMBL2185',     # mTOR
}


def map_protein_to_chembl_target(protein_id: str, 
                                protein_mapping_cache: Optional[Dict] = None) -> Optional[str]:
    """
    Map protein identifier to ChEMBL target ID.
    
    This is a pure utility function that does the mapping logic.
    Cache management should be handled by the caller (processor).
    
    Args:
        protein_id: Protein identifier (UniProt, PDB, gene name)
        protein_mapping_cache: Optional cache dictionary
        
    Returns:
        ChEMBL target ID, or None if not found
    """
    chembl = _get_chembl_client()
    if not chembl:
        logger.error("ChEMBL client not available")
        return None
    
    # Check cache if provided
    if protein_mapping_cache and protein_id in protein_mapping_cache:
        return protein_mapping_cache[protein_id]
    
    # Check gene aliases
    if protein_id.upper() in GENE_ALIASES:
        return GENE_ALIASES[protein_id.upper()]
    
    # Check PDB aliases
    if protein_id.upper() in PDB_ALIASES:
        return PDB_ALIASES[protein_id.upper()]
    
    # Extract identifier type
    id_info = extract_protein_mapping(protein_id)
    
    try:
        # Search based on identifier type
        if id_info['type'] == 'uniprot_acc':
            targets = chembl.target.filter(
                target_components__accession=id_info['uniprot_id']
            )
        elif id_info['type'] == 'uniprot_name':
            # Try gene name first
            gene = id_info.get('gene_name', '')
            targets = chembl.target.filter(pref_name__iexact=gene)
            
            # If no results, try without exact match
            if not list(targets):
                targets = chembl.target.filter(pref_name__icontains=gene)
        elif id_info['type'] == 'pdb':
            # PDB codes are harder to map directly in ChEMBL
            # Try to search by PDB code in target name/description
            pdb_code = id_info['pdb_id'].upper()
            
            # First try exact match in pref_name
            targets = chembl.target.filter(pref_name__icontains=pdb_code)
            
            # If no results, try searching in synonyms
            if not list(targets):
                targets = chembl.target.filter(synonyms__icontains=pdb_code)
            
            # If still no results, we might need to use a different approach
            # For now, log a warning
            if not list(targets):
                logger.warning(f"Could not find ChEMBL target for PDB {pdb_code}")
        elif id_info['type'] == 'chembl':
            # Already a ChEMBL ID
            return id_info['chembl_id']
        else:
            # Try as gene name
            targets = chembl.target.filter(pref_name__iexact=protein_id)
            
            # If no results, try without exact match
            if not list(targets):
                targets = chembl.target.filter(pref_name__icontains=protein_id)
        
        # Get best target - prefer SINGLE PROTEIN over others
        best_target = None
        for target in targets:
            if target['target_chembl_id']:
                # Prefer SINGLE PROTEIN targets
                if target.get('target_type') == 'SINGLE PROTEIN':
                    return target['target_chembl_id']
                elif best_target is None:
                    best_target = target
        
        # If no single protein found, use the first valid target
        if best_target and best_target['target_chembl_id']:
            return best_target['target_chembl_id']
                
    except Exception as e:
        logger.error(f"Failed to map {protein_id} to ChEMBL: {e}")
    
    return None


def query_protein_ligands(protein_id: str, 
                         activity_types: Optional[List[str]] = None,
                         min_pchembl: float = 5.0,
                         max_value_nm: float = 10000,
                         limit: Optional[int] = None) -> List[Dict]:
    """
    Query ligands for a protein target from ChEMBL.
    
    This is a pure utility function that returns raw data.
    Storage and entity management should be handled by the caller.
    
    Args:
        protein_id: Protein identifier (UniProt/PDB/gene name)
        activity_types: Types of bioactivities to retrieve (default: IC50, Ki, Kd)
        min_pchembl: Minimum pChEMBL value (negative log of activity)
        max_value_nm: Maximum activity value in nM
        limit: Maximum number of compounds
        
    Returns:
        List of ligand data dictionaries
    """
    chembl = _get_chembl_client()
    if not chembl:
        logger.error("ChEMBL client not available")
        return []
    
    # Map to ChEMBL target
    chembl_target = map_protein_to_chembl_target(protein_id)
    if not chembl_target:
        logger.warning(f"Could not map {protein_id} to ChEMBL target")
        return []
    
    # Default activity types
    if activity_types is None:
        activity_types = ['IC50', 'Ki', 'Kd', 'EC50']
    
    ligand_data = []
    
    try:
        # Query activities
        activities = chembl.activity.filter(
            target_chembl_id=chembl_target,
            type__in=activity_types,
            pchembl_value__gte=min_pchembl,
            assay_type='B'  # Binding assays
        )
        
        count = 0
        
        for activity in activities:
            # Apply limit if set
            if limit and count >= limit:
                break
            
            # Extract compound info
            compound_id = activity.get('molecule_chembl_id')
            if not compound_id:
                continue
            
            # Get compound details
            try:
                compound = chembl.molecule.get(compound_id)
                smiles = compound.get('molecule_structures', {}).get('canonical_smiles')
                
                if not smiles:
                    continue
                
                # Parse activity value
                value = activity.get('value')
                units = activity.get('units', 'nM')
                value_nm = parse_activity_value(value, units)
                
                if value_nm and value_nm <= max_value_nm:
                    ligand_info = {
                        'smiles': smiles,
                        'chembl_id': compound_id,
                        'activity_type': activity.get('type'),
                        'value': value,
                        'units': units,
                        'value_nm': value_nm,
                        'pchembl_value': activity.get('pchembl_value'),
                        'assay_id': activity.get('assay_chembl_id'),
                        'target_id': chembl_target,
                        'protein_id': protein_id
                    }
                    
                    ligand_data.append(ligand_info)
                    count += 1
                    
            except Exception as e:
                logger.debug(f"Failed to process compound {compound_id}: {e}")
                continue
        
        logger.info(f"Retrieved {len(ligand_data)} ligands for {protein_id}")
        
    except Exception as e:
        logger.error(f"Failed to query ligands for {protein_id}: {e}")
    
    return ligand_data


def get_compound_by_chembl_id(chembl_id: str) -> Optional[Dict]:
    """
    Get compound information by ChEMBL ID.
    
    Args:
        chembl_id: ChEMBL compound ID
        
    Returns:
        Dictionary with compound information, or None if not found
    """
    chembl = _get_chembl_client()
    if not chembl:
        logger.error("ChEMBL client not available")
        return None
    
    try:
        compound = chembl.molecule.get(chembl_id)
        smiles = compound.get('molecule_structures', {}).get('canonical_smiles')
        
        if smiles:
            return {
                'smiles': smiles,
                'chembl_id': chembl_id,
                'name': compound.get('pref_name', ''),
                'max_phase': compound.get('max_phase', 0),
                'therapeutic_flags': compound.get('therapeutic_flags', []),
                'molecule_type': compound.get('molecule_type', '')
            }
    except Exception as e:
        logger.error(f"Failed to get compound {chembl_id}: {e}")
    
    return None


def search_similar_compounds_chembl(smiles: str, similarity: float = 0.8, 
                                  limit: int = 100) -> List[Dict]:
    """
    Search for similar compounds in ChEMBL.
    
    Args:
        smiles: Query SMILES string
        similarity: Similarity threshold (0-1)
        limit: Maximum number of results
        
    Returns:
        List of similar compound dictionaries
    """
    chembl = _get_chembl_client()
    if not chembl:
        logger.error("ChEMBL client not available")
        return []
    
    try:
        # Search similar molecules
        similar = chembl.similarity.filter(
            smiles=smiles,
            similarity=int(similarity * 100)  # ChEMBL uses percentage
        ).only(['molecule_chembl_id', 'similarity', 'molecule_structures'])[:limit]
        
        results = []
        for mol in similar:
            mol_smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
            if mol_smiles:
                result = {
                    'chembl_id': mol['molecule_chembl_id'],
                    'smiles': mol_smiles,
                    'similarity': mol['similarity'] / 100.0
                }
                results.append(result)
        
        return results
        
    except Exception as e:
        logger.error(f"Similarity search failed: {e}")
        return []
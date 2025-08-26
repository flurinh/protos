"""
Structure-ligand analysis functions for Protos.

This module provides functions for analyzing protein-ligand interactions
using data from CifBaseProcessor without modifying the processor itself.
These functions can be used standalone or integrated into workflows.
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix, cKDTree
from typing import Dict, List, Tuple, Optional, Set, Union
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Try to import RDKit for molecule handling
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some conversion functions will be limited.")
    HAS_RDKIT = False


def extract_all_ligands(cif_processor, pdb_id: str, 
                       exclude_common: bool = True,
                       min_atoms: int = 3) -> List[Dict]:
    """
    Extract all ligands from a structure using CifBaseProcessor data.
    
    Args:
        cif_processor: CifBaseProcessor instance with loaded data
        pdb_id: Structure identifier
        exclude_common: Exclude water, ions, and common molecules
        min_atoms: Minimum number of atoms for a ligand
        
    Returns:
        List of ligand dictionaries with structure and metadata
    """
    # Common molecules to exclude
    common_hetero = {'HOH', 'WAT', 'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 
                     'SO4', 'PO4', 'IOD', 'BR', 'F', 'FE', 'CU', 'MN',
                     'CD', 'NI', 'HG', 'CO', 'SR', 'BA', 'RB', 'CS',
                     'GOL', 'EDO', 'PEG', 'ACT', 'DMS', 'FMT'}
    
    # Get structure data
    structure_data = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id].copy()
    hetatoms = structure_data[structure_data['group'] == 'HETATM']
    
    ligands = []
    
    # Group by residue
    for (res_name, chain_id, res_id), atoms in hetatoms.groupby(
        ['res_name3l', 'auth_chain_id', 'auth_seq_id']
    ):
        # Apply filters
        if exclude_common and res_name in common_hetero:
            continue
        if len(atoms) < min_atoms:
            continue
            
        # Calculate properties
        coords = atoms[['x', 'y', 'z']].values
        centroid = coords.mean(axis=0)
        
        # Create unique identifier
        ligand_id = f"{pdb_id}_{chain_id}_{res_name}_{res_id}"
        
        ligand_info = {
            'ligand_id': ligand_id,
            'pdb_id': pdb_id,
            'res_name3l': res_name,
            'chain_id': chain_id,
            'res_id': res_id,
            'num_atoms': len(atoms),
            'atoms': atoms,
            'centroid': centroid,
            'coords': coords,
            'atom_names': atoms['atom_name'].tolist(),
            'elements': atoms['atom_name'].str.extract(r'^([A-Z][a-z]?)')[0].tolist()
        }
        
        ligands.append(ligand_info)
        
    logger.info(f"Found {len(ligands)} ligands in {pdb_id}")
    return ligands


def get_ligand_by_id(cif_processor, pdb_id: str, res_name: str,
                     chain_id: Optional[str] = None,
                     res_id: Optional[int] = None) -> Optional[pd.DataFrame]:
    """
    Get specific ligand atoms from structure.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        res_name: Three-letter ligand code (e.g., 'ATP')
        chain_id: Optional chain specification
        res_id: Optional residue number
        
    Returns:
        DataFrame of ligand atoms or None if not found
    """
    structure_data = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id]
    
    # Build filter
    ligand_filter = (structure_data['group'] == 'HETATM') & \
                   (structure_data['res_name3l'] == res_name)
    
    if chain_id is not None:
        ligand_filter &= (structure_data['auth_chain_id'] == chain_id)
    if res_id is not None:
        ligand_filter &= (structure_data['auth_seq_id'] == res_id)
        
    ligand_atoms = structure_data[ligand_filter]
    
    if ligand_atoms.empty:
        return None
        
    return ligand_atoms.copy()


def get_binding_site(cif_processor, pdb_id: str, ligand_atoms: pd.DataFrame,
                     cutoff: float = 5.0, 
                     include_backbone: bool = True) -> Dict[str, pd.DataFrame]:
    """
    Get binding site residues around a ligand.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        ligand_atoms: DataFrame of ligand atoms
        cutoff: Distance cutoff in Angstroms
        include_backbone: Include backbone atoms in analysis
        
    Returns:
        Dictionary with 'residues' DataFrame and 'atoms' DataFrame
    """
    # Get protein atoms
    structure_data = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id]
    protein_atoms = structure_data[structure_data['group'] == 'ATOM'].copy()
    
    if not include_backbone:
        # Exclude backbone atoms (N, CA, C, O)
        protein_atoms = protein_atoms[~protein_atoms['atom_name'].isin(['N', 'CA', 'C', 'O'])]
    
    # Build KDTree for efficient neighbor search
    if protein_atoms.empty:
        return {'residues': pd.DataFrame(), 'atoms': pd.DataFrame()}
        
    protein_tree = cKDTree(protein_atoms[['x', 'y', 'z']].values)
    ligand_coords = ligand_atoms[['x', 'y', 'z']].values
    
    # Find all protein atoms within cutoff
    nearby_indices = []
    for coord in ligand_coords:
        indices = protein_tree.query_ball_point(coord, cutoff)
        nearby_indices.extend(indices)
    
    nearby_indices = list(set(nearby_indices))
    
    if not nearby_indices:
        return {'residues': pd.DataFrame(), 'atoms': pd.DataFrame()}
        
    # Get the atoms
    binding_atoms = protein_atoms.iloc[nearby_indices].copy()
    
    # Calculate distances
    binding_atoms['min_distance'] = np.inf
    for idx, atom in binding_atoms.iterrows():
        atom_coord = atom[['x', 'y', 'z']].values.reshape(1, -1)
        distances = distance_matrix(atom_coord, ligand_coords)
        binding_atoms.loc[idx, 'min_distance'] = distances.min()
    
    # Group by residue
    residue_data = []
    for (chain, res_id, res_name), res_atoms in binding_atoms.groupby(
        ['auth_chain_id', 'auth_seq_id', 'res_name3l']
    ):
        residue_data.append({
            'chain_id': chain,
            'res_id': res_id,
            'res_name': res_name,
            'res_name1l': res_atoms.iloc[0]['res_name1l'] if 'res_name1l' in res_atoms.columns else '',
            'min_distance': res_atoms['min_distance'].min(),
            'num_atoms': len(res_atoms),
            'contact_atoms': res_atoms['atom_name'].tolist()
        })
    
    residues_df = pd.DataFrame(residue_data)
    
    # Sort by distance
    if not residues_df.empty:
        residues_df = residues_df.sort_values('min_distance')
    
    return {
        'residues': residues_df,
        'atoms': binding_atoms
    }


def calculate_ligand_interactions(cif_processor, pdb_id: str, 
                                ligand_atoms: pd.DataFrame,
                                detailed: bool = True) -> Dict:
    """
    Calculate detailed interactions between ligand and protein.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        ligand_atoms: DataFrame of ligand atoms
        detailed: Whether to include detailed interaction lists
        
    Returns:
        Dictionary of interaction analysis
    """
    from protos.processing.structure.ligand_interactions import LigandInteractionAnalyzer
    
    # Get structure data
    structure_data = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id]
    
    # Initialize analyzer
    analyzer = LigandInteractionAnalyzer(structure_data)
    
    # Get interactions
    interactions = analyzer.get_all_interactions(ligand_atoms)
    
    # Add binding site analysis
    binding_site = get_binding_site(cif_processor, pdb_id, ligand_atoms)
    interactions['binding_site'] = {
        'num_residues': len(binding_site['residues']),
        'residue_list': binding_site['residues']['res_name'].tolist() if not binding_site['residues'].empty else []
    }
    
    if not detailed:
        # Return only summary
        return interactions['summary']
        
    return interactions


def ligand_to_smiles(ligand_atoms: pd.DataFrame, 
                    ligand_name: str,
                    use_template: bool = True,
                    use_pdb_ccd: bool = True) -> Optional[str]:
    """
    Convert ligand coordinates to SMILES string.
    
    Args:
        ligand_atoms: DataFrame with ligand atom coordinates
        ligand_name: Three-letter ligand code
        use_template: Try to use template from ligand databases
        use_pdb_ccd: Try to fetch from PDB Chemical Component Dictionary
        
    Returns:
        SMILES string or None if conversion fails
    """
    if not HAS_RDKIT:
        logger.warning("RDKit required for SMILES conversion")
        return None
        
    # First try template matching for common ligands
    if use_template:
        template_smiles = get_ligand_template_smiles(ligand_name)
        if template_smiles:
            return template_smiles
    
    # Try PDB Chemical Component Dictionary
    if use_pdb_ccd:
        ccd_smiles = fetch_smiles_from_pdb_ccd(ligand_name)
        if ccd_smiles:
            return ccd_smiles
    
    # Otherwise try to reconstruct from coordinates
    # This is complex and would require bond inference
    logger.warning(f"Coordinate-based SMILES reconstruction not implemented for {ligand_name}")
    return None


def get_ligand_template_smiles(ligand_name: str) -> Optional[str]:
    """
    Get SMILES template for known ligands.
    
    This would connect to databases like:
    - PDB Chemical Component Dictionary
    - ChEMBL
    - PubChem
    """
    # Common ligands for demonstration
    known_ligands = {
        'ATP': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N',
        'ADP': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N',
        'AMP': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N',
        'GTP': 'C1=NC2=C(N1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N',
        'NAD': 'C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N',
        'FAD': 'CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)CC(C(C(COP(=O)(O)OP(=O)(O)OCC4C(C(C(O4)N5C=NC6=C(N=CN=C65)N)O)O)O)O)O',
        'HEM': 'CC1=C(C2=CC3=C(C(=C(N3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1N2)C)C=C)C)C=C)C)CCC(=O)O)CCC(=O)O',
        'SEP': 'C(C(C(=O)O)N)OP(=O)(O)O',  # Phosphoserine
        'TPO': 'CC(C)(COP(=O)(O)O)C(C(=O)NCCC(=O)O)O',  # Phosphothreonine
        'PTR': 'CC1=CC=C(C=C1)OP(=O)(O)O',  # Phosphotyrosine
        'SO4': 'O=S(=O)([O-])[O-]',  # Sulfate
        'PO4': 'O=P([O-])([O-])[O-]',  # Phosphate
        'CIT': 'C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O',  # Citrate
        'GOL': 'C(CO)O',  # Glycerol
        'PEG': 'C(CO)O'  # Polyethylene glycol fragment
    }
    
    return known_ligands.get(ligand_name)


def export_ligand_sdf(cif_processor, pdb_id: str, ligand_atoms: pd.DataFrame,
                     output_path: str, ligand_name: str = 'LIG',
                     include_hydrogens: bool = False) -> bool:
    """
    Export ligand in SDF format with 3D coordinates.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        ligand_atoms: DataFrame of ligand atoms
        output_path: Path to save SDF file
        ligand_name: Ligand name for the file
        include_hydrogens: Whether to include hydrogen atoms
        
    Returns:
        Success status
    """
    if not HAS_RDKIT:
        logger.error("RDKit required for SDF export")
        return False
        
    try:
        # Try to create molecule with proper connectivity
        mol = create_mol_from_atoms(ligand_atoms, ligand_name)
        if mol:
            # Successfully created with connectivity
            pass
        else:
            # Fall back to template-based approach
            smiles = ligand_to_smiles(ligand_atoms, ligand_name)
            if not smiles:
                logger.error(f"Could not determine structure for {ligand_name}")
                return False
                
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False
                
            # Map template atoms to structure coordinates
            mol = map_coords_to_mol(mol, ligand_atoms)
        
        # Add properties
        mol.SetProp("_Name", f"{pdb_id}_{ligand_name}")
        mol.SetProp("PDB_ID", pdb_id)
        mol.SetProp("LIGAND_NAME", ligand_name)
        mol.SetProp("NUM_ATOMS", str(len(ligand_atoms)))
        
        # Add structure metadata
        if 'auth_chain_id' in ligand_atoms.columns:
            chains = ligand_atoms['auth_chain_id'].unique()
            mol.SetProp("CHAIN_ID", ','.join(chains))
        
        # Write SDF
        writer = Chem.SDWriter(output_path)
        writer.write(mol)
        writer.close()
        
        logger.info(f"Exported {ligand_name} to {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to export SDF: {e}")
        return False


def compare_ligand_binding_sites(cif_processor, 
                                ligand_list: List[Tuple[str, str, Optional[str]]],
                                cutoff: float = 5.0) -> pd.DataFrame:
    """
    Compare binding sites across multiple ligand-structure pairs.
    
    Args:
        cif_processor: CifBaseProcessor instance
        ligand_list: List of (pdb_id, ligand_name, chain_id) tuples
        cutoff: Distance cutoff for binding site
        
    Returns:
        Similarity matrix as DataFrame
    """
    binding_sites = {}
    
    # Get binding sites for each ligand
    for pdb_id, ligand_name, chain_id in ligand_list:
        ligand_atoms = get_ligand_by_id(cif_processor, pdb_id, ligand_name, chain_id)
        if ligand_atoms is None:
            continue
            
        site = get_binding_site(cif_processor, pdb_id, ligand_atoms, cutoff)
        if not site['residues'].empty:
            # Create residue signature
            residue_set = set(
                f"{r['res_name']}" 
                for _, r in site['residues'].iterrows()
            )
            binding_sites[f"{pdb_id}_{ligand_name}"] = residue_set
    
    # Calculate Jaccard similarity between binding sites
    site_names = list(binding_sites.keys())
    n_sites = len(site_names)
    similarity_matrix = np.zeros((n_sites, n_sites))
    
    for i, site1 in enumerate(site_names):
        for j, site2 in enumerate(site_names):
            if i == j:
                similarity_matrix[i, j] = 1.0
            else:
                set1 = binding_sites[site1]
                set2 = binding_sites[site2]
                
                if len(set1) == 0 or len(set2) == 0:
                    similarity_matrix[i, j] = 0.0
                else:
                    jaccard = len(set1 & set2) / len(set1 | set2)
                    similarity_matrix[i, j] = jaccard
                    similarity_matrix[j, i] = jaccard
    
    return pd.DataFrame(similarity_matrix, index=site_names, columns=site_names)


def find_conserved_interactions(cif_processor,
                               ligand_list: List[Tuple[str, str, Optional[str]]],
                               min_conservation: float = 0.8) -> Dict:
    """
    Find interactions conserved across multiple structures.
    
    Args:
        cif_processor: CifBaseProcessor instance
        ligand_list: List of (pdb_id, ligand_name, chain_id) tuples
        min_conservation: Minimum fraction of structures with interaction
        
    Returns:
        Dictionary of conserved interactions
    """
    all_interactions = []
    
    # Get interactions for each structure
    for pdb_id, ligand_name, chain_id in ligand_list:
        ligand_atoms = get_ligand_by_id(cif_processor, pdb_id, ligand_name, chain_id)
        if ligand_atoms is None:
            continue
            
        interactions = calculate_ligand_interactions(cif_processor, pdb_id, ligand_atoms, detailed=True)
        all_interactions.append({
            'structure': f"{pdb_id}_{ligand_name}",
            'interactions': interactions
        })
    
    n_structures = len(all_interactions)
    if n_structures < 2:
        return {}
    
    # Analyze conservation
    conserved = {
        'binding_residues': {},
        'interaction_types': {}
    }
    
    # Track residue occurrence
    residue_counts = {}
    for struct_data in all_interactions:
        binding_res = struct_data['interactions'].get('binding_residues', [])
        for res in binding_res:
            res_key = f"{res['res_name']}"
            residue_counts[res_key] = residue_counts.get(res_key, 0) + 1
    
    # Find conserved residues
    for res_key, count in residue_counts.items():
        conservation = count / n_structures
        if conservation >= min_conservation:
            conserved['binding_residues'][res_key] = {
                'conservation': conservation,
                'num_structures': count
            }
    
    return conserved


# Workflow helper functions

def analyze_all_ligands_in_structure(cif_processor, pdb_id: str,
                                   exclude_common: bool = True) -> List[Dict]:
    """
    Complete analysis of all ligands in a structure.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        exclude_common: Exclude water and ions
        
    Returns:
        List of analysis results for each ligand
    """
    results = []
    
    # Extract all ligands
    ligands = extract_all_ligands(cif_processor, pdb_id, exclude_common)
    
    for ligand in ligands:
        # Get interactions
        interactions = calculate_ligand_interactions(
            cif_processor, pdb_id, ligand['atoms']
        )
        
        # Try to get SMILES
        smiles = ligand_to_smiles(ligand['atoms'], ligand['res_name3l'])
        
        result = {
            'ligand_id': ligand['ligand_id'],
            'res_name': ligand['res_name3l'],
            'chain_id': ligand['chain_id'],
            'num_atoms': ligand['num_atoms'],
            'smiles': smiles,
            'interactions': interactions['summary'] if 'summary' in interactions else {},
            'binding_residues': len(interactions.get('binding_residues', []))
        }
        
        results.append(result)
    
    return results


def create_ligand_interaction_report(cif_processor, pdb_id: str,
                                   ligand_name: str, chain_id: Optional[str] = None) -> Dict:
    """
    Create comprehensive interaction report for a ligand.
    
    Args:
        cif_processor: CifBaseProcessor instance
        pdb_id: Structure identifier
        ligand_name: Three-letter ligand code
        chain_id: Optional chain specification
        
    Returns:
        Comprehensive report dictionary
    """
    # Get ligand
    ligand_atoms = get_ligand_by_id(cif_processor, pdb_id, ligand_name, chain_id)
    if ligand_atoms is None:
        return {'error': f'Ligand {ligand_name} not found in {pdb_id}'}
    
    # Get binding site
    binding_site = get_binding_site(cif_processor, pdb_id, ligand_atoms)
    
    # Get interactions
    interactions = calculate_ligand_interactions(cif_processor, pdb_id, ligand_atoms, detailed=True)
    
    # Create report
    report = {
        'pdb_id': pdb_id,
        'ligand': {
            'name': ligand_name,
            'chain': chain_id or 'all',
            'num_atoms': len(ligand_atoms),
            'smiles': ligand_to_smiles(ligand_atoms, ligand_name)
        },
        'binding_site': {
            'num_residues': len(binding_site['residues']),
            'residues': binding_site['residues'].to_dict('records') if not binding_site['residues'].empty else [],
            'volume_estimate': estimate_binding_site_volume(binding_site['atoms']) if not binding_site['atoms'].empty else None
        },
        'interactions': interactions,
        'summary': {
            'total_contacts': sum([
                len(interactions.get('hydrogen_bonds', [])),
                len(interactions.get('hydrophobic', [])),
                len(interactions.get('water_mediated', []))
            ]),
            'key_residues': identify_key_residues(interactions, binding_site['residues']) if not binding_site['residues'].empty else []
        }
    }
    
    return report


def estimate_binding_site_volume(binding_atoms: pd.DataFrame) -> float:
    """Estimate binding site volume using convex hull."""
    try:
        from scipy.spatial import ConvexHull
        coords = binding_atoms[['x', 'y', 'z']].values
        if len(coords) >= 4:  # Need at least 4 points for 3D hull
            hull = ConvexHull(coords)
            return hull.volume
    except:
        pass
    return 0.0


def identify_key_residues(interactions: Dict, binding_residues: pd.DataFrame) -> List[str]:
    """Identify key residues based on interaction patterns."""
    key_residues = []
    
    # Residues with multiple interaction types
    interaction_residues = set()
    
    for hbond in interactions.get('hydrogen_bonds', []):
        if 'acceptor' in hbond:
            interaction_residues.add(hbond['acceptor']['res'])
    
    for hydrophobic in interactions.get('hydrophobic', []):
        interaction_residues.add(hydrophobic['residue'])
    
    # Add residues very close to ligand
    if not binding_residues.empty:
        close_residues = binding_residues[binding_residues['min_distance'] < 3.5]
        for _, res in close_residues.iterrows():
            key_residues.append(f"{res['res_name']}{res['res_id']}")
    
    return list(set(key_residues))
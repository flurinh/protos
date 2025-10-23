"""
Ligand interaction analysis for structure processor.

This module provides methods for analyzing protein-ligand interactions
in structural data, including detection of various interaction types
and binding site characterization.
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix, cKDTree
from typing import Dict, List, Tuple, Optional, Set
import logging

logger = logging.getLogger(__name__)

# Interaction distance cutoffs (in Angstroms)
HBOND_DIST_CUTOFF = 3.5
HYDROPHOBIC_DIST_CUTOFF = 4.5
PI_STACK_DIST_CUTOFF = 5.5
SALT_BRIDGE_CUTOFF = 4.0
WATER_BRIDGE_CUTOFF = 3.5

# Atom type definitions
HBOND_DONORS = {'N', 'O'}  # Simplified - would expand for full implementation
HBOND_ACCEPTORS = {'O', 'N'}
HYDROPHOBIC_ATOMS = {'C'}  # Carbon atoms not in polar groups
AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}
CHARGED_RESIDUES = {
    'positive': {'ARG', 'LYS', 'HIS'},
    'negative': {'ASP', 'GLU'}
}

# Common solvent/ion codes to exclude
COMMON_HETERO = {'HOH', 'WAT', 'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'SO4', 'PO4'}
WATER_CODES = {'HOH', 'WAT'}


class LigandInteractionAnalyzer:
    """Analyzes protein-ligand interactions in structural data."""
    
    def __init__(self, structure_data: pd.DataFrame):
        """
        Initialize analyzer with structure data.
        
        Args:
            structure_data: DataFrame from CifBaseProcessor with structure info
        """
        self.structure = structure_data
        self.protein_atoms = structure_data[structure_data['group'] == 'ATOM'].copy()
        self.hetatoms = structure_data[structure_data['group'] == 'HETATM'].copy()
        
        # Build spatial indices for efficient neighbor searches
        if not self.protein_atoms.empty:
            self.protein_tree = cKDTree(self.protein_atoms[['x', 'y', 'z']].values)
        else:
            self.protein_tree = None
            
    def extract_ligands(
        self,
        exclude_common: bool = True,
        *,
        include_waters: bool = False,
        allowed_res_names: Optional[Set[str]] = None,
    ) -> List[Dict]:
        """
        Extract all ligands from the structure.
        
        Args:
            exclude_common: Whether to exclude water and common ions
            
        Returns:
            List of ligand information dictionaries
        """
        ligands = []
        
        # Group by residue name and chain
        ligand_groups = self.hetatoms.groupby(['res_name3l', 'auth_chain_id', 'auth_seq_id'])

        if exclude_common:
            excluded = COMMON_HETERO.copy()
            if include_waters:
                excluded -= WATER_CODES
        else:
            excluded = set()

        for (res_name, chain_id, res_id), atoms in ligand_groups:
            # Skip common molecules if requested
            if allowed_res_names is not None and res_name not in allowed_res_names:
                continue
            if exclude_common and res_name in excluded:
                continue
                
            # Calculate ligand properties
            coords = atoms[['x', 'y', 'z']].values
            centroid = coords.mean(axis=0)
            
            ligand_info = {
                'res_name3l': res_name,
                'chain_id': chain_id,
                'res_id': res_id,
                'num_atoms': len(atoms),
                'atoms': atoms,
                'centroid': centroid,
                'coords': coords
            }
            
            ligands.append(ligand_info)
            
        return ligands
    
    def get_binding_site_residues(self, ligand_atoms: pd.DataFrame, 
                                 cutoff: float = 5.0) -> pd.DataFrame:
        """
        Get all protein residues within cutoff distance of ligand.
        
        Args:
            ligand_atoms: DataFrame of ligand atoms
            cutoff: Distance cutoff in Angstroms
            
        Returns:
            DataFrame with binding site residue information
        """
        if self.protein_tree is None:
            return pd.DataFrame()
            
        # Find all protein atoms within cutoff of any ligand atom
        ligand_coords = ligand_atoms[['x', 'y', 'z']].values
        
        # Get indices of nearby protein atoms for each ligand atom
        nearby_indices = []
        for coord in ligand_coords:
            indices = self.protein_tree.query_ball_point(coord, cutoff)
            nearby_indices.extend(indices)
            
        # Remove duplicates
        nearby_indices = list(set(nearby_indices))
        
        if not nearby_indices:
            return pd.DataFrame()
            
        # Get the protein atoms
        nearby_atoms = self.protein_atoms.iloc[nearby_indices]
        
        # Group by residue
        residue_info = []
        for (chain, res_id, res_name), res_atoms in nearby_atoms.groupby(
            ['auth_chain_id', 'auth_seq_id', 'res_name3l']
        ):
            # Calculate minimum distance to ligand
            res_coords = res_atoms[['x', 'y', 'z']].values
            distances = distance_matrix(res_coords, ligand_coords)
            min_dist = distances.min()
            
            # Find which atoms are actually contacting
            contact_mask = distances.min(axis=1) <= cutoff
            contact_atoms = res_atoms[contact_mask]['atom_name'].tolist()
            
            residue_info.append({
                'chain_id': chain,
                'res_id': res_id,
                'res_name': res_name,
                'min_distance': min_dist,
                'num_contacts': sum(contact_mask),
                'contact_atoms': contact_atoms
            })
            
        return pd.DataFrame(residue_info)
    
    def detect_hydrogen_bonds(self, ligand_atoms: pd.DataFrame,
                            protein_atoms: Optional[pd.DataFrame] = None) -> List[Dict]:
        """
        Detect potential hydrogen bonds between ligand and protein.
        
        Simplified version - full implementation would include:
        - Proper donor/acceptor identification
        - Angle criteria
        - Protonation state consideration
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms
            
        hbonds = []
        
        # Get potential donors and acceptors
        # This is simplified - real implementation would parse atom types properly
        ligand_donors = ligand_atoms[ligand_atoms['atom_name'].str.contains('N|O', na=False)]
        protein_acceptors = protein_atoms[protein_atoms['atom_name'].str.contains('N|O', na=False)]
        
        if ligand_donors.empty or protein_acceptors.empty:
            return hbonds
            
        # Calculate distances
        donor_coords = ligand_donors[['x', 'y', 'z']].values
        acceptor_coords = protein_acceptors[['x', 'y', 'z']].values
        distances = distance_matrix(donor_coords, acceptor_coords)
        
        # Find pairs within H-bond distance
        donor_idx, acceptor_idx = np.where(distances <= HBOND_DIST_CUTOFF)
        
        for d_idx, a_idx in zip(donor_idx, acceptor_idx):
            donor = ligand_donors.iloc[d_idx]
            acceptor = protein_acceptors.iloc[a_idx]
            
            hbonds.append({
                'donor': {
                    'res': f"{donor['res_name3l']}{donor['auth_seq_id']}",
                    'atom': donor['atom_name'],
                    'chain': donor['auth_chain_id']
                },
                'acceptor': {
                    'res': f"{acceptor['res_name3l']}{acceptor['auth_seq_id']}",
                    'atom': acceptor['atom_name'],
                    'chain': acceptor['auth_chain_id']
                },
                'distance': distances[d_idx, a_idx]
            })
            
        return hbonds
    
    def detect_hydrophobic_contacts(self, ligand_atoms: pd.DataFrame,
                                  protein_atoms: Optional[pd.DataFrame] = None) -> List[Dict]:
        """
        Detect hydrophobic interactions between ligand and protein.
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms
            
        contacts = []
        
        # Get carbon atoms (simplified - would filter out polar carbons)
        ligand_carbons = ligand_atoms[ligand_atoms['atom_name'].str.contains('C', na=False)]
        
        # Get hydrophobic residues
        hydrophobic_res = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
        protein_hydrophobic = protein_atoms[
            protein_atoms['res_name3l'].isin(hydrophobic_res) &
            protein_atoms['atom_name'].str.contains('C', na=False)
        ]
        
        if ligand_carbons.empty or protein_hydrophobic.empty:
            return contacts
            
        # Calculate distances
        ligand_coords = ligand_carbons[['x', 'y', 'z']].values
        protein_coords = protein_hydrophobic[['x', 'y', 'z']].values
        distances = distance_matrix(ligand_coords, protein_coords)
        
        # Find contacts within cutoff
        lig_idx, prot_idx = np.where(distances <= HYDROPHOBIC_DIST_CUTOFF)
        
        # Group by residue to avoid listing every atom pair
        residue_contacts = {}
        for l_idx, p_idx in zip(lig_idx, prot_idx):
            prot_atom = protein_hydrophobic.iloc[p_idx]
            res_key = (prot_atom['auth_chain_id'], prot_atom['auth_seq_id'], 
                      prot_atom['res_name3l'])
            
            if res_key not in residue_contacts:
                residue_contacts[res_key] = {
                    'min_distance': distances[l_idx, p_idx],
                    'num_contacts': 1
                }
            else:
                residue_contacts[res_key]['min_distance'] = min(
                    residue_contacts[res_key]['min_distance'],
                    distances[l_idx, p_idx]
                )
                residue_contacts[res_key]['num_contacts'] += 1
                
        # Convert to list
        for (chain, res_id, res_name), contact_data in residue_contacts.items():
            contacts.append({
                'residue': f"{res_name}{res_id}",
                'chain': chain,
                'distance': contact_data['min_distance'],
                'num_contacts': contact_data['num_contacts']
            })
            
        return contacts
    
    def get_water_mediated_contacts(self, ligand_atoms: pd.DataFrame,
                                   protein_atoms: Optional[pd.DataFrame] = None) -> List[Dict]:
        """
        Find water molecules that bridge ligand and protein.
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms
            
        water_bridges = []
        
        # Get water molecules
        waters = self.hetatoms[self.hetatoms['res_name3l'].isin(['HOH', 'WAT'])]
        
        if waters.empty:
            return water_bridges
            
        # For each water, check if it's close to both ligand and protein
        ligand_coords = ligand_atoms[['x', 'y', 'z']].values
        protein_coords = protein_atoms[['x', 'y', 'z']].values
        
        for _, water in waters.iterrows():
            water_coord = water[['x', 'y', 'z']].values.reshape(1, -1)
            
            # Check distance to ligand
            lig_distances = distance_matrix(water_coord, ligand_coords)
            min_lig_dist = lig_distances.min()
            
            if min_lig_dist > WATER_BRIDGE_CUTOFF:
                continue
                
            # Check distance to protein
            prot_distances = distance_matrix(water_coord, protein_coords)
            min_prot_dist = prot_distances.min()
            
            if min_prot_dist > WATER_BRIDGE_CUTOFF:
                continue
                
            # Find the closest protein residue
            closest_prot_idx = prot_distances.argmin()
            closest_prot = protein_atoms.iloc[closest_prot_idx]
            
            water_bridges.append({
                'water_id': f"{water['res_name3l']}{water['auth_seq_id']}",
                'chain': water['auth_chain_id'],
                'protein_res': f"{closest_prot['res_name3l']}{closest_prot['auth_seq_id']}",
                'protein_chain': closest_prot['auth_chain_id'],
                'distance_to_ligand': min_lig_dist,
                'distance_to_protein': min_prot_dist
            })
            
        return water_bridges
    
    def get_all_interactions(self, ligand_atoms: pd.DataFrame,
                           cutoff: float = 5.0) -> Dict:
        """
        Get comprehensive interaction analysis for a ligand.
        
        Returns:
            Dictionary containing all interaction types and binding residues
        """
        # Get binding site residues first
        binding_residues = self.get_binding_site_residues(ligand_atoms, cutoff)
        
        # Get specific interactions
        interactions = {
            'binding_residues': binding_residues.to_dict('records') if not binding_residues.empty else [],
            'hydrogen_bonds': self.detect_hydrogen_bonds(ligand_atoms),
            'hydrophobic': self.detect_hydrophobic_contacts(ligand_atoms),
            'water_mediated': self.get_water_mediated_contacts(ligand_atoms),
            # Pi-stacking and salt bridges would be implemented similarly
            'pi_stacking': [],
            'salt_bridges': []
        }
        
        # Add summary statistics
        interactions['summary'] = {
            'num_binding_residues': len(binding_residues),
            'num_hydrogen_bonds': len(interactions['hydrogen_bonds']),
            'num_hydrophobic': len(interactions['hydrophobic']),
            'num_water_bridges': len(interactions['water_mediated']),
            'binding_site_chains': binding_residues['chain_id'].unique().tolist() if not binding_residues.empty else []
        }
        
        return interactions


def analyze_ligand_interactions(structure_data: pd.DataFrame, 
                               pdb_id: str,
                               ligand_id: str,
                               chain_id: Optional[str] = None) -> Dict:
    """
    Convenience function to analyze interactions for a specific ligand.
    
    Args:
        structure_data: Full structure DataFrame
        pdb_id: PDB identifier
        ligand_id: Three-letter ligand code
        chain_id: Optional chain specification
        
    Returns:
        Interaction analysis dictionary
    """
    # Filter for specific PDB
    pdb_data = structure_data[structure_data['pdb_id'] == pdb_id]
    
    # Initialize analyzer
    analyzer = LigandInteractionAnalyzer(pdb_data)
    
    # Get ligand atoms
    ligand_filter = (pdb_data['group'] == 'HETATM') & (pdb_data['res_name3l'] == ligand_id)
    if chain_id:
        ligand_filter &= (pdb_data['auth_chain_id'] == chain_id)
        
    ligand_atoms = pdb_data[ligand_filter]
    
    if ligand_atoms.empty:
        logger.warning(f"No ligand {ligand_id} found in {pdb_id}" + 
                      (f" chain {chain_id}" if chain_id else ""))
        return {}
        
    # Get all interactions
    return analyzer.get_all_interactions(ligand_atoms)

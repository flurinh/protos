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
HBOND_ANGLE_CUTOFF = 120.0  # Minimum D-H-A angle in degrees
HYDROPHOBIC_DIST_CUTOFF = 4.5
PI_STACK_DIST_CUTOFF = 5.5
PI_STACK_ANGLE_CUTOFF = 30.0  # Max angle deviation from parallel/perpendicular
SALT_BRIDGE_CUTOFF = 4.0
WATER_BRIDGE_CUTOFF = 3.5
CATION_PI_CUTOFF = 6.0

# Atom type definitions
HBOND_DONORS = {'N', 'O'}  # Simplified - would expand for full implementation
HBOND_ACCEPTORS = {'O', 'N'}
HYDROPHOBIC_ATOMS = {'C'}  # Carbon atoms not in polar groups
AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}

# Aromatic ring atom definitions for each residue type
AROMATIC_RING_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # 6-membered ring
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # 6-membered ring
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],  # 6-membered ring (use larger)
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],  # 5-membered ring
}

# Charged residue definitions with specific charged atoms
CHARGED_RESIDUES = {
    'positive': {'ARG', 'LYS', 'HIS'},
    'negative': {'ASP', 'GLU'}
}

# Atoms carrying charge for salt bridge detection
CHARGED_ATOMS = {
    'ARG': ['NH1', 'NH2', 'NE'],  # Guanidinium group
    'LYS': ['NZ'],  # Terminal amine
    'HIS': ['ND1', 'NE2'],  # Imidazole (protonated)
    'ASP': ['OD1', 'OD2'],  # Carboxylate
    'GLU': ['OE1', 'OE2'],  # Carboxylate
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

    def _get_ring_centroid_and_normal(self, atoms: pd.DataFrame,
                                       ring_atom_names: List[str]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Calculate the centroid and normal vector of an aromatic ring.

        Args:
            atoms: DataFrame of residue atoms
            ring_atom_names: List of atom names that form the ring

        Returns:
            Tuple of (centroid, normal_vector) or None if insufficient atoms
        """
        ring_atoms = atoms[atoms['atom_name'].isin(ring_atom_names)]

        if len(ring_atoms) < 3:
            return None

        coords = ring_atoms[['x', 'y', 'z']].values
        centroid = coords.mean(axis=0)

        # Calculate normal using SVD (best-fit plane)
        centered = coords - centroid
        _, _, vh = np.linalg.svd(centered)
        normal = vh[-1]  # Last row is normal to best-fit plane

        # Normalize
        normal = normal / np.linalg.norm(normal)

        return centroid, normal

    def _find_ligand_aromatic_rings(self, ligand_atoms: pd.DataFrame, max_rings: int = 5) -> List[Dict]:
        """
        Identify aromatic rings in the ligand based on atom connectivity patterns.

        This is a simplified heuristic that looks for planar ring patterns.
        For full accuracy, would need bond information or RDKit.

        Args:
            ligand_atoms: Ligand atom DataFrame
            max_rings: Maximum number of rings to return (prevents combinatorial explosion)
        """
        rings = []

        # Get all carbon and nitrogen atoms (aromatic rings include N in heterocycles)
        aromatic_atoms = ligand_atoms[
            ligand_atoms['atom_name'].str.match(r'^[CN]', na=False)
        ]

        if len(aromatic_atoms) < 5:
            return rings

        coords = aromatic_atoms[['x', 'y', 'z']].values

        # Build simple connectivity based on distance (1.2-1.6 Å for aromatic bonds)
        n_atoms = len(coords)
        dist_matrix = distance_matrix(coords, coords)
        adjacency = (dist_matrix > 1.2) & (dist_matrix < 1.6)

        # Find connected ring-like structures using simple BFS from each atom
        visited_rings = set()

        for ring_size in [6, 5]:  # Try 6-membered first, then 5-membered
            if len(rings) >= max_rings:
                break

            for start_idx in range(n_atoms):
                if len(rings) >= max_rings:
                    break

                # Try to find a ring starting from this atom
                ring_indices = self._find_ring_from_atom(
                    start_idx, ring_size, adjacency, n_atoms
                )

                if ring_indices is None:
                    continue

                # Check if we've already found this ring (same atoms, different order)
                ring_key = tuple(sorted(ring_indices))
                if ring_key in visited_rings:
                    continue

                ring_coords = coords[list(ring_indices)]
                centroid = ring_coords.mean(axis=0)
                centered = ring_coords - centroid

                try:
                    _, s, vh = np.linalg.svd(centered)

                    # If third singular value is small, points are coplanar
                    if len(s) >= 3 and s[2] < 0.5:  # Threshold for planarity
                        normal = vh[-1]
                        normal = normal / np.linalg.norm(normal)

                        # Check if ring-like (atoms roughly equidistant from centroid)
                        distances = np.linalg.norm(centered, axis=1)
                        if np.std(distances) < 0.5:  # Low variance = ring-like
                            visited_rings.add(ring_key)
                            rings.append({
                                'centroid': centroid,
                                'normal': normal,
                                'atom_indices': list(ring_indices),
                                'size': ring_size
                            })
                except np.linalg.LinAlgError:
                    continue

        return rings

    def _find_ring_from_atom(self, start: int, ring_size: int,
                             adjacency: np.ndarray, n_atoms: int) -> Optional[List[int]]:
        """
        Try to find a ring of given size starting from a specific atom using DFS.

        Returns list of atom indices forming the ring, or None if no ring found.
        """
        def dfs(current: int, path: List[int], target_len: int) -> Optional[List[int]]:
            if len(path) == target_len:
                # Check if we can close the ring
                if adjacency[current, start]:
                    return path
                return None

            for neighbor in range(n_atoms):
                if neighbor == start and len(path) < target_len - 1:
                    continue  # Don't return to start too early
                if neighbor in path:
                    continue  # Don't revisit
                if not adjacency[current, neighbor]:
                    continue  # Must be connected

                result = dfs(neighbor, path + [neighbor], target_len)
                if result is not None:
                    return result

            return None

        return dfs(start, [start], ring_size)

    def detect_pi_stacking(self, ligand_atoms: pd.DataFrame,
                          protein_atoms: Optional[pd.DataFrame] = None) -> List[Dict]:
        """
        Detect pi-stacking interactions between aromatic rings.

        Detects both parallel (sandwich) and T-shaped (edge-to-face) stacking.

        Criteria:
        - Centroid distance: <= 5.5 Å
        - Angle between normals: < 30° (parallel) or 60-90° (T-shaped)
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms

        pi_stacks = []

        # Find aromatic rings in ligand
        ligand_rings = self._find_ligand_aromatic_rings(ligand_atoms)

        if not ligand_rings:
            # Fall back: treat ligand as having potential aromatic character
            # if it has planar carbons
            carbons = ligand_atoms[ligand_atoms['atom_name'].str.startswith('C')]
            if len(carbons) >= 5:
                coords = carbons[['x', 'y', 'z']].values
                centroid = coords.mean(axis=0)
                # Estimate normal from SVD
                centered = coords - centroid
                try:
                    _, _, vh = np.linalg.svd(centered)
                    normal = vh[-1]
                    normal = normal / np.linalg.norm(normal)
                    ligand_rings = [{'centroid': centroid, 'normal': normal, 'size': len(carbons)}]
                except:
                    return pi_stacks

        # Find aromatic residues in protein
        aromatic_protein = protein_atoms[protein_atoms['res_name3l'].isin(AROMATIC_RESIDUES)]

        if aromatic_protein.empty:
            return pi_stacks

        # Group by residue
        for (chain, res_id, res_name), res_atoms in aromatic_protein.groupby(
            ['auth_chain_id', 'auth_seq_id', 'res_name3l']
        ):
            ring_atom_names = AROMATIC_RING_ATOMS.get(res_name, [])
            if not ring_atom_names:
                continue

            result = self._get_ring_centroid_and_normal(res_atoms, ring_atom_names)
            if result is None:
                continue

            prot_centroid, prot_normal = result

            # Check against each ligand ring
            for lig_ring in ligand_rings:
                lig_centroid = lig_ring['centroid']
                lig_normal = lig_ring['normal']

                # Calculate centroid distance
                centroid_dist = np.linalg.norm(prot_centroid - lig_centroid)

                if centroid_dist > PI_STACK_DIST_CUTOFF:
                    continue

                # Calculate angle between normals
                dot_product = abs(np.dot(prot_normal, lig_normal))
                angle = np.degrees(np.arccos(np.clip(dot_product, -1, 1)))

                # Classify interaction type
                if angle < PI_STACK_ANGLE_CUTOFF:
                    stack_type = 'parallel'  # Sandwich stacking
                elif angle > (90 - PI_STACK_ANGLE_CUTOFF):
                    stack_type = 'T-shaped'  # Edge-to-face
                else:
                    stack_type = 'offset'  # Offset parallel

                pi_stacks.append({
                    'residue': f"{res_name}{res_id}",
                    'chain': chain,
                    'distance': float(centroid_dist),
                    'angle': float(angle),
                    'type': stack_type,
                    'ligand_ring_size': lig_ring.get('size', 'unknown')
                })

        return pi_stacks

    def detect_salt_bridges(self, ligand_atoms: pd.DataFrame,
                           protein_atoms: Optional[pd.DataFrame] = None) -> List[Dict]:
        """
        Detect salt bridge interactions between charged groups.

        Criteria:
        - Distance between charged atoms: <= 4.0 Å
        - Opposite charges (positive-negative pairing)
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms

        salt_bridges = []

        # Identify charged atoms in ligand
        # Look for carboxylates (O with specific naming) and amines (N)
        ligand_negative = ligand_atoms[
            ligand_atoms['atom_name'].str.contains('O', na=False) &
            ~ligand_atoms['atom_name'].str.contains('OH', na=False)  # Exclude hydroxyls
        ]
        ligand_positive = ligand_atoms[
            ligand_atoms['atom_name'].str.contains('N', na=False)
        ]

        # Get charged residues from protein
        positive_residues = protein_atoms[protein_atoms['res_name3l'].isin(CHARGED_RESIDUES['positive'])]
        negative_residues = protein_atoms[protein_atoms['res_name3l'].isin(CHARGED_RESIDUES['negative'])]

        # Check ligand negative vs protein positive
        if not ligand_negative.empty and not positive_residues.empty:
            for res_name in CHARGED_RESIDUES['positive']:
                res_subset = positive_residues[positive_residues['res_name3l'] == res_name]
                charged_atoms_names = CHARGED_ATOMS.get(res_name, [])

                for (chain, res_id), res_atoms in res_subset.groupby(['auth_chain_id', 'auth_seq_id']):
                    charged_atoms = res_atoms[res_atoms['atom_name'].isin(charged_atoms_names)]

                    if charged_atoms.empty:
                        continue

                    prot_coords = charged_atoms[['x', 'y', 'z']].values
                    lig_coords = ligand_negative[['x', 'y', 'z']].values

                    distances = distance_matrix(prot_coords, lig_coords)
                    min_dist = distances.min()

                    if min_dist <= SALT_BRIDGE_CUTOFF:
                        salt_bridges.append({
                            'protein_residue': f"{res_name}{res_id}",
                            'protein_chain': chain,
                            'protein_charge': 'positive',
                            'ligand_charge': 'negative',
                            'distance': float(min_dist),
                            'protein_atoms': charged_atoms['atom_name'].tolist()
                        })

        # Check ligand positive vs protein negative
        if not ligand_positive.empty and not negative_residues.empty:
            for res_name in CHARGED_RESIDUES['negative']:
                res_subset = negative_residues[negative_residues['res_name3l'] == res_name]
                charged_atoms_names = CHARGED_ATOMS.get(res_name, [])

                for (chain, res_id), res_atoms in res_subset.groupby(['auth_chain_id', 'auth_seq_id']):
                    charged_atoms = res_atoms[res_atoms['atom_name'].isin(charged_atoms_names)]

                    if charged_atoms.empty:
                        continue

                    prot_coords = charged_atoms[['x', 'y', 'z']].values
                    lig_coords = ligand_positive[['x', 'y', 'z']].values

                    distances = distance_matrix(prot_coords, lig_coords)
                    min_dist = distances.min()

                    if min_dist <= SALT_BRIDGE_CUTOFF:
                        salt_bridges.append({
                            'protein_residue': f"{res_name}{res_id}",
                            'protein_chain': chain,
                            'protein_charge': 'negative',
                            'ligand_charge': 'positive',
                            'distance': float(min_dist),
                            'protein_atoms': charged_atoms['atom_name'].tolist()
                        })

        return salt_bridges

    def detect_hydrogen_bonds_with_angles(self, ligand_atoms: pd.DataFrame,
                                          protein_atoms: Optional[pd.DataFrame] = None,
                                          check_angles: bool = True) -> List[Dict]:
        """
        Detect hydrogen bonds with optional angle criteria.

        Full implementation including:
        - Distance criteria (D-A <= 3.5 Å)
        - Angle criteria (D-H-A >= 120°) when hydrogen positions available

        Args:
            ligand_atoms: Ligand atom DataFrame
            protein_atoms: Protein atom DataFrame (optional)
            check_angles: Whether to apply angle criteria

        Returns:
            List of H-bond dictionaries with distance and angle info
        """
        if protein_atoms is None:
            protein_atoms = self.protein_atoms

        hbonds = []

        # Get potential donors (N, O with attached H) and acceptors (N, O)
        ligand_polar = ligand_atoms[ligand_atoms['atom_name'].str.contains('N|O', na=False)]
        protein_polar = protein_atoms[protein_atoms['atom_name'].str.contains('N|O', na=False)]

        if ligand_polar.empty or protein_polar.empty:
            return hbonds

        # Get hydrogen atoms for angle calculation
        ligand_hydrogens = ligand_atoms[ligand_atoms['atom_name'].str.startswith('H')]
        protein_hydrogens = protein_atoms[protein_atoms['atom_name'].str.startswith('H')]

        # Calculate donor-acceptor distances
        lig_coords = ligand_polar[['x', 'y', 'z']].values
        prot_coords = protein_polar[['x', 'y', 'z']].values
        distances = distance_matrix(lig_coords, prot_coords)

        # Find pairs within H-bond distance
        lig_idx, prot_idx = np.where(distances <= HBOND_DIST_CUTOFF)

        for l_idx, p_idx in zip(lig_idx, prot_idx):
            lig_atom = ligand_polar.iloc[l_idx]
            prot_atom = protein_polar.iloc[p_idx]
            dist = distances[l_idx, p_idx]

            # Determine donor/acceptor (simplified: assume donor has H nearby)
            lig_coord = lig_atom[['x', 'y', 'z']].values
            prot_coord = prot_atom[['x', 'y', 'z']].values

            angle = None
            donor_type = 'unknown'

            if check_angles and (not ligand_hydrogens.empty or not protein_hydrogens.empty):
                # Check if ligand atom is donor (has H attached)
                if not ligand_hydrogens.empty:
                    h_coords = ligand_hydrogens[['x', 'y', 'z']].values
                    h_distances = np.linalg.norm(h_coords - lig_coord, axis=1)
                    nearby_h = h_distances < 1.2  # H bonded to heavy atom

                    if nearby_h.any():
                        donor_type = 'ligand'
                        h_idx = np.argmin(h_distances)
                        h_coord = h_coords[h_idx]

                        # Calculate D-H-A angle
                        vec_dh = h_coord - lig_coord
                        vec_ha = prot_coord - h_coord

                        cos_angle = np.dot(vec_dh, vec_ha) / (
                            np.linalg.norm(vec_dh) * np.linalg.norm(vec_ha) + 1e-10
                        )
                        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

                        # Skip if angle too small
                        if angle < HBOND_ANGLE_CUTOFF:
                            continue

                # Check if protein atom is donor
                if donor_type == 'unknown' and not protein_hydrogens.empty:
                    h_coords = protein_hydrogens[['x', 'y', 'z']].values
                    h_distances = np.linalg.norm(h_coords - prot_coord, axis=1)
                    nearby_h = h_distances < 1.2

                    if nearby_h.any():
                        donor_type = 'protein'
                        h_idx = np.argmin(h_distances)
                        h_coord = h_coords[h_idx]

                        vec_dh = h_coord - prot_coord
                        vec_ha = lig_coord - h_coord

                        cos_angle = np.dot(vec_dh, vec_ha) / (
                            np.linalg.norm(vec_dh) * np.linalg.norm(vec_ha) + 1e-10
                        )
                        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

                        if angle < HBOND_ANGLE_CUTOFF:
                            continue

            hbond_info = {
                'ligand_atom': {
                    'res': f"{lig_atom['res_name3l']}{lig_atom['auth_seq_id']}",
                    'atom': lig_atom['atom_name'],
                    'chain': lig_atom['auth_chain_id']
                },
                'protein_atom': {
                    'res': f"{prot_atom['res_name3l']}{prot_atom['auth_seq_id']}",
                    'atom': prot_atom['atom_name'],
                    'chain': prot_atom['auth_chain_id']
                },
                'distance': float(dist),
                'donor': donor_type
            }

            if angle is not None:
                hbond_info['angle'] = float(angle)

            hbonds.append(hbond_info)

        return hbonds

    def get_all_interactions(self, ligand_atoms: pd.DataFrame,
                           cutoff: float = 5.0,
                           use_angle_criteria: bool = True) -> Dict:
        """
        Get comprehensive interaction analysis for a ligand.

        Args:
            ligand_atoms: DataFrame of ligand atoms
            cutoff: Distance cutoff for binding site residues
            use_angle_criteria: Whether to apply angle criteria for H-bonds

        Returns:
            Dictionary containing all interaction types and binding residues
        """
        # Get binding site residues first
        binding_residues = self.get_binding_site_residues(ligand_atoms, cutoff)

        # Get specific interactions - now using all implemented methods
        pi_stacking = self.detect_pi_stacking(ligand_atoms)
        salt_bridges = self.detect_salt_bridges(ligand_atoms)
        hydrogen_bonds = self.detect_hydrogen_bonds_with_angles(
            ligand_atoms, check_angles=use_angle_criteria
        )
        hydrophobic = self.detect_hydrophobic_contacts(ligand_atoms)
        water_mediated = self.get_water_mediated_contacts(ligand_atoms)

        interactions = {
            'binding_residues': binding_residues.to_dict('records') if not binding_residues.empty else [],
            'hydrogen_bonds': hydrogen_bonds,
            'hydrophobic': hydrophobic,
            'water_mediated': water_mediated,
            'pi_stacking': pi_stacking,
            'salt_bridges': salt_bridges
        }

        # Add summary statistics
        interactions['summary'] = {
            'num_binding_residues': len(binding_residues),
            'num_hydrogen_bonds': len(hydrogen_bonds),
            'num_hydrophobic': len(hydrophobic),
            'num_water_bridges': len(water_mediated),
            'num_pi_stacking': len(pi_stacking),
            'num_salt_bridges': len(salt_bridges),
            'binding_site_chains': binding_residues['chain_id'].unique().tolist() if not binding_residues.empty else [],
            'pi_stacking_types': list(set(p['type'] for p in pi_stacking)) if pi_stacking else []
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

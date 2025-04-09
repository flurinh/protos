import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import json
from collections import Counter


def generate_grn_binding_domain_graph(cif_processor, r_max, selected_grns):
    """
    Generate edge indices for binding domains based on the closest atom distances between residues,
    including all residues in loops that contain any GRN from the initial list.
    Also calculates average connectivity of residues.

    Parameters:
    cp_protein_family (object): The protein family object containing structural data.
    r_max (float): Maximum distance for considering residues as neighbors.
    initial_grn_list (list): Initial list of generic residue numbers (GRNs) to consider.

    Returns:
    tuple: A tuple containing:
        - dict: A dictionary where keys are PDB IDs and values are lists of edge indices.
        - set: A set of all unique edges across all structures.
        - dict: A dictionary of average connectivity for each residue.
    """

    def is_loop_grn(grn):
        return 'x' in grn and len(grn.split('x')[0]) == 2

    def get_loop_id(grn):
        if is_loop_grn(grn):
            tm1, tm2 = int(grn[0]), int(grn[1])
            return min(tm1, tm2)
        return None

    graph_dict = {}
    all_edges = set()
    connectivity_count = Counter()
    pdb_count = 0

    for pdb_id in cif_processor.data['pdb_id'].unique():
        # Filter data for the current PDB
        pdb_data = cif_processor.data[cif_processor.data['pdb_id'] == pdb_id]

        # Identify loops to include
        loops_to_include = set(get_loop_id(grn) for grn in selected_grns if is_loop_grn(grn))

        # Create the expanded GRN list
        expanded_grn_list = set(selected_grns)
        for grn in pdb_data['grn'].unique():
            if is_loop_grn(grn) and get_loop_id(grn) in loops_to_include:
                expanded_grn_list.add(grn)

        # Filter for residues in the expanded GRN list
        pdb_data = pdb_data[pdb_data['grn'].isin(expanded_grn_list)]

        if len(pdb_data['grn'].unique()) < 2:
            print(f"Warning: Not enough residues for PDB {pdb_id}. Skipping.")
            continue

        # Group by GRN and calculate pairwise distances between all atoms of different residues
        grns = pdb_data['grn'].unique()
        n_grns = len(grns)
        min_distances = np.full((n_grns, n_grns), np.inf)

        for i, grn1 in enumerate(grns):
            coords1 = pdb_data[pdb_data['grn'] == grn1][['x', 'y', 'z']].values
            for j, grn2 in enumerate(grns[i + 1:], start=i + 1):
                coords2 = pdb_data[pdb_data['grn'] == grn2][['x', 'y', 'z']].values
                distances = cdist(coords1, coords2)
                min_distances[i, j] = min_distances[j, i] = np.min(distances)

        # Find edges based on distance threshold
        edges = np.where((min_distances <= r_max) & (min_distances > 0))

        # Convert node indices to GRNs
        edge_indices = np.array([grns[edges[0]], grns[edges[1]]])

        # Store edge indices for this PDB
        graph_dict[pdb_id] = edge_indices.T.tolist()

        # Add edges to the set of all edges and count connectivity
        for edge in edge_indices.T:
            all_edges.add(tuple(sorted(edge)))
            connectivity_count[edge[0]] += 1
            connectivity_count[edge[1]] += 1

        pdb_count += 1

    # Calculate average connectivity
    avg_connectivity = {grn: count / pdb_count for grn, count in connectivity_count.items()}

    return graph_dict, all_edges, avg_connectivity


def write_binding_pockets_to_json(all_edges_mo, all_edges_ao, all_edges_rbp2,
                                  filename='data/grn/configs/binding_domain.json'):
    """
    Write the binding pockets (edges) for microbial opsins, animal opsins, and RBP2 to a JSON file.

    Parameters:
    all_edges_mo (set): Set of edges for microbial opsins
    all_edges_ao (set): Set of edges for animal opsins
    all_edges_rbp2 (set): Set of edges for RBP2
    filename (str): Path to the output JSON file
    """

    # Convert sets of tuples to lists of lists for JSON serialization
    binding_pockets = {
        'microbial_opsins': [list(edge) for edge in all_edges_mo],
        'gpcr_a': [list(edge) for edge in all_edges_ao],
        'iLBP': [list(edge) for edge in all_edges_rbp2]
    }

    # Write to JSON file
    with open(filename, 'w') as f:
        json.dump(binding_pockets, f, indent=2)

    print(f"Binding pockets have been written to {filename}")


def read_binding_pockets(filename='data/grn/configs/binding_domain.json'):
    """
    Read the binding pockets (edges) for microbial opsins, animal opsins, and RBP2 from a JSON file.

    Parameters:
    filename (str): Path to the input JSON file

    Returns:
    BindingPockets: A named tuple containing:
        - MO: Set of edges for microbial opsins
        - AO: Set of edges for animal opsins
        - RBP2: Set of edges for RBP2
    """

    with open(filename, 'r') as f:
        binding_pockets = json.load(f)

    all_edges_mo = {tuple(edge) for edge in binding_pockets['microbial_opsins']}
    all_edges_ao = {tuple(edge) for edge in binding_pockets['gpcr_a']}
    all_edges_rbp2 = {tuple(edge) for edge in binding_pockets['iLBP']}

    # Return the named tuple
    return {'MO': all_edges_mo, 'AO': all_edges_ao, 'RBP2': all_edges_rbp2}
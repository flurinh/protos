"""
Protos analysis modules.

This package contains analysis functions that work with processor data
without modifying the processors themselves.
"""

from .structure_ligand_analysis import (
    extract_all_ligands,
    extract_water_molecules,
    get_ligand_by_id,
    get_binding_site,
    calculate_ligand_interactions,
    ligand_to_smiles,
    export_ligand_sdf,
    compare_ligand_binding_sites,
    find_conserved_interactions,
    analyze_all_ligands_in_structure,
    create_ligand_interaction_report
)
from .structure_water_networks import (
    WaterNetworkAnalyzer,
    analyze_water_networks,
    summarize_structure_water_network,
    summarize_water_networks,
)
from .structure_embedding_similarity import (
    ChainSelection,
    compute_structure_embedding_similarity,
)

__all__ = [
    'extract_all_ligands',
    'extract_water_molecules',
    'get_ligand_by_id', 
    'get_binding_site',
    'calculate_ligand_interactions',
    'ligand_to_smiles',
    'export_ligand_sdf',
    'compare_ligand_binding_sites',
    'find_conserved_interactions',
    'analyze_all_ligands_in_structure',
    'create_ligand_interaction_report',
    'WaterNetworkAnalyzer',
    'analyze_water_networks',
    'summarize_structure_water_network',
    'summarize_water_networks',
    'ChainSelection',
    'compute_structure_embedding_similarity',
]

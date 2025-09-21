"""
Structure analysis functions.

This module provides analysis functions for protein structures including:
- Geometric calculations (distances, rotations, contacts)
- Structure quality assessment and validation
- Secondary structure and property calculations
- Sequence extraction and analysis
- Membrane protein analysis
- Structure alignment and comparison
- Annotation functions
"""

# Import commonly used functions for convenience
from .geometry import (
    calculate_distance_matrix,
    calculate_rotation_matrix,
    apply_rotation,
    apply_translation,
    flip_structure,
    calculate_distance,
    find_contacts,
    calculate_center_of_mass,
    calculate_radius_of_gyration
)

from .quality import (
    check_missing_atoms,
    validate_bond_lengths,
    check_clashes,
    calculate_b_factor_statistics,
    validate_chirality,
    check_chain_breaks,
    validate_structure_integrity
)

from .properties import (
    calculate_secondary_structure,
    calculate_solvent_accessibility,
    calculate_hydrophobic_moment,
    identify_binding_sites,
    calculate_electrostatic_potential,
    identify_surface_residues
)

from .sequence import (
    extract_sequence,
    extract_all_sequences,
    map_structure_to_sequence,
    identify_missing_residues,
    get_sequence_segments,
    compare_chain_sequences,
    annotate_sequence_conservation
)

from .membrane import (
    calculate_membrane_normal,
    orient_in_membrane,
    orient_n_terminus_up,
    annotate_transmembrane_helices
)

from .alignment import (
    align_structures,
    get_structure_alignment,
    align_on_retinal,
    extract_ca_and_ligand_coords,
    calculate_rmsd,
    kabsch_alignment,
    simple_align_structures
)

from .comparison import (
    compare_all_vs_all,
    compare_one_vs_all,
    normalize_structures
)

__all__ = [
    # Geometry
    'calculate_distance_matrix',
    'calculate_rotation_matrix',
    'apply_rotation',
    'apply_translation',
    'flip_structure',
    'calculate_distance',
    'find_contacts',
    'calculate_center_of_mass',
    'calculate_radius_of_gyration',
    # Quality
    'check_missing_atoms',
    'validate_bond_lengths',
    'check_clashes',
    'calculate_b_factor_statistics',
    'validate_chirality',
    'check_chain_breaks',
    'validate_structure_integrity',
    # Properties
    'calculate_secondary_structure',
    'calculate_solvent_accessibility', 
    'calculate_hydrophobic_moment',
    'identify_binding_sites',
    'calculate_electrostatic_potential',
    'identify_surface_residues',
    # Sequence
    'extract_sequence',
    'extract_all_sequences',
    'map_structure_to_sequence',
    'identify_missing_residues',
    'get_sequence_segments',
    'compare_chain_sequences',
    'annotate_sequence_conservation',
    # Membrane
    'calculate_membrane_normal',
    'orient_in_membrane',
    'orient_n_terminus_up',
    'annotate_transmembrane_helices',
    # Alignment
    'align_structures',
    'get_structure_alignment',
    'align_on_retinal',
    'extract_ca_and_ligand_coords',
    'calculate_rmsd',
    'kabsch_alignment',
    'simple_align_structures',
    # Comparison
    'compare_all_vs_all',
    'compare_one_vs_all',
    'normalize_structures'
]
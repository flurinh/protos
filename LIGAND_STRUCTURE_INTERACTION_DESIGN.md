# Ligand-Structure Interaction Design for Protos

## Overview
This document outlines the comprehensive design for ligand-structure interactions between LigandProcessor and CifBaseProcessor, enabling extraction, conversion, and detailed characterization of protein-ligand interactions.

## Current State Analysis

### Structure Data Format (CifBaseProcessor)
The CifBaseProcessor stores structural data as pandas DataFrames with the following columns:
- `pdb_id`: Structure identifier
- `group`: 'ATOM' (protein) or 'HETATM' (ligands, waters, ions)
- `auth_chain_id`: Chain identifier
- `res_name3l`: Three-letter residue/ligand code (e.g., 'ATP', 'NAD', 'ALA')
- `res_name1l`: One-letter code (for amino acids)
- `auth_seq_id`: Residue sequence number
- `atom_name`: Atom name (e.g., 'CA', 'N', 'O')
- `x`, `y`, `z`: 3D coordinates
- `atom_id`: Unique atom identifier

### Key Observations
1. Ligands are identified by `group == 'HETATM'` and specific `res_name3l` codes
2. Spatial relationships can be calculated using x,y,z coordinates
3. Chain assignment for ligands can be ambiguous (handled by `fix_ligand_chain`)

## Proposed Architecture

### 1. Ligand Extraction and Format Conversion

#### CifBaseProcessor Methods (New)
```python
def extract_ligands(self, pdb_id: str, exclude_common: bool = True) -> List[Dict]:
    """
    Extract all ligands from a structure.
    
    Args:
        pdb_id: Structure identifier
        exclude_common: Exclude water (HOH), common ions (NA, CL, MG, etc.)
    
    Returns:
        List of dicts with ligand info:
        - res_name3l: Three-letter code
        - chain_id: Chain assignment
        - res_id: Residue number
        - atoms: DataFrame of atom coordinates
        - centroid: Geometric center
        - smiles: SMILES string (if convertible)
        - inchi: InChI string (if convertible)
    """

def ligand_to_mol(self, pdb_id: str, ligand_id: str, chain_id: str) -> RDKitMol:
    """
    Convert ligand from structure to RDKit molecule object.
    
    Uses either:
    1. Template matching from PDB Chemical Component Dictionary
    2. Coordinate-based bond inference
    3. SMILES lookup from ligand databases
    """

def export_ligand_sdf(self, pdb_id: str, ligand_id: str, output_path: str):
    """Export ligand in SDF format with 3D coordinates from structure."""
```

#### LigandProcessor Methods (Enhanced)
```python
def import_from_structure(self, cif_processor, pdb_id: str, ligand_id: str) -> str:
    """
    Import ligand from structure and register as entity.
    
    Returns:
        SMILES string (entity ID) for the imported ligand
    """

def get_structure_ligands(self, pdb_id: str) -> Dict[str, List[str]]:
    """
    Get all ligands from a structure mapped to their SMILES.
    
    Returns:
        Dict mapping ligand_id to list of SMILES (multiple for different conformations)
    """
```

### 2. Interaction Characterization

#### CifBaseProcessor Methods (New)
```python
def get_ligand_interactions(self, pdb_id: str, ligand_id: str, 
                           chain_id: str, cutoff: float = 4.0) -> Dict:
    """
    Characterize all interactions for a specific ligand.
    
    Args:
        pdb_id: Structure identifier
        ligand_id: Three-letter ligand code
        chain_id: Chain containing ligand
        cutoff: Distance cutoff in Angstroms
    
    Returns:
        Dict containing:
        - hydrogen_bonds: List of H-bond interactions
        - hydrophobic: List of hydrophobic contacts
        - pi_stacking: List of pi-stacking interactions
        - salt_bridges: List of ionic interactions
        - water_mediated: List of water-mediated contacts
        - binding_residues: List of all interacting residues
    """

def get_binding_site_residues(self, pdb_id: str, ligand_id: str, 
                             chain_id: str, cutoff: float = 5.0) -> pd.DataFrame:
    """
    Get all residues within cutoff distance of ligand.
    
    Returns:
        DataFrame with columns:
        - chain_id, res_id, res_name
        - min_distance: Closest atom-atom distance
        - num_contacts: Number of atom pairs within cutoff
        - contact_atoms: List of contacting atom names
    """

def calculate_interaction_fingerprint(self, pdb_id: str, ligand_id: str) -> np.array:
    """
    Generate interaction fingerprint for ML applications.
    
    Binary vector encoding:
    - Residue types in binding site
    - Interaction types present
    - Geometric descriptors
    """
```

### 3. Advanced Interaction Analysis

#### Interaction Types to Detect

1. **Hydrogen Bonds**
   - Donor-Acceptor distance: 2.5-3.5 Å
   - Angle criteria: >120°
   - Distinguish backbone vs sidechain

2. **Hydrophobic Interactions**
   - C-C distance: 3.5-4.5 Å
   - Consider aromatic vs aliphatic

3. **Pi-Stacking**
   - Aromatic ring centroids: 3.5-5.5 Å
   - Angle between planes: <30° (face-to-face) or 60-90° (edge-to-face)

4. **Salt Bridges/Ionic**
   - Charged group distance: <4.0 Å
   - Consider pKa states

5. **Water-Mediated**
   - Ligand-Water-Protein bridges
   - Conserved water molecules

#### Implementation Details
```python
class InteractionCalculator:
    """Specialized class for interaction detection."""
    
    def __init__(self, structure_data: pd.DataFrame):
        self.structure = structure_data
        self._build_spatial_index()
    
    def _build_spatial_index(self):
        """Build KDTree for efficient neighbor searches."""
        
    def detect_hydrogen_bonds(self, ligand_atoms: pd.DataFrame, 
                             protein_atoms: pd.DataFrame) -> List[HBond]:
        """Detect H-bonds using geometric criteria."""
        
    def detect_pi_stacking(self, ligand_atoms: pd.DataFrame,
                          protein_atoms: pd.DataFrame) -> List[PiStack]:
        """Detect aromatic interactions."""
```

### 4. Binding Site Comparison

```python
def compare_binding_sites(self, ligand_structures: List[Tuple[str, str]]) -> pd.DataFrame:
    """
    Compare binding sites across multiple structures.
    
    Args:
        ligand_structures: List of (pdb_id, ligand_id) tuples
    
    Returns:
        Similarity matrix of binding sites
    """

def get_conserved_interactions(self, ligand_structures: List[Tuple[str, str]]) -> Dict:
    """
    Identify conserved interactions across structures.
    
    Useful for:
    - Pharmacophore modeling
    - Understanding key interactions
    - Drug design
    """
```

### 5. Integration with Other Processors

#### PropertyProcessor Integration
```python
# Store calculated interactions as properties
prop_processor.create_property_table(
    "pdb_ligand_interactions",
    {
        "4HHB_HEM": {
            "num_hbonds": 5,
            "num_hydrophobic": 12,
            "binding_energy_estimate": -8.5,
            "buried_surface_area": 450.2
        }
    }
)
```

#### Workflow Examples

1. **Structure-Based Drug Design Pipeline**
```python
# 1. Load structure
cif_proc = CifBaseProcessor()
cif_proc.load_dataset("kinase_structures")

# 2. Extract and analyze ligands
for pdb_id in cif_proc.pdb_ids:
    ligands = cif_proc.extract_ligands(pdb_id)
    
    for ligand in ligands:
        # Get interactions
        interactions = cif_proc.get_ligand_interactions(
            pdb_id, ligand['res_name3l'], ligand['chain_id']
        )
        
        # Convert to SMILES and register
        smiles = lig_proc.import_from_structure(
            cif_proc, pdb_id, ligand['res_name3l']
        )
        
        # Store interaction data
        prop_proc.assign_property(
            smiles,
            f"interactions_{pdb_id}",
            interactions,
            dataset="structure_interactions"
        )
```

2. **Binding Site Similarity Analysis**
```python
# Find all ATP binding sites
atp_structures = cif_proc.find_structures_with_ligand("ATP")

# Compare binding sites
similarity_matrix = cif_proc.compare_binding_sites(
    [(pdb, "ATP") for pdb in atp_structures]
)

# Find conserved interactions
conserved = cif_proc.get_conserved_interactions(
    [(pdb, "ATP") for pdb in atp_structures]
)
```

### 6. Data Schema Extensions

#### Interaction Data Format
```json
{
    "ligand_id": "ATP",
    "chain_id": "A",
    "interactions": {
        "hydrogen_bonds": [
            {
                "donor": {"res": "SER195", "atom": "OG"},
                "acceptor": {"res": "ATP", "atom": "O2A"},
                "distance": 2.8,
                "angle": 165.2
            }
        ],
        "hydrophobic": [
            {
                "residue": "VAL201",
                "atoms": ["CG1", "CG2"],
                "distance": 3.9
            }
        ],
        "water_mediated": [
            {
                "water_id": "HOH403",
                "protein_res": "ASN197",
                "distance_to_ligand": 2.7,
                "distance_to_protein": 2.9
            }
        ]
    }
}
```

### 7. Visualization Support

```python
def export_interaction_pymol_script(self, pdb_id: str, ligand_id: str, 
                                  output_path: str):
    """Generate PyMOL script to visualize interactions."""

def export_interaction_report(self, pdb_id: str, ligand_id: str,
                            output_path: str, format: str = 'html'):
    """Generate human-readable interaction report with 2D diagrams."""
```

## Additional Capabilities

### 1. Conformational Analysis
- Extract multiple conformations of same ligand
- Compare binding modes across structures
- Identify flexible vs rigid portions

### 2. Allosteric Site Detection
- Identify ligands not in active site
- Characterize allosteric pockets
- Compare to known functional sites

### 3. Fragment-Based Analysis
- Decompose ligands into fragments
- Track fragment binding modes
- Support fragment-based drug design

### 4. Machine Learning Features
- Generate ML-ready features from interactions
- Support for GNN on protein-ligand graphs
- Interaction fingerprints for similarity

### 5. Dynamics Considerations
- B-factor analysis of binding site
- Identify flexible loops near ligand
- Water network stability

## Implementation Priority

1. **Phase 1: Core Extraction**
   - Extract ligands with coordinates
   - Convert to SMILES/SDF
   - Basic distance-based interactions

2. **Phase 2: Detailed Interactions**
   - Implement all interaction types
   - Binding site characterization
   - Integration with properties

3. **Phase 3: Advanced Analysis**
   - Binding site comparison
   - Conformational analysis
   - ML feature generation

4. **Phase 4: Visualization**
   - PyMOL scripts
   - 2D interaction diagrams
   - Web-based viewers

## Dependencies

### Required
- scipy.spatial: Distance calculations
- pandas/numpy: Data manipulation

### Optional but Recommended
- RDKit: Molecule manipulation, SMILES conversion
- BioPython: Additional structure tools
- ProLIF: Protein-ligand interaction fingerprints
- PLIP: Protein-ligand interaction profiler

## Conclusion

This design enables comprehensive analysis of protein-ligand interactions while maintaining Protos' principles of separation of concerns. The CifBaseProcessor handles all spatial/structural calculations, while the LigandProcessor manages molecular representations and the PropertyProcessor stores interaction data for cross-structure analysis.
# Protos Ligand Processor TODO

## Overview
This document outlines the required functionalities to be implemented in the Protos ligand processor module, following Protos' data management principles and architecture.

## Architecture Integration

### Core Requirements
1. **LigandProcessor must inherit from BaseProcessor**
2. **Use ProtosPaths for all file management** - no hardcoded paths
3. **Integrate with EntityRegistry** - track ligands as entities with SMILES as primary identifier
4. **Support DatasetManager** - enable ligand collections and datasets
5. **Follow human-readable naming** - no hash IDs in filenames
6. **Separation of Concerns** - LigandProcessor handles molecular structures only; PropertyProcessor handles ALL metadata including bioactivity data

### Data Organization Structure

Following Protos' separation of concerns principle:

```
data/
├── ligand/                  # LigandProcessor domain
│   ├── sdf/                # SDF/MOL files (molecular structures only)
│   ├── cache/              # Cached fingerprints
│   ├── datasets/           # Ligand dataset definitions
│   └── registry.json       # Ligand processor registry
│
└── property/               # PropertyProcessor domain
    ├── tables/             # Bioactivity data, experimental properties
    │   ├── EGFR_chembl_activities.csv
    │   ├── kinase_inhibitor_properties.csv
    │   └── drug_metabolism_data.csv
    ├── datasets/           # Property dataset definitions
    └── cache/              # Property cache
```

### Architectural Principle: Separation of Concerns

**LigandProcessor** handles:
- Molecular structures (SMILES, SDF, MOL files)
- Structural fingerprints and similarity calculations
- Basic molecular descriptors (MW, LogP as calculated properties)
- Ligand entity registration and tracking

**PropertyProcessor** handles:
- Bioactivity data (IC50, Ki, Kd values)
- Experimental properties (solubility, permeability)
- Target associations
- Any metadata or annotations about ligands
- Cross-entity property relationships

This separation allows:
- Properties to be assigned to ANY entity type (protein, ligand, sequence)
- Consistent property management across all data types
- Flexible property schemas without modifying core processors
- Clean separation between structure and metadata

## Implementation Plan

## 0. Ligand Data Loaders - Foundation Layer

Before implementing the processor functionalities, we need to create loaders that follow Protos' established patterns. These loaders will handle downloading, caching, and registration of ligand data from various sources.

### Loader Architecture Overview

Following the patterns from existing loaders (UniprotDL, GPCRdbDL, download_structures), the ligand loaders will:

1. **Use ProtosPaths exclusively** - No hardcoded paths
2. **Support both individual and bulk operations**
3. **Integrate with EntityRegistry** for tracking
4. **Handle caching and deduplication**
5. **Provide clear error handling and logging**
6. **Support resume/retry for failed downloads**

### ChEMBL Loader Implementation

```python
class ChEMBLDL:
    """
    Downloads ligand data from ChEMBL following Protos loader patterns.
    
    This loader handles:
    - Protein target to ligand mapping
    - Activity data retrieval
    - Compound structure downloads
    - Automatic entity registration
    """
    
    def __init__(self, data_root=None, reload=False, limit=None):
        """
        Initialize ChEMBL loader with ProtosPaths.
        
        Args:
            data_root: Root directory for data. If None, uses default.
            reload: Whether to reload data even if cached.
            limit: Maximum number of compounds to download.
        """
        # Initialize ProtosPaths
        self.paths = ProtosPaths(user_data_root=data_root)
        
        # Get standard directories
        self.ligand_dir = self.paths.get_processor_path('ligand')
        self.sdf_dir = self.paths.get_subdir_path('ligand', 'sdf')
        self.chembl_dir = self.ligand_dir / 'chembl'
        self.cache_dir = self.ligand_dir / 'cache'
        
        # Create directories
        for dir_path in [self.sdf_dir, self.chembl_dir, self.cache_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize entity registry
        self.entity_registry = EntityRegistry(self.paths)
        
        # Settings
        self.reload = reload
        self.limit = limit
        
        # Initialize ChEMBL client
        from chembl_webresource_client.new_client import new_client
        self.chembl = new_client
```

#### Key Methods to Implement:

1. **download_protein_ligands()**
   ```python
   def download_protein_ligands(self, protein_ids, activity_types=['IC50', 'Ki', 'Kd'], 
                               min_pchembl=5.0, save_sdf=True):
       """
       Download ligands for given proteins from ChEMBL.
       
       Args:
           protein_ids: List of protein identifiers (UniProt/PDB)
           activity_types: Types of bioactivities to retrieve
           min_pchembl: Minimum pChEMBL value (negative log of activity)
           save_sdf: Whether to save SDF files
           
       Returns:
           Dict mapping protein_id to list of ligand data
       """
       # Implementation following UniprotDL pattern
   ```

2. **save_ligand_table()**
   ```python
   def save_ligand_table(self, protein_id, ligand_data, force_overwrite=True):
       """
       Save ligand activity data to CSV table.
       
       Table format:
       smiles | chembl_id | activity_type | value | units | pchembl_value
       
       Saves to: ligand/tables/{protein_id}_chembl_activities.csv
       """
   ```

3. **register_ligand()**
   ```python
   def register_ligand(self, smiles, chembl_id, metadata):
       """
       Register ligand in entity system with SMILES as primary ID.
       
       Handles:
       - SMILES sanitization for filenames
       - Alias registration (ChEMBL ID, InChI)
       - Target linkage in metadata
       """
   ```

4. **download_compound_structures()**
   ```python
   def download_compound_structures(self, chembl_ids, format='sdf'):
       """
       Download molecular structures from ChEMBL.
       
       Saves to: ligand/sdf/{sanitized_smiles}.sdf
       """
   ```

### PubChem Loader Implementation

```python
class PubChemDL:
    """
    Downloads ligand data from PubChem as supplementary source.
    
    Features:
    - Bioassay data retrieval
    - Similar compound search
    - Patent information
    - Cross-references to other databases
    """
    
    def __init__(self, data_root=None):
        # Similar initialization to ChEMBLDL
        self.pubchem_dir = self.ligand_dir / 'pubchem'
```

### Enamine Loader Implementation

```python
class EnamineREALDL:
    """
    Interface to Enamine REAL database for commercial compound availability.
    
    Features:
    - Similarity search API
    - Price and availability data
    - Bulk order integration
    """
    
    def __init__(self, data_root=None, api_key=None):
        # API authentication
        # Cache management for search results
```

### Integration with download_with_registration.py

Create a unified download script for ligands:

```python
# src/protos/loaders/download_ligands.py
def download_ligands_with_registration(
    source='chembl',
    identifiers=None,
    data_root=None,
    max_workers=4,
    activity_threshold=100,  # nM
    reload=False
):
    """
    Download ligands from various sources with parallel processing.
    
    Args:
        source: 'chembl', 'pubchem', or 'enamine'
        identifiers: List of protein/compound identifiers
        data_root: Root data directory
        max_workers: Number of parallel download threads
        activity_threshold: Activity cutoff in nM
        reload: Force re-download
        
    Returns:
        Tuple of (successful, failed) downloads
    """
    # Use ThreadPoolExecutor for parallel downloads
    # Register each downloaded ligand
    # Handle rate limiting and retries
```

### Common Utilities

Create shared utilities for all ligand loaders:

```python
# src/protos/loaders/ligand_utils.py

def sanitize_smiles_filename(smiles):
    """
    Convert SMILES to safe filename.
    
    Example:
    CC(=O)OC1=CC=CC=C1C(=O)O → CC_eq_O_OC1_eq_CC_eq_CC_eq_C1C_eq_O_O
    """
    
def validate_smiles(smiles):
    """Validate SMILES string using RDKit."""
    
def smiles_to_inchi(smiles):
    """Convert SMILES to InChI for aliasing."""
    
def extract_protein_mapping(identifier):
    """
    Map various protein identifiers.
    
    Handles: UniProt ID, PDB ID, gene name, ChEMBL target ID
    """
```

### Error Handling and Logging

Following existing patterns:

```python
import logging
logger = logging.getLogger(__name__)

try:
    response = self.chembl.activity.filter(
        target_chembl_id=target_id,
        pchembl_value__gte=min_pchembl
    )
except Exception as e:
    logger.warning(f"Failed to query ChEMBL for {target_id}: {e}")
    failed.append(target_id)
```

### Testing Strategy for Loaders

1. **Mock external APIs**
   ```python
   @patch('chembl_webresource_client.new_client')
   def test_chembl_download(mock_client):
       # Test with mock data
   ```

2. **Test path management**
   ```python
   def test_loader_paths():
       # Verify ProtosPaths integration
       # Check directory creation
   ```

3. **Test entity registration**
   ```python
   def test_ligand_registration():
       # Verify SMILES sanitization
       # Check alias handling
   ```

## 1. ChEMBL Integration - Active Ligands Retrieval (UPDATED)

### Architectural Approach: Loader → LigandProcessor → PropertyProcessor

The ChEMBL integration follows Protos' separation of concerns:

1. **ChEMBLDL (Loader)**: Downloads and caches raw data
2. **LigandProcessor**: Manages molecular structures
3. **PropertyProcessor**: Manages bioactivity data and metadata

### Detailed Implementation Steps:

1. **ChEMBL Loader downloads data**
   ```python
   chembl_loader = ChEMBLDL()
   raw_data = chembl_loader.download_protein_ligands("EGFR")
   # Returns: [{smiles, chembl_id, activity_type, value, units}, ...]
   ```

2. **LigandProcessor handles structures**
   ```python
   lig_proc = LigandProcessor()
   for compound in raw_data:
       # Register ligand entity (structure only)
       lig_proc.save_entity(
           compound['smiles'],
           {
               'smiles': compound['smiles'],
               'chembl_id': compound['chembl_id']  # Basic ID only
           }
       )
   ```

3. **PropertyProcessor handles bioactivity data**
   ```python
   prop_proc = PropertyProcessor()
   
   # Create bioactivity dataset
   for compound in raw_data:
       prop_proc.assign_property(
           entity_name=compound['smiles'],
           property_name=f"{compound['activity_type']}_nM",
           value=compound['value_nm'],
           dataset=f"EGFR_chembl_activities"
       )
   ```

4. **Result: Clean separation**
   - Ligand structures in: `data/ligand/sdf/`
   - Bioactivity data in: `data/property/tables/EGFR_chembl_activities.csv`
   - Both linked through entity registry using SMILES as key

### Benefits of this approach:
- Bioactivity data can be queried across all proteins
- Properties are not locked to ligand format
- Same property system works for protein properties, sequence properties, etc.
- Enables complex cross-entity queries

## 2. Protein Druggability Assessment

### Detailed Implementation Steps:

1. **Integrate with CifBaseProcessor**
   ```python
   def assess_druggability(self, protein_name: str):
       # Use entity registry to find structure
       struct_processor = CifBaseProcessor(paths=self.paths)
       if not struct_processor.entity_exists(protein_name):
           raise ValueError(f"No structure found for {protein_name}")
   ```

2. **Implement pocket detection wrapper**
   - Create abstract `PocketDetector` base class
   - Implement `P2RankDetector` and `CASTpDetector`
   - Store results in `ligand/pockets/{protein_id}_pockets.json`

3. **Pocket characterization pipeline**
   - Calculate volume, surface area, hydrophobicity
   - Store in `ligand/tables/{protein_id}_pocket_properties.csv`
   - Register pockets as sub-entities linked to protein

4. **Create druggability dataset**
   - Aggregate pocket data across proteins
   - Save to `ligand/datasets/druggability_assessment.json`
   - Link to structure entities via registry

## 3. Target Prediction from SMILES

### Detailed Implementation Steps:

1. **Implement ECFP fingerprint generation**
   ```python
   class FingerprintGenerator:
       def __init__(self, paths: ProtosPaths):
           self.cache_dir = paths.get_processor_path("ligand") / "cache" / "fingerprints"
       
       def get_ecfp(self, smiles: str, radius: int = 2):
           # Check cache first
           cache_file = self.cache_dir / f"{sanitize_smiles(smiles)}_ecfp{radius}.pkl"
   ```

2. **Build target prediction system**
   - Pre-compute fingerprints for ChEMBL compounds
   - Store in `ligand/cache/chembl_fingerprints.pkl`
   - Implement fast similarity search with indexed fingerprints

3. **Create prediction results structure**
   - Save to `ligand/similarity/{query_smiles}_targets.json`
   - Include: target_id, similarity_score, known_actives
   - Register as entity linked to query compound

4. **Build target prediction dataset**
   - Track all predictions in `ligand/datasets/target_predictions.json`
   - Enable batch processing and result aggregation

## 4. Enamine REAL Database Search

### Detailed Implementation Steps:

1. **Create Enamine API client**
   - Handle authentication and rate limiting
   - Cache API responses in `ligand/cache/enamine/`

2. **Implement similarity search**
   - Use same fingerprint system as target prediction
   - Store results in `ligand/similarity/{query_smiles}_enamine.json`

3. **Track commercial availability**
   - Create `ligand/tables/enamine_catalog.csv`
   - Columns: smiles, enamine_id, price, availability
   - Link to query compounds via properties

## 5. Generative Chemistry Integration

### Detailed Implementation Steps:

1. **Create abstract generator interface**
   ```python
   class CompoundGenerator(ABC):
       @abstractmethod
       def generate(self, protein: str, pocket: dict, n_compounds: int) -> List[str]:
   ```

2. **Implement Dragonfly wrapper**
   - Use pocket information from druggability assessment
   - Store generated compounds in `ligand/generated/{protein_id}_{timestamp}.json`

3. **Validate and register generated compounds**
   - Check validity with RDKit
   - Register each as entity with metadata indicating origin
   - Create dataset for each generation run

## 6. Binding Affinity Filtering

### Detailed Implementation Steps:

1. **Implement property calculator**
   - MW, LogP, HBA, HBD, TPSA, rotatable bonds
   - Cache in `ligand/cache/properties/{smiles}_props.json`

2. **Create filtering pipeline**
   ```python
   class AffinityFilter:
       def __init__(self, paths: ProtosPaths):
           self.rules = {
               'mw': (150, 500),
               'logp': (-0.4, 5.6),
               'hbd': (0, 5),
               'hba': (0, 10)
           }
   ```

3. **Track filtering results**
   - Store in `ligand/tables/filtered_compounds.csv`
   - Include: smiles, pass/fail, failed_rules
   - Create datasets for passed/failed compounds

## 7. QSAR Model Development

### Detailed Implementation Steps:

1. **Create QSAR pipeline class**
   ```python
   class QSARPipeline:
       def __init__(self, paths: ProtosPaths):
           self.model_dir = paths.get_processor_path("ligand") / "models"
           self.feature_generator = FingerprintGenerator(paths)
   ```

2. **Implement training workflow**
   - Load activity data from ChEMBL results
   - Generate ECFP features (use cached when available)
   - Train XGBoost model with cross-validation
   - Save model to `ligand/models/{target}_{timestamp}.pkl`

3. **Create prediction interface**
   - Load trained models by target
   - Predict on new compounds (single or batch)
   - Store predictions in `ligand/tables/{target}_qsar_predictions.csv`

4. **Track model performance**
   - Save metrics to `ligand/models/{target}_performance.json`
   - Create dataset linking training compounds to predictions

## Cross-Processor Integration (UPDATED)

### Workflows Following Separation of Concerns:

1. **Complete ChEMBL Download → Registration Workflow**
   ```python
   # Step 1: Download ligand data
   chembl_loader = ChEMBLDL()
   compounds = chembl_loader.download_protein_ligands("EGFR", min_pchembl=6.0)
   
   # Step 2: Register structures with LigandProcessor
   lig_proc = LigandProcessor()
   for compound in compounds["EGFR"]:
       lig_proc.save_entity(compound['smiles'], {
           'smiles': compound['smiles'],
           'chembl_id': compound['chembl_id']
       })
   
   # Step 3: Register bioactivity with PropertyProcessor
   prop_proc = PropertyProcessor()
   prop_proc.create_property_dataset(
       dataset_name="EGFR_bioactivity",
       entity_names=[c['smiles'] for c in compounds["EGFR"]],
       properties=['IC50_nM', 'Ki_nM', 'pChEMBL'],
       data=compounds["EGFR"]
   )
   ```

2. **Structure-Based Ligand Discovery**
   ```python
   # Find all ligands for a structure
   struct_proc = CifBaseProcessor()
   structure = struct_proc.load_entity("1M17")
   
   # Get ligands from structure
   lig_proc = LigandProcessor(cif_processor=struct_proc)
   structure_ligands = lig_proc.find_ligand_in_structures("ATP")
   
   # Get bioactivity data for these ligands
   prop_proc = PropertyProcessor()
   for pdb_id in structure_ligands:
       activities = prop_proc.get_entity_properties(f"{pdb_id}_ATP")
   ```

3. **Cross-Entity Property Queries**
   ```python
   # Find all entities with specific property range
   prop_proc = PropertyProcessor()
   
   # Get all molecules with IC50 < 100 nM for any target
   potent_compounds = prop_proc.query_properties(
       property_name="IC50_nM",
       value_range=(0, 100)
   )
   
   # Get their structures
   lig_proc = LigandProcessor()
   structures = [lig_proc.load_entity(smiles) for smiles in potent_compounds]
   
   # Check which ones are drug-like
   drug_like = lig_proc.filter_drug_like(potent_compounds, strict=True)
   ```

4. **Multi-Target Analysis**
   ```python
   # Compare ligands across multiple proteins
   proteins = ["EGFR", "ERBB2", "ERBB3"]
   
   # Download all ligands
   all_compounds = {}
   for protein in proteins:
       compounds = chembl_loader.download_protein_ligands(protein)
       all_compounds[protein] = compounds[protein]
   
   # Find promiscuous ligands
   smiles_to_targets = {}
   for protein, compounds in all_compounds.items():
       for compound in compounds:
           smiles = compound['smiles']
           if smiles not in smiles_to_targets:
               smiles_to_targets[smiles] = []
           smiles_to_targets[smiles].append(protein)
   
   # Register multi-target compounds
   for smiles, targets in smiles_to_targets.items():
       if len(targets) > 1:
           prop_proc.assign_property(
               smiles,
               "target_count",
               len(targets),
               dataset="promiscuous_compounds"
           )
   ```

## Implementation Priority (UPDATED)

### Phase 1: Foundation ✓ COMPLETED
- [x] Basic LigandProcessor with entity management
- [x] ChEMBL loader with proper separation
- [x] Integration with existing processors
- [x] Test scripts demonstrating functionality

### Phase 2: Structure-Ligand Integration ✓ COMPLETED
- [x] Extract ligands from structure files
- [x] Analyze binding pockets and residues
- [x] Calculate detailed interactions (H-bonds, hydrophobic, water bridges)
- [x] Create comprehensive interaction reports

### Phase 3: Ligand Format Conversion (IN PROGRESS)
- [ ] Extract ligand 3D coordinates to separate molecules
- [ ] Convert structure ligands to SMILES using templates
- [ ] Generate SDF files with proper bond orders
- [ ] Register structure ligands in LigandProcessor
- [ ] Link structure ligands to ChEMBL entries

### Phase 4: Interaction Data Management (NEXT)
- [ ] Define interaction data schema and storage format
- [ ] Implement interaction dataset types in PropertyProcessor
- [ ] Create interaction fingerprints for ML
- [ ] Enable cross-structure interaction queries
- [ ] Build interaction comparison tools

### Phase 5: Advanced Structure Features
- [ ] Protein druggability assessment with pocket detection
- [ ] Binding site similarity metrics
- [ ] Conserved interaction patterns
- [ ] Allosteric site identification

### Phase 6: Property Integration Enhancement
- [ ] Update ChEMBL loader for full property pipeline
- [ ] Batch processing for large compound sets
- [ ] Property-based ligand filtering and queries
- [ ] Activity cliff detection

### Phase 7: Machine Learning Features
- [ ] ECFP fingerprint generation and caching
- [ ] Target prediction using similarity
- [ ] QSAR model training pipeline
- [ ] Interaction-based ML features

### Phase 8: External Integrations
- [ ] PDB Chemical Component Dictionary integration
- [ ] PubChem data integration
- [ ] Enamine REAL database search
- [ ] Patent database connections

### Phase 9: Advanced Applications
- [ ] Virtual screening workflows
- [ ] Lead optimization pipelines
- [ ] Generative model integration
- [ ] Fragment-based drug design support

## Interaction Data Storage Design

### Current State
We successfully extract and analyze interactions, but need to formalize storage within Protos' framework.

### Storage Options Analysis

#### Option 1: PropertyProcessor Tables (Recommended)
Store interactions as property tables where:
- Each row is a structure-ligand pair (e.g., "1ATP_ATP")
- Columns represent interaction features
- Enables cross-structure queries

```python
# Example interaction property table
interaction_df = pd.DataFrame({
    'entity_id': ['1ATP_ATP', '4HHB_HEM', '1HVH_CGP'],
    'num_hbonds': [22, 15, 18],
    'num_hydrophobic': [7, 12, 9],
    'num_water_bridges': [4, 2, 3],
    'binding_site_volume': [450.2, 520.1, 380.5],
    'key_residues': ['LYS168,LYS72,GLU127', 'HIS87,VAL68', 'ASP25,GLY27']
})
```

#### Option 2: JSON Datasets in LigandProcessor
Store detailed interactions as JSON with rich structure:
```json
{
  "1ATP_ATP": {
    "structure": "1ATP",
    "ligand": "ATP",
    "interactions": {
      "hydrogen_bonds": [...],
      "hydrophobic": [...],
      "water_mediated": [...]
    }
  }
}
```

#### Option 3: Hybrid Approach (Best of Both)
- Summary statistics in PropertyProcessor (for queries)
- Detailed interactions in LigandProcessor datasets (for visualization)
- Cross-references via entity registry

### Proposed Implementation

1. **Interaction Properties Table**
   - Location: `data/property/tables/structure_interactions.csv`
   - Schema: entity_id, num_hbonds, num_hydrophobic, binding_affinity_estimate
   - Enables: Fast queries, ML features, cross-structure analysis

2. **Detailed Interaction Datasets**
   - Location: `data/ligand/datasets/interaction_details/`
   - Format: JSON files with full interaction geometry
   - Enables: Visualization, detailed analysis, reproducibility

3. **Interaction Fingerprints**
   - Location: `data/ligand/cache/interaction_fingerprints/`
   - Format: Binary numpy arrays
   - Enables: Fast similarity searches, ML training

### Workflow Integration

```python
# After calculating interactions
interactions = calculate_ligand_interactions(struct_proc, pdb_id, ligand_atoms)

# 1. Store summary in PropertyProcessor
prop_proc.assign_property(
    entity_identifier=f"{pdb_id}_{ligand_name}",
    property_name="num_total_interactions",
    property_value=len(interactions['hydrogen_bonds']) + len(interactions['hydrophobic']),
    dataset_name="interaction_summary"
)

# 2. Store details in LigandProcessor dataset
lig_proc.save_interaction_details(
    f"{pdb_id}_{ligand_name}",
    interactions,
    dataset="structure_interactions"
)

# 3. Generate and cache fingerprint
fingerprint = generate_interaction_fingerprint(interactions)
lig_proc.cache_fingerprint(f"{pdb_id}_{ligand_name}", fingerprint)
```

## Dependencies to Add
```toml
[project.optional-dependencies]
ligand = [
    "chembl_webresource_client>=0.10",
    "rdkit>=2023.3.1",
    "scikit-learn>=1.3",
    "xgboost>=2.0",
    "requests>=2.31",
    "p2rank",  # If available as package
    "numpy>=1.24",
    "pandas>=2.0"
]
```

## Testing Strategy

1. **Unit tests for each component**
   - Mock external APIs (ChEMBL, Enamine)
   - Test fingerprint caching
   - Validate entity registration

2. **Integration tests**
   - Full workflow from protein to ligands
   - Cross-processor data flow
   - Dataset creation and loading

3. **Test data organization**
   ```
   tests/test-data/ligand/
   ├── test_compounds.sdf
   ├── test_activities.csv
   ├── test_pockets.json
   └── test_models/
   ```
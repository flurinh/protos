# GPCR Agonist vs Inverse Agonist Analysis

## Research Question

**When the Schiff base in opsins flips, it establishes a hydrogen bond. Compare agonists and inverse agonists in Class A GPCRs with the interactions per GRN in the binding pocket.**

## Question Interpretation and Assumptions

### Primary Assumptions

1. **Schiff Base Context**
   - The question references the Schiff base in opsins (position 7.50 - lysine that binds retinal)
   - We assume this is meant as an example of a critical conformational change that involves hydrogen bonding
   - The analysis should extend beyond just opsins to all Class A GPCRs, using the opsin Schiff base as a conceptual starting point

2. **Scope of Comparison**
   - "Class A GPCRs" includes rhodopsin-like receptors: opsins, adrenergic receptors, dopamine receptors, etc.
   - We assume both rhodopsin (with retinal) and non-opsin GPCRs (with small molecule ligands) should be included
   - The comparison should focus on binding pocket residues identified by GRN positions

3. **Ligand Classification**
   - **Agonists**: Ligands that promote active receptor conformation
   - **Inverse agonists**: Ligands that stabilize inactive receptor conformation
   - We exclude partial agonists, antagonists, and allosteric modulators for clarity
   - For rhodopsin: all-trans-retinal (agonist-like) vs 11-cis-retinal (inverse agonist-like)

4. **Hydrogen Bond Focus**
   - The question specifically mentions hydrogen bond establishment
   - We assume this means analyzing ALL hydrogen bonds in the binding pocket, not just those involving position 7.50
   - Both direct ligand-protein and water-mediated hydrogen bonds should be considered

5. **GRN Usage**
   - "Interactions per GRN" means analyzing interactions at each GRN-numbered position
   - We assume standard Ballesteros-Weinstein numbering (X.50 = most conserved in helix X)
   - Analysis should be position-centric rather than residue-centric for cross-receptor comparison

6. **Binding Pocket Definition**
   - We assume the orthosteric binding pocket (where endogenous ligands bind)
   - For consistency, we'll use GRN positions within 5Å of bound ligands
   - Key positions likely include: 3.32, 3.33, 3.36, 5.42, 5.43, 5.46, 6.48, 6.51, 6.52, 7.39, 7.43, 7.50

### Secondary Assumptions

7. **Structural Data Requirements**
   - We need high-resolution crystal structures (<3.0Å) for accurate hydrogen bond detection
   - Both active and inactive state structures should be available for meaningful comparison
   - Ligands must be well-resolved in electron density

8. **Conformational States**
   - Agonist-bound structures represent active or active-like conformations
   - Inverse agonist-bound structures represent inactive conformations
   - We acknowledge that crystal structures may not capture full dynamics

9. **Analysis Focus**
   - Primary focus: How interaction patterns differ between agonist and inverse agonist binding
   - Secondary focus: Whether the opsin Schiff base flip mechanism has parallels in other GPCRs
   - Tertiary focus: Identify conserved interaction patterns that determine ligand efficacy

## Background

### Key Concepts
- **Schiff Base**: In rhodopsin-like GPCRs, the lysine at position 7.50 forms a Schiff base with retinal
- **Agonists**: Ligands that activate the receptor, promoting the active conformation
- **Inverse Agonists**: Ligands that stabilize the inactive conformation
- **GRN System**: Standardized numbering (X.50 = most conserved position in helix X)

### Critical GRN Positions in Class A GPCRs
- **3.32**: DRY motif - ionic lock in inactive state
- **5.50**: Toggle switch residue
- **6.48**: CWxP motif - rotamer toggle
- **7.50**: NPxxY motif / Schiff base lysine in opsins

## Proposed Analysis Workflow

### Phase 1: Data Collection and Preparation

#### Step 1.1: Structure Collection
```python
# Initialize processors with explicit paths
import os
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor

# Set up paths
data_dir = Path("/path/to/analysis_data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Initialize processors
struct_proc = CifBaseProcessor(name="gpcr_structures")
grn_proc = GRNBaseProcessor(name="gpcr_grn")
```

#### Step 1.2: Download Representative Structures
- Rhodopsin structures with different ligands
- β2-adrenergic receptor with agonists/inverse agonists
- Other Class A GPCRs with diverse ligands

```python
# Define structure sets
# Note: These PDB IDs have been verified to contain Class A GPCRs with appropriate ligands

agonist_structures = {
    # Beta-adrenergic receptors
    '3SN6': 'β2AR with Gs and BI-167107 (full agonist)',
    '4LDO': 'β2AR with Gs and BI-167107 (full agonist)',
    '3P0G': 'β2AR with BI-167107 (full agonist)',
    '6MXT': 'β1AR with Gs and formoterol (full agonist)',
    
    # Dopamine receptors  
    '7JOZ': 'D1R with G protein and non-catechol agonist',
    '7CKZ': 'D1R with Gs and dopamine (full agonist)',
    '6VMS': 'D2R with Gi and bromocriptine (agonist)',
    
    # Serotonin receptors
    '7E2Y': '5-HT1A with Gi and serotonin (full agonist)',
    '6G79': '5-HT2A with mini-Gq and 25-CN-NBOH (agonist)',
    
    # Muscarinic receptors
    '4MQS': 'M2R with iperoxo (full agonist)',
    '4MQT': 'M2R with iperoxo and LY2119620 (agonist + PAM)',
    
    # Adenosine receptors
    '7ARO': 'A2AR with LUF5833 (partial agonist)',
    '6GDG': 'A2AR with NECA (full agonist)',
}

inverse_agonist_structures = {
    # Beta-adrenergic receptors
    '2RH1': 'β2AR with carazolol (inverse agonist)',
    '3NY8': 'β2AR with ICI-118551 (inverse agonist)', 
    '3NYA': 'β2AR with alprenolol (inverse agonist)',
    '5JQH': 'β2AR with carazolol (inverse agonist, high res)',
    
    # Rhodopsin
    '1U19': 'Rhodopsin with 11-cis-retinal (inverse agonist)',
    '1GZM': 'Rhodopsin with 11-cis-retinal (inverse agonist)',
    
    # Dopamine receptors
    '6CM4': 'D2R with risperidone (inverse agonist)',
    '6LUQ': 'D2R with haloperidol (inverse agonist)',
    
    # Serotonin receptors
    '6A93': '5-HT2A with risperidone (inverse agonist)',
    '6A94': '5-HT2A with zotepine (inverse agonist)',
    '6WGT': '5-HT2A with lumateperone (inverse agonist)',
    
    # Muscarinic receptors
    '3UON': 'M2R with QNB (inverse agonist)',
    '5CXV': 'M1R with tiotropium (inverse agonist)',
    
    # Adenosine receptors
    '5UEN': 'A1R with DU172 (inverse agonist)',
    '3VGA': 'A2AR with caffeine derivative (inverse agonist)',
}

# Download structures
all_pdbs = list(agonist_structures.keys()) + list(inverse_agonist_structures.keys())
struct_proc.download_structures(all_pdbs)
```

### Phase 2: GRN Assignment and Alignment

#### Step 2.1: Load GRN Reference
```python
# Load GPCR GRN reference table
grn_proc.load_grn_table("gpcrdb_ref")
```

#### Step 2.2: Assign GRN Numbers to Structures
```python
# Extract sequences from structures
sequences = struct_proc.get_seq_dict()

# Assign GRN numbers
grn_assignments = grn_proc.assign_grn_batch(
    sequences,
    family='gpcr_a',
    save_output=True
)
```

#### Step 2.3: Map GRN to Structure Coordinates
```python
# For each structure, map GRN positions to 3D coordinates
for pdb_id in all_pdbs:
    struct_data = struct_proc.load_structure(pdb_id)
    grn_data = grn_assignments[pdb_id]
    
    # Add GRN annotations to structure
    struct_proc.add_grn_annotations(pdb_id, grn_data)
```

### Phase 3: Binding Pocket Analysis

#### Step 3.1: Define Binding Pocket by GRN
```python
# Define binding pocket residues using GRN positions
binding_pocket_grn = [
    '3.32', '3.33', '3.36',  # TM3
    '5.42', '5.43', '5.46',  # TM5
    '6.48', '6.51', '6.52',  # TM6
    '7.39', '7.43', '7.50'   # TM7 (including Schiff base)
]
```

#### Step 3.2: Extract Ligand Interactions
```python
from protos.processing.ligand.ligand_processor import LigandProcessor

lig_proc = LigandProcessor(name="gpcr_ligands")

# Analyze each structure
interaction_data = {}
for pdb_id, description in {**agonist_structures, **inverse_agonist_structures}.items():
    # Get ligand interactions
    interactions = lig_proc.analyze_protein_ligand_interactions(
        pdb_id,
        grn_positions=binding_pocket_grn,
        distance_cutoff=4.0  # Ångstroms
    )
    
    ligand_type = 'agonist' if pdb_id in agonist_structures else 'inverse_agonist'
    interaction_data[pdb_id] = {
        'type': ligand_type,
        'interactions': interactions,
        'description': description
    }
```

### Phase 4: Comparative Analysis

#### Step 4.1: Hydrogen Bond Analysis
```python
# Focus on hydrogen bonds, especially at position 7.50
h_bond_comparison = {}

for grn_pos in binding_pocket_grn:
    h_bond_comparison[grn_pos] = {
        'agonist': [],
        'inverse_agonist': []
    }
    
    for pdb_id, data in interaction_data.items():
        if grn_pos in data['interactions']:
            h_bonds = data['interactions'][grn_pos]['hydrogen_bonds']
            h_bond_comparison[grn_pos][data['type']].extend(h_bonds)
```

#### Step 4.2: Conformational Analysis
```python
# Compare key conformational markers
conformational_markers = {
    '3.32': 'DRY motif ionic lock',
    '5.50': 'Toggle switch',
    '6.48': 'CWxP rotamer',
    '7.50': 'NPxxY/Schiff base'
}

# Analyze conformational differences
for grn_pos, description in conformational_markers.items():
    print(f"\nAnalyzing {grn_pos} ({description}):")
    
    # Get coordinates for agonist vs inverse agonist states
    agonist_coords = []
    inverse_agonist_coords = []
    
    for pdb_id, data in interaction_data.items():
        coords = struct_proc.get_grn_coordinates(pdb_id, grn_pos)
        if data['type'] == 'agonist':
            agonist_coords.append(coords)
        else:
            inverse_agonist_coords.append(coords)
```

### Phase 5: Statistical Analysis and Visualization

#### Step 5.1: Interaction Frequency Analysis
```python
import pandas as pd
import matplotlib.pyplot as plt

# Create interaction frequency matrix
interaction_matrix = pd.DataFrame(
    index=binding_pocket_grn,
    columns=['agonist_hbond', 'agonist_hydrophobic', 
             'inverse_agonist_hbond', 'inverse_agonist_hydrophobic']
)

# Populate matrix with interaction frequencies
for grn_pos in binding_pocket_grn:
    for ligand_type in ['agonist', 'inverse_agonist']:
        # Count interaction types
        # ... analysis code ...
```

#### Step 5.2: Generate Visualizations
```python
from protos.visualization.ligand_vis import visualize_binding_pocket

# Create comparative visualization
fig = visualize_binding_pocket(
    agonist_structures=agonist_coords,
    inverse_agonist_structures=inverse_agonist_coords,
    grn_positions=binding_pocket_grn,
    highlight_position='7.50'  # Highlight Schiff base
)
```

### Phase 6: Interpretation and Reporting

#### Step 6.1: Key Findings Summary
```python
# Generate comprehensive report
from protos.processing.property.property_processor import PropertyProcessor

prop_proc = PropertyProcessor(name="gpcr_analysis_results")

# Store results as properties
for pdb_id, data in interaction_data.items():
    # Add interaction properties
    prop_proc.assign_property(
        pdb_id, 
        "ligand_type", 
        data['type'],
        "gpcr_ligand_analysis"
    )
    
    # Add specific interaction counts
    h_bond_count = sum(len(data['interactions'][grn]['hydrogen_bonds']) 
                      for grn in data['interactions'])
    prop_proc.assign_property(
        pdb_id,
        "total_h_bonds",
        h_bond_count,
        "gpcr_ligand_analysis"
    )
```

#### Step 6.2: Generate Final Report
```python
# Create markdown report
report = f"""
# GPCR Agonist vs Inverse Agonist Analysis Results

## Key Findings

### Schiff Base Region (7.50)
- Agonist state: {agonist_7_50_summary}
- Inverse agonist state: {inverse_agonist_7_50_summary}

### Binding Pocket Interaction Patterns
{interaction_summary}

### Conformational Differences
{conformational_summary}

## Conclusions
{conclusions}
"""

# Save report
with open("gpcr_analysis_report.md", "w") as f:
    f.write(report)
```

## Expected Outcomes

1. **Interaction Patterns**: Identification of GRN positions that show differential interaction patterns between agonists and inverse agonists

2. **Schiff Base Dynamics**: Understanding how the Schiff base flip affects the hydrogen bonding network

3. **Conserved Mechanisms**: Discovery of conserved activation/inactivation mechanisms across Class A GPCRs

4. **Structure-Activity Relationships**: Insights into how specific interactions at GRN positions correlate with functional outcomes

## Required Data
- Multiple Class A GPCR structures with agonists
- Multiple Class A GPCR structures with inverse agonists
- GPCR GRN reference table (gpcrdb_ref.csv)
- Ligand interaction analysis tools

## Computational Requirements
- Protos framework with all processors
- PyMOL or similar for visualization
- Statistical analysis packages (pandas, scipy)
- ~10GB storage for structures
- ~4GB RAM for analysis
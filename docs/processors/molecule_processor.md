# MoleculeProcessor

Manages small molecule descriptors (SMILES, InChI) and related metadata.

**Location:** `protos.processing.molecule.molecule_processor`

**Processor Type:** `molecule`

## Overview

The MoleculeProcessor provides:
- Storage of molecule records (SMILES, InChI)
- Molecular property calculation (requires RDKit)
- Drug-likeness filtering (Lipinski's Rule of Five)
- Integration with structure entities via relationships

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(name)` | Load molecule record by name. Returns dict with SMILES, metadata. |
| `save_entity(name, data, metadata)` | Save molecule from SMILES string or record dict. |
| `list_ligands()` | List all registered molecule entities. |

### Batch Operations

| Method | Description |
|--------|-------------|
| `register_smiles_map(smiles_map, dataset_name, metadata)` | Register multiple molecules from SMILES dict. |

### Property Calculation

| Method | Description |
|--------|-------------|
| `calculate_properties(smiles)` | Calculate molecular properties (MW, LogP, etc.). Requires RDKit. |
| `filter_drug_like(entity_names, strict)` | Filter molecules by Lipinski's Rule of Five. |

### Path Properties

| Property | Description |
|----------|-------------|
| `records_dir` | JSON molecule records directory |
| `datasets_dir` | Dataset definitions directory |

---

## Data Format

Molecules are stored as JSON records:

```json
{
    "name": "aspirin",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "inchi": "InChI=1S/C9H8O4/...",
    "kind": "smiles",
    "metadata": {
        "source": "pubchem",
        "cid": 2244
    }
}
```

### Record Fields

| Field | Type | Description |
|-------|------|-------------|
| `name` | str | Molecule identifier |
| `smiles` | str | SMILES representation |
| `inchi` | str | InChI representation (optional) |
| `kind` | str | Record type: "smiles", "inchi", "structure_record" |
| `metadata` | dict | Additional information |

---

## Molecular Properties

Properties calculated by `calculate_properties()` (requires RDKit):

| Property | Description |
|----------|-------------|
| `molecular_weight` | Molecular weight (Da) |
| `logp` | Octanol-water partition coefficient |
| `hbd` | Hydrogen bond donors |
| `hba` | Hydrogen bond acceptors |
| `tpsa` | Topological polar surface area |
| `rotatable_bonds` | Number of rotatable bonds |
| `num_rings` | Number of rings |
| `num_aromatic_rings` | Number of aromatic rings |

---

## Usage Examples

### Basic Operations

```python
from protos.processing.molecule import MoleculeProcessor

proc = MoleculeProcessor()

# Save from SMILES string
proc.save_entity("aspirin", "CC(=O)Oc1ccccc1C(=O)O")

# Save with full record
proc.save_entity("caffeine", {
    "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "kind": "smiles",
    "metadata": {"common_name": "caffeine", "cas": "58-08-2"}
})

# Load molecule
record = proc.load_entity("aspirin")
print(f"SMILES: {record['smiles']}")
print(f"Metadata: {record.get('metadata', {})}")

# List all molecules
molecules = proc.list_ligands()
```

### Batch Registration

```python
# Register multiple molecules
smiles_map = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "acetaminophen": "CC(=O)Nc1ccc(O)cc1"
}

dataset_id = proc.register_smiles_map(
    smiles_map,
    dataset_name="common_drugs",
    metadata={"source": "manual_curation"}
)

print(f"Created dataset: {dataset_id}")
```

### Property Calculation

```python
# Calculate properties (requires RDKit)
props = proc.calculate_properties("CC(=O)Oc1ccccc1C(=O)O")

if props:
    print(f"MW: {props['molecular_weight']:.2f} Da")
    print(f"LogP: {props['logp']:.2f}")
    print(f"HBD: {props['hbd']}")
    print(f"HBA: {props['hba']}")
    print(f"TPSA: {props['tpsa']:.2f}")
    print(f"Rotatable bonds: {props['rotatable_bonds']}")
```

### Drug-Likeness Filtering

```python
molecules = ["aspirin", "caffeine", "ibuprofen", "large_peptide"]

# Standard Lipinski filtering
drug_like = proc.filter_drug_like(molecules)
print(f"Drug-like: {drug_like}")

# Strict filtering
strict_drug_like = proc.filter_drug_like(molecules, strict=True)
print(f"Strict drug-like: {strict_drug_like}")
```

### Integration with LigandLoader

```python
from protos import LigandLoader
from protos.processing.molecule import MoleculeProcessor

# Import ligands (creates both structure and molecule entities)
loader = LigandLoader()
dataset_id, names = loader.import_smiles({
    "ligand_1": "CC(=O)Oc1ccccc1C(=O)O",
    "ligand_2": "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
}, dataset_name="my_ligands", generate_3d=True)

# Access molecule records
mol_proc = MoleculeProcessor()
for name in names:
    record = mol_proc.load_entity(name)
    print(f"{name}: {record['smiles']}")
```

### Finding Related Structures

```python
from protos.io.core import get_registry

registry = get_registry()

# Molecules from LigandLoader have structure relationships
related = registry.get_related_entities(
    "ligand_1",
    rel_type="has_structure"
)

for struct_name in related:
    print(f"Structure entity: {struct_name}")
```

### Analysis Pipeline

```python
proc = MoleculeProcessor()

# Register compound library
compounds = {
    "cmpd_001": "CC(=O)Oc1ccccc1C(=O)O",
    "cmpd_002": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "cmpd_003": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "cmpd_004": "CC(=O)Nc1ccc(O)cc1",
    "cmpd_005": "c1ccccc1"  # benzene
}

proc.register_smiles_map(compounds, dataset_name="screening_library")

# Analyze all compounds
results = []
for name in compounds.keys():
    record = proc.load_entity(name)
    props = proc.calculate_properties(record["smiles"])

    if props:
        results.append({
            "name": name,
            "smiles": record["smiles"],
            **props
        })

# Filter by properties
import pandas as pd
df = pd.DataFrame(results)
drug_like_df = df[(df["molecular_weight"] < 500) & (df["logp"] < 5)]
print(f"Drug-like compounds: {len(drug_like_df)}/{len(df)}")
```

### Dataset Operations

```python
# Create dataset manually
proc.create_dataset(
    "fragment_library",
    ["frag_001", "frag_002", "frag_003"]
)

# Access via dataset manager
entities = proc.dataset_manager.get_dataset_entities("fragment_library")

# Delete dataset
proc.dataset_manager.delete_dataset("old_library")
```

### Name Sanitization

Entity names are automatically sanitized:

```python
# Spaces and special characters handled
proc.save_entity("My Compound (v2)", "CCO")

# Stored as "My_Compound_v2"
record = proc.load_entity("My_Compound_v2")
```

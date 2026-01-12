# StructureProcessor

Manages protein structures as pandas DataFrames with PKL canonical storage.

**Location:** `protos.processing.structure.structure_processor`

**Processor Type:** `structure`

## Overview

The StructureProcessor provides:
- Loading structures from CIF files and PKL cache
- DataFrame representation with MultiIndex `(structure_id, atom_id)`
- Automatic entity registration
- Ligand and contact analysis
- Structure alignment and transformation
- GRN annotation
- Structure export to CIF/PDB formats

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(structure_id)` | Load structure from registry or PKL cache. Returns DataFrame or None. |
| `save_entity(structure_id, df, metadata)` | Save structure as PKL and register with EntityRegistry. |
| `delete_entity(name)` | Remove structure from registry and delete PKL file. |
| `export_entity(name, out_path, format)` | Export structure to CIF or PDB file. |
| `export_dataset(dataset_name, output_dir, format)` | Export all structures in a dataset. |

### Dataset Operations

| Method | Description |
|--------|-------------|
| `load_dataset(dataset_name, return_format)` | Load all structures in dataset. Returns stacked DataFrame or dict. |
| `save_dataset(dataset_name, entities, metadata)` | Create dataset from structure entities. |
| `get_dataset(dataset_name)` | Get dataset configuration. |
| `get_dataset_entities(dataset_name)` | Get list of entity names in dataset. |
| `load_structures(identifiers)` | Load multiple structures by name. |

### Chain and Sequence Operations

| Method | Description |
|--------|-------------|
| `get_chains(structure_id)` | Get list of chain IDs in structure. |
| `get_sequence(structure_id, chain_id)` | Extract amino acid sequence from chain. |
| `get_all_sequences(structure_id)` | Get sequences for all chains as dict. |
| `collect_chain_sequences(structure_id, chains)` | Collect sequences from specified chains. |
| `register_chain_sequences(structure_id, chains, dataset_name)` | Register chain sequences as separate entities. |
| `filter_by_chain(structure_id, chain_ids, invert)` | Filter structure to specific chains. |

### Structure Filtering

| Method | Description |
|--------|-------------|
| `filter_structure(structure_id, chains, residue_range, atom_names, remove_hetatm)` | Apply multiple filters to structure. |
| `filter_by_residue_range(structure_id, chain_id, start, end)` | Filter to residue number range. |
| `remove_hetatm(structure_id)` | Remove all HETATM records (water, ligands, ions). |

### Structure Modification

| Method | Description |
|--------|-------------|
| `delete_atoms(structure_id, atom_ids)` | Remove specific atoms by ID. |
| `delete_residues(structure_id, chain_id, residue_ids)` | Remove specific residues. |
| `delete_chain(structure_id, chain_id)` | Remove entire chain from structure. |
| `add_atoms(structure_id, atoms_df)` | Add new atoms to structure. |
| `add_ligand(structure_id, ligand_df, chain_id)` | Add ligand to structure. |
| `reindex_atom_ids(structure_id)` | Renumber atom IDs sequentially. |
| `renumber_residues(structure_id, chain_id, start)` | Renumber residues starting from value. |

### Annotation

| Method | Description |
|--------|-------------|
| `annotate_structure(structure_id, column, values, selector)` | Add custom annotation column. |
| `annotate_chain(structure_id, chain_id, column, values)` | Annotate specific chain. |
| `annotate_with_grn(structure_id, grn_table, sequence_id)` | Add GRN positions from table. |
| `annotate_structures_with_grn(structure_ids, grn_table)` | Batch GRN annotation. |
| `assign_grns(structure_id, grn_mapping)` | Assign GRN from position mapping. |

### Ligand Analysis

| Method | Description |
|--------|-------------|
| `list_ligands(structure_id, exclude_common, min_atoms)` | List all ligands in structure. |
| `summarize_ligands(structure_id, group_by, include_chains, min_atoms)` | Get detailed ligand summary. |
| `get_ligand_interactions(structure_id, ligand_id, cutoff, include_water_bridges)` | Find protein-ligand contacts. |

### Contact Analysis

| Method | Description |
|--------|-------------|
| `get_water_contacts(structure_id, cutoff, protein_only)` | Find water-mediated contacts. |
| `get_ion_contacts(structure_id, ion_codes, cutoff)` | Find ion coordination sites. |
| `compute_water_networks(structure_id, distance_cutoff)` | Analyze water hydrogen bond networks. |

### Structure Alignment

| Method | Description |
|--------|-------------|
| `align_structures(structure_ids, reference, method)` | Align multiple structures to reference. |
| `align_and_record(query_id, reference_id, method)` | Align pair and store result. |
| `align_pair(structure_id_1, structure_id_2)` | Align two structures, return RMSD. |
| `align_all_to_reference(structure_ids, reference_id)` | Align all to single reference. |
| `align_one_vs_all(query_id, target_ids)` | Align one against multiple targets. |
| `align_all_vs_all(structure_ids)` | Pairwise alignment of all structures. |
| `export_aligned_structures(alignment_results, output_dir)` | Export aligned coordinates. |

### Coordinate Transformation

| Method | Description |
|--------|-------------|
| `apply_transformation(structure_id, rotation, translation)` | Apply rotation/translation matrix. |
| `orient_structure(structure_id, method)` | Reorient structure (e.g., membrane alignment). |
| `get_ca_coordinates(structure_id, chain_id)` | Extract CA atom coordinates as numpy array. |

### Utility Methods

| Method | Description |
|--------|-------------|
| `compute_content_hash(df)` | Generate SHA256 hash of core structural content. |
| `structure_ids` | Property: list of all loaded structure IDs. |
| `pdb_ids` | Property: alias for structure_ids. |
| `alignment_engine` | Property: lazy-loaded StructureAlignmentEngine. |
| `list_related_sequences(structure_ids)` | Find sequence entities related to structures. |
| `list_dataset_related_sequences(dataset_name)` | Find sequences related to dataset structures. |

### Path Properties

| Property | Description |
|----------|-------------|
| `path_pkl_dir` | PKL cache directory |
| `path_cif_dir` | CIF structure files |
| `path_pdb_dir` | PDB export directory |
| `path_sdf_dir` | SDF ligand files |
| `path_dataset_dir` | Dataset PKL files |
| `path_temp_dir` | Temporary files |

---

## Data Format

Structures are stored as DataFrames with MultiIndex `(structure_id, atom_id)`:

| Column | Type | Description |
|--------|------|-------------|
| `structure_id` | str | Entity identifier |
| `atom_id` | Int64 | Unique atom ID |
| `auth_chain_id` | str | Chain identifier |
| `auth_seq_id` | Int64 | Residue number |
| `res_name3l` | str | 3-letter residue name |
| `atom_name` | str | Atom name (CA, CB, etc.) |
| `element` | str | Element symbol |
| `x`, `y`, `z` | float | Coordinates (Angstroms) |
| `b_factor` | float | B-factor / temperature factor |
| `occupancy` | float | Occupancy (0-1) |
| `group` | str | ATOM or HETATM |
| `grn` | str | GRN position (if annotated) |

---

## Usage Examples

### Basic Operations

```python
from protos import StructureProcessor

proc = StructureProcessor()

# Load and inspect
df = proc.load_entity("1ubq")
print(f"Atoms: {len(df)}")
print(f"Chains: {proc.get_chains('1ubq')}")

# Save modified structure
proc.save_entity("1ubq_modified", df, metadata={"source": "modified"})

# Export
proc.export_entity("1ubq", Path("output.cif"), format="cif")
```

### Filtering and Modification

```python
# Filter to chain A only
df_chain_a = proc.filter_by_chain("3sn6", chain_ids=["A"])

# Filter to residue range
df_range = proc.filter_by_residue_range("1ubq", "A", start=10, end=50)

# Remove waters and ligands
df_protein = proc.remove_hetatm("1ubq")

# Delete specific residues
proc.delete_residues("1ubq", chain_id="A", residue_ids=[1, 2, 3])
```

### Ligand Analysis

```python
# List ligands
ligands = proc.list_ligands("1ubq", exclude_common=True, min_atoms=5)

# Get binding site residues
interactions = proc.get_ligand_interactions(
    "1ubq",
    ligand_id="ATP",
    cutoff=4.5
)

# Water network analysis
networks = proc.compute_water_networks("1ubq", distance_cutoff=3.5)
```

### Structure Alignment

```python
# Align structures
results = proc.align_structures(
    ["struct_1", "struct_2", "struct_3"],
    reference="struct_1"
)

# Get pairwise RMSD
rmsd = proc.align_pair("struct_1", "struct_2")
print(f"RMSD: {rmsd:.2f} Å")

# All-vs-all alignment
matrix = proc.align_all_vs_all(["s1", "s2", "s3", "s4"])
```

### GRN Annotation

```python
# Annotate with GRN table
proc.annotate_with_grn(
    "3sn6",
    grn_table=grn_df,
    sequence_id="receptor_seq"
)

# Batch annotation
proc.annotate_structures_with_grn(
    ["s1", "s2", "s3"],
    grn_table=grn_df
)
```

### Chain Sequence Extraction

```python
# Get sequence from chain
seq = proc.get_sequence("1ubq", chain_id="A")

# Get all chain sequences
seqs = proc.get_all_sequences("3sn6")

# Register chains as separate entities
proc.register_chain_sequences(
    "3sn6",
    chains=["A", "B"],
    dataset_name="chains_dataset"
)
```

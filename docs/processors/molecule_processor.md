# MoleculeProcessor

Import: `from protos.processing.molecule import MoleculeProcessor`

`MoleculeProcessor` stores small-molecule descriptor records as JSON. Saving a
string treats it as SMILES; saving a mapping preserves its fields and defaults
`kind` to `smiles`.

```python
from protos.processing.molecule import MoleculeProcessor

molecules = MoleculeProcessor()
dataset = molecules.register_smiles_map(
    {"ethanol": "CCO", "benzene": "c1ccccc1"},
    dataset_name="demo_molecules",
)

assert dataset == "demo_molecules"
assert molecules.load_entity("ethanol")["smiles"] == "CCO"
assert set(molecules.load_dataset(dataset)) == {"ethanol", "benzene"}
```

| Method | Current behavior |
| --- | --- |
| `save_entity(name, data, metadata=None)` | write a JSON record and register it |
| `load_entity(name)` | return the JSON record plus registry metadata, or `None` |
| `register_smiles_map(mapping, dataset_name=None, metadata=None)` | replace/create a dataset containing saved records |
| `list_ligands()` | list molecule entity names |
| `calculate_properties(smiles)` | return RDKit descriptors, or `None` without RDKit/for invalid input |
| `filter_drug_like(names, strict=False)` | retain entities passing the implemented RDKit filters |

Descriptor storage does not validate SMILES. RDKit is optional and is required
for actual molecular validation, property calculation, and drug-likeness
filtering. Without RDKit, `calculate_properties()` returns `None` and
`filter_drug_like()` retains nothing.

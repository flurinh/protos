# Unified Data Access

Protos provides unified data access through specialized loaders that acquire data from external sources and automatically register entities with the central registry. Loaders handle data acquisition while processors handle data manipulation.

## Design Principles

1. **Zero Configuration** - Loaders work out of the box with no path parameters required
2. **Human-Readable Names** - All operations use names, not UUIDs
3. **Automatic Registration** - Downloaded data is immediately registered and available
4. **Source Abstraction** - Unified interface across different data sources

## BaseLoader Architecture

All loaders inherit from `BaseLoader` (`protos.io.core.base_loader`), which provides:

- Automatic integration with `ProtosPaths` for directory management
- Automatic integration with `EntityRegistry` for entity tracking
- Automatic integration with `DatasetManager` for bulk operations
- Logging and error handling

```python
from protos.io.core.base_loader import BaseLoader

class CustomLoader(BaseLoader):
    loader_type = "custom"

    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        """Parse and validate identifier."""
        ...

    def fetch_entity(self, identifier: str, **kwargs) -> Optional[Path]:
        """Fetch entity from external source."""
        ...
```

### Core Methods

| Method | Description |
|--------|-------------|
| `fetch_entity(identifier, **kwargs)` | Fetch single entity from external source |
| `parse_identifier(identifier)` | Parse and validate identifier format |
| `download_and_register(identifier, name, metadata)` | Download and register with entity registry |
| `download_batch(identifiers, dataset_name)` | Bulk download with optional dataset creation |
| `validate_identifier(identifier)` | Check if identifier is valid for this loader |
| `list_sources()` | List available data sources |

---

## StructureLoader

Fetches protein structures from RCSB PDB, AlphaFold Database, or local files.

**Location:** `protos.io.ingest.structure_loader`

### Supported Sources

| Source | Identifier Format | Example |
|--------|-------------------|---------|
| RCSB PDB | 4-character alphanumeric | `1ubq`, `3SN6` |
| AlphaFold | UniProt ID or AF-formatted | `P00533`, `AF-P00533-F1` |
| Local | File path | `/path/to/structure.cif` |

### Basic Usage

```python
from protos import StructureLoader

loader = StructureLoader()

# Download from RCSB PDB
name = loader.download_and_register("1ubq")

# Download from AlphaFold
name = loader.download_and_register("P00533", source="alphafold")

# With custom name
name = loader.download_and_register("P00533", name="EGFR_HUMAN", source="alphafold")

# Import local file
name = loader.download_and_register("/path/to/structure.cif", source="local")
```

### Batch Downloads

```python
# Download multiple structures
successful, failed = loader.download_batch(
    ["1ubq", "2w9s", "3sn6"],
    dataset_name="gpcr_structures"
)

print(f"Downloaded: {len(successful)}, Failed: {len(failed)}")
```

### AlphaFold Convenience Method

```python
# Specify AlphaFold version explicitly
name = loader.download_and_register_alphafold(
    "P00533",
    name="EGFR_HUMAN",
    version=4
)
```

### Import Local Files

```python
# Import multiple local files
successful, failed = loader.import_local_structures(
    ["/path/to/file1.cif", "/path/to/file2.pdb"],
    names=["protein_a", "protein_b"]
)

# Register from input folder
successful, failed = loader.register_from_input_folder()
```

### Identifier Parsing

The loader automatically determines the source from the identifier format:

```python
info = loader.parse_identifier("1ubq")
# {'id': '1ubq', 'source': 'rcsb', 'type': 'experimental', 'original_id': '1ubq'}

info = loader.parse_identifier("AF-P00533-F1")
# {'id': 'P00533', 'source': 'alphafold', 'version': 4, 'original_id': 'AF-P00533-F1'}

info = loader.parse_identifier("P00533")
# {'id': 'P00533', 'source': 'alphafold', 'version': 4, 'original_id': 'P00533'}
```

### Source Aliases

Multiple names map to the same canonical source:

| Aliases | Canonical Source |
|---------|-----------------|
| `rcsb`, `pdb`, `cif`, `mmcif` | `rcsb` |
| `alphafold`, `alpha_fold`, `alphafold_db`, `af` | `alphafold` |
| `local`, `file`, `filesystem` | `local` |

---

## SequenceLoader

Fetches protein sequences from UniProt or local FASTA files.

**Location:** `protos.io.ingest.sequence_loader`

### Supported Sources

| Source | Identifier Format | Example |
|--------|-------------------|---------|
| Local FASTA | File path | `/path/to/sequences.fasta` |
| UniProt | UniProt ID or prefixed | `P00533`, `uniprot:P00533` |

### Basic Usage

```python
from protos import SequenceLoader

loader = SequenceLoader()

# Load from local FASTA
name = loader.download_and_register("/path/to/protein.fasta")

# Download from UniProt
name = loader.download_and_register("P00533")

# With explicit prefix
name = loader.download_and_register("uniprot:P00533")

# Multiple UniProt IDs
name = loader.download_and_register("uniprot:P00533,P12345,Q67890")
```

### Register In-Memory Sequences

```python
# Register sequences from Python objects
records = [
    {"name": "protein_a", "sequence": "MKTAYIAKQR...", "metadata": {"source": "custom"}},
    {"name": "protein_b", "sequence": "MVLSPADKTN...", "metadata": {"source": "custom"}},
]

result = loader.register_sequence_records(
    records,
    dataset_name="custom_sequences",
    overwrite=False
)

print(result)
# {'entities': ['protein_a', 'protein_b'], 'dataset': 'custom_sequences'}
```

### Options

| Parameter | Description |
|-----------|-------------|
| `materialize_entities` | Save each sequence as individual entity (default: False for multi-sequence FASTA) |
| `overwrite` | Overwrite existing entities (default: False) |

```python
# Single sequence - materialize as individual entity
name = loader.download_and_register(
    "/path/to/single.fasta",
    materialize_entities=True
)
```

---

## LigandLoader

Imports small molecules from SDF files or SMILES strings, creating both structure and molecule entities.

**Location:** `protos.io.ingest.ligand_loader`

**Dependencies:** RDKit (optional, required for import operations)

### Import from SDF

```python
from protos import LigandLoader

loader = LigandLoader()

# Import all molecules from SDF file
dataset_id, entity_names = loader.import_sdf(
    "/path/to/ligands.sdf",
    dataset_name="my_ligands",
    chain_id="L"  # Chain identifier for structure representation
)

print(f"Created dataset '{dataset_id}' with {len(entity_names)} ligands")
```

### Import from SMILES

```python
# Define molecules as name -> SMILES mapping
smiles_map = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}

dataset_id, entity_names = loader.import_smiles(
    smiles_map,
    dataset_name="common_drugs",
    generate_3d=True,  # Generate 3D coordinates
    chain_id="L"
)
```

### Dual Registration

The LigandLoader creates two types of entities:

1. **Structure entities** - 3D atomic coordinates (via StructureProcessor)
2. **Molecule entities** - Chemical descriptors like SMILES (via MoleculeProcessor)

```python
# These are linked via relationships
from protos import get_registry

registry = get_registry()
related = registry.get_related_entities("aspirin_1", rel_type="has_structure")
```

### Metadata Stored

For each imported ligand:

| Field | Description |
|-------|-------------|
| `entity_type` | Always "ligand" |
| `source_format` | "sdf" or "smiles" |
| `canonical_smiles` | SMILES representation |
| `mol_block` | (SDF only) Original MOL block |
| `source_file` | (SDF only) Original filename |

---

## Other Available Loaders

Protos includes additional specialized loaders:

| Loader | Purpose | Location |
|--------|---------|----------|
| `NCBILoader` | NCBI sequence data | `protos.io.ingest.ncbi_loader` |
| `UniProtLoader` | UniProt protein sequences | `protos.io.ingest.uniprot_loader` |
| `GPCRDBLoader` | GPCR-specific structures | `protos.io.ingest.gpcrdb_loader` |
| `ChemBLLoader` | Chemical compound data | `protos.io.ingest.chembl_loader` |
| `CCDLoader` | Chemical Component Dictionary | `protos.io.ingest.ccd_loader` |
| `EnamineLoader` | Enamine molecule collections | `protos.io.ingest.enamine_loader` |

---

## Complete Example

```python
from protos import (
    StructureLoader,
    SequenceLoader,
    LigandLoader,
    StructureProcessor,
    get_registry
)

# Initialize loaders
struct_loader = StructureLoader()
seq_loader = SequenceLoader()
ligand_loader = LigandLoader()

# Download protein structure and sequence
struct_loader.download_and_register("1ubq", name="ubiquitin")
seq_loader.download_and_register("P0CG48", name="ubiquitin")

# Import ligands
ligand_loader.import_smiles({
    "mg_ion": "[Mg+2]",
    "atp": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)[C@@H](O)[C@H]1O"
}, dataset_name="cofactors")

# All entities are now registered and accessible
registry = get_registry()
print(registry.list_entities())

# Load data via processors
proc = StructureProcessor()
df = proc.load_entity("ubiquitin")
print(df.head())
```

---

## Error Handling

Loaders provide informative error messages and return `None` on failure:

```python
name = loader.download_and_register("invalid_id")
if name is None:
    print("Download failed - check logs for details")

# Batch operations return failed identifiers
successful, failed = loader.download_batch(["1ubq", "XXXX", "2w9s"])
if failed:
    print(f"Failed to download: {failed}")
```

All loaders use Python's logging module:

```python
import logging
logging.basicConfig(level=logging.INFO)

# Loader operations will now print detailed logs
loader = StructureLoader()
loader.download_and_register("1ubq")
```

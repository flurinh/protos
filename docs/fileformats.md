# File Formats

Protos supports a comprehensive range of file formats commonly used in structural biology research. Each format is handled by a specialized processor that manages parsing, validation, and conversion operations.

## Supported Formats

### Structure Formats

#### PDB/mmCIF (.cif, .pdb)
- **Processor**: `CifBaseProcessor`
- **Description**: Protein Data Bank formats for 3D molecular structures
- **Features**:
  - Atomic coordinates
  - Chain and residue information
  - Experimental metadata
  - B-factors and occupancy
  - Secondary structure annotations

```python
# Loading structures
structure = cp.load_structure("1ubq")  # Automatically detects format
structure_cif = cp.load_structure("1ubq.cif")
structure_pdb = cp.load_structure("1ubq.pdb")
```

#### Structure Datasets (.json)
- **Processor**: `CifBaseProcessor`
- **Description**: Collections of related structures
- **Format**:
```json
{
  "dataset_id": "kinases",
  "name": "Human kinase structures",
  "content": ["1atp", "2src", "3erk"],
  "metadata": {
    "organism": "Homo sapiens",
    "family": "protein kinase"
  }
}
```

### Sequence Formats

#### FASTA (.fasta, .fa, .faa)
- **Processor**: `SeqProcessor`
- **Description**: Standard sequence format
- **Features**:
  - Single or multiple sequences
  - Header parsing
  - Automatic entity registration

```python
# FASTA format example
>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEIT
YVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELP
```

#### Sequence Alignments (.aln, .clustal)
- **Processor**: `SeqProcessor`
- **Description**: Multiple sequence alignments
- **Formats**: CLUSTAL, Stockholm, PHYLIP

### GRN Formats

#### GRN Tables (.csv, .tsv)
- **Processor**: `GRNBaseProcessor`
- **Description**: Generic Residue Numbering annotations
- **Format**: Tabular with sequences as rows, positions as columns

```csv
sequence_id,1.50,2.50,3.50,4.50,5.50,6.50,7.50
BACR_HALSA,M62,V90,L129,W171,T205,P238,K257
CHR2_CHLRE,M1,L65,E107,G149,E195,I236,T280
```

#### GRN Configurations (.json)
- **Processor**: `GRNBaseProcessor`
- **Description**: GRN scheme definitions
- **Content**:
  - Position definitions
  - Helix boundaries
  - Functional annotations

### Property Formats

#### Property Tables (.csv, .xlsx)
- **Processor**: `PropertyProcessor`
- **Description**: Experimental properties and metadata
- **Requirements**:
  - Entity identifier column
  - Property columns with consistent types
  - Optional metadata columns

```csv
protein_id,lambda_max,expression_level,thermal_stability
BACR_HALSA,568,high,65.2
CHR2_CHLRE,470,medium,55.8
```

#### Property Datasets (.json)
- **Processor**: `PropertyProcessor`
- **Description**: Structured property collections
- **Features**:
  - Multiple property types
  - Complex metadata
  - Relationships between entities

### Embedding Formats

#### Pickle Files (.pkl)
- **Processor**: `EmbeddingProcessor`
- **Description**: Serialized numpy arrays or tensors
- **Content**: Pre-computed embeddings from various models

#### HDF5 Files (.h5, .hdf5)
- **Processor**: `EmbeddingProcessor`
- **Description**: Hierarchical data format for large embeddings
- **Features**:
  - Compression support
  - Partial loading
  - Metadata storage

## Format Specifications

### Structure Data Schema

CIF/PDB files are parsed into standardized DataFrames:

```python
# Structure DataFrame columns
structure_df = pd.DataFrame({
    'record_type': ['ATOM', 'ATOM', ...],
    'atom_serial': [1, 2, ...],
    'atom_name': ['N', 'CA', ...],
    'residue_name': ['MET', 'MET', ...],
    'chain_id': ['A', 'A', ...],
    'residue_number': [1, 1, ...],
    'x': [10.123, 11.456, ...],
    'y': [20.789, 21.234, ...],
    'z': [30.567, 31.890, ...],
    'occupancy': [1.00, 1.00, ...],
    'b_factor': [20.5, 21.3, ...],
    'element': ['N', 'C', ...]
})
```

### GRN Table Schema

GRN tables follow a specific format:

```python
# GRN DataFrame structure
grn_df = pd.DataFrame({
    # Index: sequence identifiers
    # Columns: GRN positions (helix.position format)
    '1.50': ['M62', 'M1', 'L21', ...],  # Position 50 in helix 1
    '2.50': ['V90', 'L65', 'A65', ...],  # Position 50 in helix 2
    '3.50': ['L129', 'E107', 'E107', ...],  # Position 50 in helix 3
    # ... more positions
}, index=['BACR_HALSA', 'CHR2_CHLRE', 'PROTA_HUMAN', ...])
```

### Property Data Schema

Property data requires entity linkage:

```python
# Property DataFrame requirements
property_df = pd.DataFrame({
    'entity_column': ['PROT1', 'PROT2', ...],  # Links to entities
    'numeric_property': [1.5, 2.3, ...],
    'categorical_property': ['high', 'low', ...],
    'boolean_property': [True, False, ...],
    # Additional columns as needed
})
```

## Format Conversion

### Automatic Conversions

Protos handles many conversions automatically:

```python
# Structure to sequence
structure = cp.load_structure("1ubq")
sequence = cp.extract_sequence(structure)

# Sequence to embeddings
embeddings = emb_proc.embed_sequences({"1ubq": sequence})

# Properties from any source
prop_proc.import_from_csv("experimental_data.csv", entity_column="protein_id")
```

### Cross-Format Operations

```python
# Load structure, extract sequence, generate embedding
structure = struct_proc.load_structure("1ubq")
sequence = struct_proc.extract_sequence(structure, chain="A")
embedding = emb_proc.embed_sequences({"1ubq_A": sequence})

# All linked through entity system
entity_id = registry.resolve_identifier("1ubq_A")
```

## Format Validation

### Structure Validation

```python
# CifBaseProcessor validates:
- Coordinate data types (numeric)
- Required columns present
- Chain ID consistency
- Residue numbering continuity
```

### GRN Validation

```python
# GRNBaseProcessor validates:
- Position format (helix.position)
- Residue+position format in cells
- No invalid helix numbers (0, 8+)
- Consistent sequence lengths
```

### Property Validation

```python
# PropertyProcessor validates:
- Entity column exists
- Data types consistent per column
- No duplicate entity entries
- Numeric ranges reasonable
```

## Custom Format Support

### Adding New Formats

Extend existing processors to support new formats:

```python
class ExtendedSeqProcessor(SeqProcessor):
    def load_genbank(self, file_path: str) -> Dict[str, str]:
        """Load sequences from GenBank format."""
        sequences = {}
        # Custom parsing logic
        with open(file_path) as f:
            # Parse GenBank format
            pass
        
        # Register as entities
        for name, seq in sequences.items():
            self.save_entity(name, seq)
        
        return sequences
```

### Format Detection

Implement automatic format detection:

```python
def detect_format(file_path: Path) -> str:
    """Detect file format from content."""
    with open(file_path) as f:
        first_line = f.readline().strip()
        
    if first_line.startswith('>'):
        return 'fasta'
    elif first_line.startswith('HEADER'):
        return 'pdb'
    elif 'data_' in first_line:
        return 'cif'
    else:
        # Check extension
        return file_path.suffix.lstrip('.')
```

## Format Best Practices

### 1. Use Standard Formats

Prefer widely-supported formats:
- Structures: mmCIF over PDB
- Sequences: FASTA for single sequences
- Tables: CSV with headers
- Metadata: JSON for complex data

### 2. Include Metadata

Embed metadata in files when possible:
```fasta
>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor OS=Homo sapiens GN=EGFR
MRPSGTAGAALLALLAALCPASRA...
```

### 3. Validate Before Processing

```python
# Check file before processing
if not file_path.exists():
    raise FileNotFoundError(f"File not found: {file_path}")

if file_path.stat().st_size == 0:
    raise ValueError(f"Empty file: {file_path}")

# Validate content
try:
    data = processor.load_entity(name)
except Exception as e:
    logger.error(f"Invalid format in {name}: {e}")
```

### 4. Handle Large Files

```python
# Stream large files
def read_large_fasta(file_path: Path):
    """Read FASTA file in chunks."""
    with open(file_path) as f:
        header = None
        sequence_parts = []
        
        for line in f:
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence_parts)
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
        
        if header:
            yield header, ''.join(sequence_parts)
```

## Format Examples

### Complete Workflow Example

```python
# 1. Load structure
structure = cp.load_structure("1ubq")

# 2. Extract sequence
sequence = cp.extract_sequence(structure, chain="A")

# 3. Save sequence
sp.save_sequences({"1ubq_A": sequence}, "extracted.fasta")

# 4. Create GRN annotation
grn_data = pd.DataFrame({
    '3.50': ['I30'],
    '4.50': ['W40'],
    '5.50': ['Y50']
}, index=['1ubq_A'])
gp.save_grn_table("1ubq_annotation", grn_data)

# 5. Add properties
prop_data = pd.DataFrame({
    'protein_id': ['1ubq_A'],
    'stability': [75.5],
    'expression': ['high']
})
pp.import_properties(prop_data, entity_column='protein_id')

# All formats now linked through entity system
```

## Summary

Protos supports:
- Major structural biology file formats
- Automatic format detection and conversion
- Validation and error handling
- Cross-format integration through entities
- Extensible format support

The format abstraction allows researchers to work with familiar file types while leveraging Protos' powerful data management capabilities.
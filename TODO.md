# PROTOS TODO

## 🎯 CORE PHILOSOPHY & CRITICAL PRINCIPLES

### The Protos Promise
**Users work with names, never paths. Protos handles ALL file system complexity.**

### What This Means:
1. **Complete Abstraction** - Users NEVER see or manipulate file paths
2. **Name-Based Access** - Everything accessed by biological/dataset names
3. **Universal Entity IDs** - One protein = one hash ID across all formats
4. **Transparent Operations** - Load/save/convert without format concerns
5. **Registry as Truth** - All lookups and discovery through registry

### The Golden Rules:
```python
# ✅ CORRECT - Users work with names
processor.load_structure("1ABC")           # PDB ID
processor.load_sequence("P12345")          # UniProt ID
processor.load_grn_annotation("BR1_HUMAN") # Sequence name
processor.save_dataset("my_results", data) # Dataset name

# ❌ WRONG - Users never see paths or files
open("/path/to/file.csv")                 # NO direct file access
Path("data/structures/1ABC.cif")          # NO path construction
pd.read_csv("file.csv")                   # NO manual format handling
os.path.exists("some/path")               # NO file system checks
```

### Testing Principles:
0. **USE THE EXISTING MINICONDA ENV: 'protos'** - Activate it during development 
```bash
conda activate protos
```
1. **USE REAL TEST DATA** - Never mock biological data
2. **USE tests/test-data/** - Mirror production data structure  
3. **NO TEMPORARY FILES** - Use fixture data, not generated files
4. **TEST WITH REAL FILES** - CIF, FASTA, CSV from test-data
5. **FOLLOW PROTOS PHILOSOPHY** - Names only, never paths

### Current Violations:
- **170+ path violations in tests** - Breaking core abstraction
- **No universal entity system** - Can't track entities across formats
- **Incomplete registry** - Users can't discover available data
- **Direct file operations** - Bypassing processor abstraction

---

## 📁 PATH SYSTEM - ProtosPaths

**See [docs/PROTOSPATH.md](docs/PROTOSPATH.md) for complete documentation**

### Architecture Summary:
- **ProtosPaths** - Central path resolver (singleton pattern)
- **Global Configuration** - Set once via `ProtosPaths.set_data_root()`
- **Processor Integration** - BaseProcessor uses ProtosPaths automatically
- **Complete Abstraction** - Users never construct or see paths

### Key Components:
```python
# Global configuration (only in conftest.py for tests)
ProtosPaths.set_data_root("/path/to/data")

# Processors use it automatically
processor = CifBaseProcessor()  # No path parameters!
processor.load_structure("1ABC")  # Just names!
```

---

## 🔑 REGISTRY SYSTEM - Entity & Dataset Management

### Universal Entity ID System:
**One biological entity = One hash ID across ALL formats**

```python
# Same protein gets SAME entity ID everywhere
"P12345" → entity_id "a3f2d8c91b"
  ├── sequence format: /seq/a3f2d8c91b.fasta
  ├── structure format: /struct/a3f2d8c91b.cif  
  ├── GRN format: in table with entity_id column
  └── embedding format: /emb/a3f2d8c91b.pkl
```

### Registry Architecture:
1. **EntityRegistry** - Tracks all entities with hash IDs
   - Multi-format support (same entity, multiple formats)
   - Original ID → Hash ID lookups
   - Metadata per format
   
2. **Registry as Universal Translator**:
   ```python
   # User provides biological name
   processor.load_structure("1ABC")
   
   # Registry translates transparently
   entity_id = registry.resolve_identifier("1ABC")  # → "a3f2d8c91b"
   
   # System uses hash internally
   data = processor._load_entity_by_hash(entity_id)
   ```

3. **GRN Tables Special Case**:
   ```
   entity_id    | sequence_id | 1.50 | 2.50 | 3.50
   -------------|-------------|------|------|------
   a3f2d8c91b   | BR1_HUMAN   | L45  | V87  | I123
   b4e5f6a72c   | BR2_MOUSE   | L46  | V88  | I124
   ```

---

## 📋 TODOS BY PRIORITY

### 🔴 HIGH PRIORITY - Registry Implementation

#### Phase 1: Complete Multi-Format EntityRegistry ✅ COMPLETED
- [x] Basic EntityRegistry with hash IDs
- [x] Multi-format support per entity
- [x] Original ID lookups
- [x] Add resolve_identifier method:
  ```python
  def resolve_identifier(self, identifier: str) -> str:
      """Resolve any name to entity hash ID."""
      # Check if already hash
      # Lookup by original ID  
      # Generate new if needed
  ```

#### Phase 2: Update ALL Processors for Entity Support ✅ COMPLETED
- [✅] **CifBaseProcessor**:
  - [x] Basic entity registration
  - [x] load_structure uses resolve_identifier
  - [x] list_structures returns original IDs (not hashes!)
  - [x] Full multi-format awareness
  - [x] save_structure_as_entity method
  - [x] get_entity_id_for_pdb method
  
- [✅] **GRNBaseProcessor** (Special Case - Tables contain MULTIPLE entities):
  ```python
  def save_grn_table_with_entities(self, grn_df):
      # Add entity_id column as FIRST column
      entity_ids = [generate_entity_id(seq_id) for seq_id in grn_df.index]
      grn_df.insert(0, 'entity_id', entity_ids)
      
      # Register each row as an entity
      for seq_id, entity_id in zip(grn_df.index, entity_ids):
          register_entity(entity_id, "grn", seq_id, metadata={"in_table": table_name})
  ```
  
  **Important GRN Design Decisions**:
  - GRN tables contain MANY entities (one per row/sequence)
  - Each row gets an entity_id column based on sequence ID
  - Consider: Should single GRN entities be stored as:
    a) JSON files for individual annotated sequences?
    b) FASTA files (since GRN entity = annotated sequence)?
    c) Single-row CSV tables (current approach)?
  
- [✅] **SeqProcessor**:
  - [x] Parse FASTA headers for entity names
  - [x] Register each sequence as entity
  - [x] Support multi-sequence files
  - [x] load_sequence_entity method
  - [x] save_sequence_entity method
  - [x] list_sequence_entities method
  
- [✅] **EmbeddingProcessor**:
  - [x] Use same entity IDs as source sequences
  - [x] Register embeddings by model type
  - [x] _register_embedding_entities method
  - [x] load_embedding_entity method  
  - [x] list_embedding_entities method

#### Phase 3: Complete Processor Testing ⚡ CURRENT FOCUS
- [ ] Test resolve_identifier for all processors
- [ ] Test list operations return names not hashes
- [ ] Test multi-format entity scenarios
- [ ] Test GRN tables with entity_id column
- [ ] Test that all processors can load by both name and hash ID

### 🟡 MEDIUM PRIORITY - Cross-Format Operations

#### Phase 4: Cross-Format Workflows
- [ ] Sequence → Structure (AlphaFold)
- [ ] Structure → Sequence extraction  
- [ ] Sequence → GRN assignment
- [ ] Any format → Embeddings
- [ ] Track conversion lineage in metadata

#### Phase 5: Migration Tools
- [ ] Script to migrate existing data to entity system
- [ ] Update all test data with proper entities
- [ ] Documentation for users to migrate

### 🟢 LOW PRIORITY - Future Enhancements

#### Phase 6: Advanced Features
- [ ] Entity versioning support
- [ ] Entity relationship tracking
- [ ] Bulk operations optimization
- [ ] Export/import entity registry

#### Phase 7: Documentation & Examples
- [ ] Complete entity system user guide
- [ ] Cross-format workflow examples
- [ ] Migration guide from old system
- [ ] Best practices documentation

---

## 🐛 OTHER ISSUES TO FIX

### Test Issues:
- [ ] Fix 170+ path violations in tests
- [ ] Remove all Path() usage from tests
- [ ] Fix GRN test data setup (column formats)
- [ ] Add proper test markers (@pytest.mark.slow, etc.)
- [ ] Fix embedding test timeouts (use small models)

### CLI Updates:
- [ ] Update CLI embed.py to use new EmbeddingProcessor
- [ ] Ensure all CLI commands use entity system
- [ ] Add entity discovery commands

### Documentation:
- [ ] Update README.md with correct usage
- [ ] Document entity system
- [ ] Add workflow examples

---

## 📊 PROGRESS TRACKING

### Completed ✅:
- ProtosPaths architecture and documentation
- Basic EntityRegistry implementation  
- Multi-format entity support
- Hash-based entity ID generation
- BaseProcessor entity methods

### In Progress 🔄:
- resolve_identifier implementation
- CifBaseProcessor entity integration
- Test updates for new system

### Blocked ❌:
- Cross-format operations (need processors done first)
- Migration tools (need final system design)

---

## 💡 IMPLEMENTATION NOTES

### Entity ID Generation:
```python
# ALWAYS hash the biological identifier, not filename
entity_id = generate_entity_id("P12345")     # UniProt
entity_id = generate_entity_id("1ABC")       # PDB
entity_id = generate_entity_id("BR1_HUMAN")  # Sequence name

# For FASTA with random names
>some_random_name
MKLSPADKTN...
entity_id = generate_entity_id("some_random_name")
```

### Processor Pattern:
```python
class AnyProcessor(BaseProcessor):
    def load_entity(self, identifier: str):
        # User provides name
        entity_id = self.registry.resolve_identifier(identifier)
        
        # Internal uses hash
        return self._load_by_hash(entity_id)
        
    def list_entities(self):
        # Return names, not hashes!
        entities = self.registry.list_entities(format_type=self.format)
        return [self.registry.get_original_id(e) for e in entities]
```

### Timeline:
- Week 1: Complete resolve_identifier and processor updates
- Week 2: Full testing suite  
- Week 3: Cross-format operations
- Week 4: Migration and documentation
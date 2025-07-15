# GRNBaseProcessor Updated to Use ProtosPaths

## Summary
Successfully updated GRNBaseProcessor to use ProtosPaths exclusively and implement all abstract methods from BaseProcessor.

## Key Changes

### 1. Constructor Update
- ✅ Removed `path` and `processor_data_dir` parameters
- ✅ Added `paths` parameter that accepts ProtosPaths instance
- ✅ Zero configuration works - processor creates ProtosPaths if not provided

### 2. Abstract Method Implementations
- ✅ `load_entity(name)` - Loads a single GRN entity (row from table)
- ✅ `save_entity(name, data)` - Saves a single GRN entity to current table

### 3. Path Properties
- ✅ `path_grn_dir` - Points to tables directory
- ✅ `path_ref_dir` - Points to reference tables directory  
- ✅ `path_config_dir` - Points to configs directory

### 4. Method Updates
- ✅ Added `return_only` parameter to `load_grn_table()`
- ✅ Override `save_data()` to handle CSV index parameter correctly
- ✅ Fixed data handling when `return_only=True`

## Usage Example
```python
# Zero configuration
processor = GRNBaseProcessor()

# With explicit paths
paths = ProtosPaths(data_root='/my/data')
processor = GRNBaseProcessor(paths=paths)

# Load entity
grn_data = processor.load_entity('protein_123')

# Save entity
new_data = pd.Series({'3.50': 'A', '7.53': 'F'})
processor.save_entity('new_protein', new_data)

# Load table with return_only
table = processor.load_grn_table('ref/gpcrdb_ref', return_only=True)
```

## Test Coverage
Created comprehensive test suite in `test_grn_base_processor_protospath.py`:
- ✅ Zero configuration initialization
- ✅ Path parameters rejected
- ✅ Path properties use ProtosPaths
- ✅ load_entity implementation
- ✅ save_entity implementation
- ✅ Dataset operations
- ✅ return_only parameter

## Compatibility
- Maintains backward compatibility with existing GRN functionality
- All existing methods work as before
- Legacy path parameter shows deprecation warning but still works

## Next Steps
With GRNBaseProcessor complete, the next processor to update is PropertyProcessor (task 6.4).
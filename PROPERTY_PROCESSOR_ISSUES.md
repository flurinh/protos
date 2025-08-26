# PropertyProcessor Critical Issues

## Summary
PropertyProcessor violates several core data management principles outlined in DATA_MANAGEMENT_UNIFIED.md:

### 1. **Hash IDs Exposed to Users** ❌
- `list_entities()` returns hash IDs like `'4a0145b85d'` instead of human names like `'protein_A'`
- The processor logs show hash IDs: "Assigned property 'lambda_max' = 568 to entity 1b911509e5"
- This violates the principle: "Hash IDs are NEVER exposed to users"

### 2. **Internal Structure Issues**
- PropertyProcessor stores entities by hash ID in `self.property_datasets[dataset_name]['entities']`
- Should store by human-readable name instead
- The `_resolve_entity_identifier` method generates hash IDs when it should preserve human names

### 3. **Inconsistent with Other Processors**
- SeqProcessor and GRNBaseProcessor correctly use human-readable names
- PropertyProcessor is the outlier that needs fixing

## Required Fixes

### Fix 1: Store Entities by Human Name
```python
# Current (WRONG):
self.property_datasets[dataset_name]['entities'][hash_id] = properties

# Should be:
self.property_datasets[dataset_name]['entities'][human_name] = properties
```

### Fix 2: Update list_entities to Return Human Names
```python
def list_entities(self, dataset_name: Optional[str] = None) -> List[str]:
    """List entities - returns human-readable names."""
    if dataset_name:
        if dataset_name in self.property_datasets:
            # Return the human names directly
            return list(self.property_datasets[dataset_name].get('entities', {}).keys())
    # ...
```

### Fix 3: Fix _resolve_entity_identifier
Instead of always generating hash IDs, it should:
1. Check if entity is already registered -> return original name
2. If not registered -> return the original identifier as-is
3. Only use hash IDs for internal registry lookups, never for storage

### Fix 4: Update All Internal Storage
- `self.entity_properties` should key by human name, not hash ID
- All methods should work with human names throughout
- Hash IDs only used when interfacing with EntityRegistry internally

## Impact
These fixes are critical because:
1. **User Experience**: Users should NEVER see hash IDs
2. **Data Integrity**: Files and datasets should be human-readable
3. **Consistency**: All processors must follow the same pattern
4. **Debugging**: Human names make debugging and data inspection easier

## Test Failures Explained
The test failures all stem from this core issue:
- `test_list_operations` fails because it gets hash IDs instead of 'protein_A', 'protein_B', etc.
- `test_integration_with_entity_registry` fails because of hash ID usage
- Other tests fail due to cascading effects of this core problem
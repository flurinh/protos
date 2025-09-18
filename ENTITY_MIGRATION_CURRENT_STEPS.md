# Entity System Migration - Current Steps

## Active Work: Phase 1 - Core Infrastructure

### What We're Building
Upgrading the entity registry from hash-based IDs to UUID-based IDs while maintaining complete backward compatibility and ensuring users only ever see human-readable names.

### Changes Made (Jan 18, 2025)

#### Phase 1: Entity Registry
1. ✅ Modified `_generate_hash_id()` to use `uuid.uuid4()` instead of SHA256 hash
2. ✅ Added `name_index` to registry JSON for O(1) name lookups
3. ✅ Updated all internal references from `hash_id` to `entity_id`
4. ✅ Renamed `_resolve_to_hash` to `_resolve_to_id`
5. ✅ All entity registry tests pass with new UUID system
6. ✅ Implemented relationship tracking with clean model:
   - Relationships stored as list of objects with consistent schema
   - Only canonical direction stored (inverse derived at query time)
   - Relationship type registry defines inverses and symmetry
   - APIs only expose human-readable names, never UUIDs
7. ✅ Added comprehensive relationship tests (11 tests all passing)

#### Phase 2: Dataset Manager
8. ✅ Updated DatasetManager to store entity IDs internally:
   - Datasets now have `entity_ids` array alongside `entities`
   - All operations work with UUIDs internally for stability
9. ✅ Maintained backward compatibility:
   - Old datasets without `entity_ids` still work
   - New datasets store both for compatibility
10. ✅ Handle entity lifecycle properly:
    - Renamed entities show new names in datasets
    - Deleted entities show historic names or placeholders
11. ✅ Added comprehensive dataset tests (7 tests passing)

#### Phase 3: BaseProcessor Integration
12. ✅ Updated BaseProcessor.load_dataset to use get_dataset_entities:
    - Correctly resolves entity IDs to current names
    - Handles renamed and deleted entities gracefully
13. ✅ Fixed BaseProcessor weaknesses identified in review:
    - `_sanitize_filename` now mandatory for all file operations
    - Error handling uses logger instead of print
    - Path resolution centralized through ProtosPaths
    - Clarified metadata scope (processor_metadata)
14. ✅ Added comprehensive processor integration tests (12 tests passing)

### Current Status (Jan 19, 2025)

#### Phase 0: Safe Registration Complete ✅
15. ✅ Created InputManager for safe file registration:
    - Replaces unsafe "drag-and-drop" with validated input folder workflow
    - Files placed in `data/input/` are validated and registered
    - Duplicate detection using content hashes
    - Conflict resolution strategies (skip, version, replace)
    - Processed files archived with timestamps
    - Rejected files moved with error logs
16. ✅ Created centralized FormatRegistry:
    - Single source of truth for file formats
    - ProcessorType and FormatCategory enums
    - Automatic format detection by extension
    - Replaced hardcoded format mappings across codebase
17. ✅ Created RegistryHealthChecker:
    - Find unregistered files in data directories
    - Check dataset integrity
    - Validate entity-file relationships
    - Generate comprehensive health reports

### Current Focus: Update Individual Processors

Next steps involve updating each processor to fully utilize the new entity system:
1. Remove deprecated generate_entity_id usage
2. Ensure all processors use EntityRegistry for ID management
3. Update save/load methods to work with new system
4. Add relationship tracking where appropriate

### Previous Work: Migration Tests ✅

Created comprehensive tests in `migration_tests/` folder to ensure:
1. ProtosPaths handles all data structure setup
2. Entity registry correctly manages UUID/name mappings
3. Users never see entity IDs (only human names)
4. File-based and table-based processors work correctly

---

## Step 1: Test Infrastructure Setup ✅ COMPLETE

### 1.1 Create Base Test Class ✅ COMPLETE
- [x] Create `test_base.py` with ProtosPaths setup
- [x] Ensure temporary test directories
- [x] Helper methods for entity creation
- [x] Cleanup after tests

### 1.2 Test Entity Registry Core ✅ COMPLETE
- [x] Test UUID generation (not hash-based) - Currently fails, needs migration
- [x] Test name-to-ID mapping
- [x] Test alias resolution
- [x] Test duplicate detection
- [x] Test multi-format entities
- [x] Test registry persistence
- [x] Test only human names returned

### 1.3 Test File Registration Safety ✅ COMPLETE
- [x] Test content hash detection
- [x] Test name conflict handling
- [x] Test safe registration workflow
- [x] Test dual format handling (CIF/PKL)
- [x] Test table-based registration
- [x] Test input folder workflow
- [x] Test registry health checks

### Test Results Summary
- ✅ 11 tests passing (current EntityRegistry functionality works)
- ⏭️ 4 tests skipped (waiting for new features)
- 🔴 1 test expects UUID implementation

The tests confirm:
1. Current system uses 10-character hash IDs
2. Human names are properly returned in all public methods
3. Multi-format entities work correctly
4. Registry persistence works

---

## Step 2: Update EntityRegistry Implementation ✅ COMPLETE

### 2.1 Core Changes
- [x] Replace hash ID with UUID generation ✅ COMPLETE
- [x] Add name-to-ID index for fast lookups ✅ COMPLETE
- [ ] Implement content hashing for duplicates (deferred to later phase)
- [x] Add relationship tracking ✅ COMPLETE

### 2.2 Registry Structure (Implemented)
```json
{
  "entities": {
    "550e8400-e29b-41d4-a716-446655440000": {
      "original_id": "1UBQ",
      "aliases": ["UBIQ_HUMAN", "P62988"],
      "formats": {
        "structure": {
          "file_path": "structure/mmcif/1UBQ.cif",
          "metadata": {"format": "cif"},
          "created": "2025-01-18T12:00:00"
        }
      },
      "relationships": [
        {
          "type": "derived_from",
          "source": "660e8400-e29b-41d4-a716-446655440001",
          "target": "550e8400-e29b-41d4-a716-446655440000",
          "metadata": {"method": "to_pkl"},
          "created": "2025-01-18T13:00:00"
        }
      ],
      "created": "2025-01-18T12:00:00",
      "modified": "2025-01-18T13:00:00"
    }
  },
  "name_index": {
    "1UBQ": "550e8400-e29b-41d4-a716-446655440000",
    "UBIQ_HUMAN": "550e8400-e29b-41d4-a716-446655440000",
    "P62988": "550e8400-e29b-41d4-a716-446655440000"
  }
}
```

---

## Step 3: Update Dataset Manager ✅ COMPLETE

### 3.1 Core Changes
- [x] Datasets store entity IDs internally for stability
- [x] All APIs continue to use human-readable names
- [x] Backward compatibility maintained with old datasets
- [x] Handle deleted entities gracefully

### 3.2 Implementation Details
- Datasets now store both `entities` (names) and `entity_ids` (UUIDs)
- `get_dataset_entities()` resolves IDs to current names
- Deleted entities show historic names or placeholder
- All operations (add/remove) work with entity IDs internally
- 7 comprehensive tests passing

---

## Step 4: Test Processor Integration ✅ COMPLETE

### 4.1 BaseProcessor Integration
- [x] Processors correctly register entities with UUID system
- [x] Dataset operations work through processors
- [x] Entity relationships are tracked
- [x] Human names used in all processor APIs
- [x] Entity renames handled correctly
- [x] Processor-specific datasets work

### 4.2 BaseProcessor Fixes Applied
- [x] `_sanitize_filename` now used consistently in save_data/load_data
- [x] Complex identifiers (sp|P12345|TEST_HUMAN) handled correctly
- [x] Error handling uses logger.warning instead of print
- [x] Path resolution improved in delete_entity
- [x] Renamed metadata to processor_metadata for clarity
- [x] 12 comprehensive tests all passing

---

## Key Principles to Maintain

1. **ProtosPaths manages ALL paths** - No hardcoded paths anywhere
2. **Users see ONLY human names** - Entity IDs are completely internal
3. **Backward compatibility** - Existing code must continue to work
4. **Safe registration** - Prevent data loss and duplicates

---

## Progress Tracking

- 🔴 Not Started
- 🟡 In Progress  
- 🟢 Complete
- ⚠️ Blocked

### Today's Goals
1. ✅ Set up test infrastructure
2. ✅ Create core registry tests
3. ✅ Implement UUID-based registry
4. ✅ Update Dataset Manager
5. 🟡 Test with real processors (Next step)

---

## Notes & Decisions

- Using UUIDs instead of content hashes for entity IDs (more stable)
- name_index provides O(1) lookup performance
- Dual format tracking (CIF/PKL) requires special handling
- Table-based processors need different entity concept
- **UUID Implementation Complete**: EntityRegistry now generates UUIDs instead of 10-char hashes
- **Name Index Implemented**: Fast O(1) lookups using name_index mapping names/aliases to UUIDs
- All existing tests pass with new UUID system
- **Relationship Tracking Complete**: Clean implementation following user's improved model:
  - Relationships as first-class objects with consistent schema
  - Only canonical direction stored, inverses derived at query time
  - RELATIONSHIP_TYPES registry defines relationship semantics
  - No entity IDs exposed in public APIs
- **Dataset Manager Updated**: Now uses entity IDs internally while maintaining human-readable APIs:
  - Datasets store both `entities` and `entity_ids` for stability and compatibility
  - Handles entity renames transparently (shows current names)
  - Gracefully handles deleted entities (shows historic names)
  - Full backward compatibility with existing datasets
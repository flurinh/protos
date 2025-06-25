# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **Note**: This file contains two sections:
> 1. **General Session Management** (lines 5-128): Apply these procedures if a `/claude` directory exists with tracking files
> 2. **Protos-Specific Guidance** (lines 131+): Always apply these when working with the Protos codebase

## 🔴 MANDATORY SESSION START PROCEDURE (If /claude directory exists)
**EVERY new session MUST begin with:**
1. **READ** `/claude/execution.txt` - Understand system setup and activation procedures
2. **READ** `/claude/content.txt` - Get current project structure
3. **READ** `/claude/logs.txt` - Check last 20-30 entries for recent context (use tail -n 30)
4. **READ** `/claude/todo.txt` - Review project-level task trees and pending work
5. **READ** `/claude/csv_tables.txt` - Understand CSV handling rules and context-aware operations

**This is NOT optional - Do this BEFORE any other action**

## 🚨 CRITICAL RULES - CHECK EVERY TIME
0. **EVERY FILE YOU CREATE IF POSSIBLE YOU SHOULD TEST IMMEDIATLY**
1. **LOG EVERY FILE CHANGE** → Append to `/claude/logs.txt` → `[YYYY-MM-DD HH:MM] [CREATE/MODIFY/DELETE] file - message`
2. **UPDATE PROJECT STRUCTURE** → `/claude/content.txt` → After any file creation/deletion
3. **TEMP FILES** → Prefix with `temp_` or `dev_`
4. **DEPRECATED FILES** → Prefix with `deprecated_`
5. **BEFORE ANY FILE OP** → Check `/claude/content.txt` first
6. **NO SUDO COMMANDS** → Claude cannot run sudo commands. User must manually run:
   - `sudo service postgresql start` (start database)
   - Any other system-level commands requiring sudo

## 🎯 PROJECT QUICK INFO

## 📁 KEY PATHS
```
/claude/logs.txt          → File change log (MANDATORY)
/claude/content.txt       → Current architecture
/claude/todo.txt          → Project task trees (XML format)
/claude/csv_tables.txt    → CSV handling rules and context
/claude/csv_helper.py     → CSV operations helper script
```

## 🔧 COMMON COMMANDS

### CSV Operations (Memory-Efficient)
```bash
# First understand table structure
python /claude/csv_helper.py schema output/rmsd_matrix.csv

# Get row by index as dictionary
python /claude/csv_helper.py get-row output/msa_table_grn.csv --index "PDB_001"

# Get specific cell value
python /claude/csv_helper.py get-cell output/distance_table_grn.csv --index "PDB_001" --column "3.50"

# Search indices
python /claude/csv_helper.py search-index output/protein_summary.csv --pattern ".*_exp$"
```

## 🔍 LOG SEARCH PATTERNS
Use Grep tool on `/claude/logs.txt` to find:
- File history: `grep "filename.py" /claude/logs.txt`
- Recent creates: `grep "\[CREATE\]" /claude/logs.txt | tail -20`
- Recent modifies: `grep "\[MODIFY\]" /claude/logs.txt | tail -20`
- Specific date: `grep "2025-06-23" /claude/logs.txt`
- Deprecated files: `grep "\[DEPRECATE\]" /claude/logs.txt`
- **View recent logs**: `tail -n 30 /claude/logs.txt`

## 📊 CONTEXT BUILDING
The logs provide:
- **Change history** for any file
- **Development timeline** 
- **Pattern detection** (repeated modifications)
- **Deprecation tracking**
- **Collaboration context** (who changed what when)

## 📝 LOG FORMAT & INSTRUCTIONS
**When logging changes:**
- Format: `[YYYY-MM-DD HH:MM] [ACTION] /path/to/file - commit message`
- Actions: CREATE, MODIFY, DELETE, RENAME (old → new), DEPRECATE
- Always APPEND to `/claude/logs.txt` (never overwrite)
- Use bash: `echo "[$(date +%Y-%m-%d' '%H:%M)] [ACTION] /path - message" >> /claude/logs.txt`

## 🎯 TASK TREE MANAGEMENT
**Project-level tasks in `/claude/todo.txt`:**
- XML format with nested task trees
- Break complex tasks into subtasks
- Status: pending | in_progress | completed | blocked
- Complete all subtasks before marking parent complete
- When solving subtasks, continue until parent task is resolved
- Update XML structure as tasks progress

### Task Tree XML Structure
Each task element can have:
- `id`: unique identifier (e.g., "1", "backend-auth", "feature.1.2")
- `status`: pending | in_progress | completed | blocked
- `priority`: high | medium | low
- `description`: what needs to be done
- `subtasks`: child tasks that must be completed first

Example format:
```xml
<task id="feature-1" status="in_progress" priority="high">
  <description>Implement user authentication system</description>
  <subtasks>
    <task id="feature-1.1" status="completed" priority="high">
      <description>Create user model with UUID</description>
    </task>
    <task id="feature-1.2" status="pending" priority="medium">
      <description>Add login/logout endpoints</description>
    </task>
  </subtasks>
</task>
```

### Task Completion Rules
1. Mark task "in_progress" when starting work
2. Complete ALL subtasks before marking parent complete
3. If blocked, add reason in description
4. When completing a task, check if parent can be completed
5. Always update todo.txt when task status changes

### Task Breakdown Guidelines
**Complex tasks should be broken down when:**
- Task requires multiple distinct operations
- Task spans multiple files or modules
- Task has dependencies on other components
- Task can be parallelized into independent subtasks

**Atomic tasks (no breakdown needed):**
- Single file edits
- Simple function implementations
- Configuration changes
- Documentation updates

---

## 🧬 PROTOS-SPECIFIC GUIDANCE

The following sections provide guidance specific to the Protos structural biology library.

## Development Commands

### Installation
```bash
# Install in development mode (required for proper path resolution)
pip install -e .

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "import protos; print(protos.__version__)"
```

### Testing
```bash
# Run all tests (uses pytest)
pytest

# Run tests from old_tests directory (legacy tests)
python run_tests.py -m test_processing/test_grn

# Run tests with verbose output
python run_tests.py -v

# Run tests with coverage
python run_tests.py -c

# Run a single test file
pytest tests/test_io/test_fasta_utils.py

# Run a single test
pytest tests/test_io/test_fasta_utils.py::test_specific_function

# Note: Tests are split between:
# - tests/        (new tests)
# - old_tests/    (legacy tests, use run_tests.py)
```

### Code Quality
```bash
# Format code with Black
black src/

# Sort imports with isort
isort src/

# Type check with mypy
mypy src/
```

## Architecture Overview

Protos is a structural biology library with a modular processor architecture designed for handling diverse biological data types.

### Core Processor System

The library uses specialized **Processors** that inherit from `BaseProcessor`:

- **CifBaseProcessor**: Handles 3D protein structures (PDB/CIF files)
- **GRNBaseProcessor**: Manages Generic Residue Numbering for consistent residue mapping
- **SeqProcessor**: Processes sequences and alignments
- **EMBProcessor**: Generates and manages protein embeddings
- **PropertyProcessor**: Integrates metadata and calculated properties
- **LigProcessor**: Handles ligand information

Processors are designed to be interoperable - outputs from one can serve as inputs for another, enabling complex analysis pipelines.

### Critical Path Management

**IMPORTANT**: Protos has a path management system that requires explicit configuration to avoid bugs:

```python
import os
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths

# ALWAYS set environment variables explicitly
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())

# Initialize paths with explicit directories
paths = ProtosPaths(
    user_data_root=str(data_dir.absolute()),
    ref_data_root=str(data_dir.absolute()),
    create_dirs=True
)

# When creating processors, always pass explicit paths
processor = CifBaseProcessor(
    name="my_processor",
    data_root=str(data_dir.absolute()),
    processor_data_dir="structure"
)
```

The path system has a known bug with `DataSource.AUTO` defaulting to reference data when it should use user data. Always use explicit paths to avoid this issue.

### Data Organization

```
protos_data/
├── structure/               # Structure data (CifBaseProcessor)
│   ├── mmcif/              # PDB/CIF files
│   ├── alignments/         # Structure alignments
│   └── structure_dataset/  # Dataset definitions
├── grn/                    # GRN data (GRNBaseProcessor)
│   ├── ref/                # Reference GRN tables
│   ├── tables/             # Calculated GRN tables
│   └── configs/            # GRN configurations
├── sequence/               # Sequence data (SeqProcessor)
│   ├── fasta/              # FASTA files
│   └── alignments/         # Sequence alignments
├── embedding/              # Embeddings (EMBProcessor)
└── property/               # Properties (PropertyProcessor)
```

### Key Implementation Patterns

1. **Multi-Stage Caching**: Production code uses pickle files for caching intermediate results
2. **Type Enforcement**: Always ensure coordinates are numeric with `pd.to_numeric()`
3. **Cross-Platform Paths**: Use `pathlib.Path` for Windows/Linux compatibility
4. **Error Handling**: Check if DataFrames are empty before operations
5. **Reference vs User Data**: Processors check user data first, then fall back to reference data

### Common Workflow Example

```python
# 1. Set up paths explicitly
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# 2. Load structures with type enforcement
cp = CifBaseProcessor(name="analysis", data_root=str(data_dir.absolute()))
cp.load_dataset("my_dataset", apply_dtypes=True)

# 3. Filter and process
for pdb_id in cp.pdb_ids:
    df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
    df_chain = df[df['auth_chain_id'] == 'A'].copy()
    
    if not df_chain.empty:
        # Ensure numeric coordinates
        for coord in ['x', 'y', 'z']:
            df_chain[coord] = pd.to_numeric(df_chain[coord], errors='coerce')

# 4. Perform analysis with GRN mapping
grnp = GRNBaseProcessor(name="grn", data_root=str(data_dir.absolute()))
grnp.load_grn_table("reference_grn_table")
```

### External Tool Integration

The library integrates with bioinformatics tools:
- **FoldMason**: Structure alignment (use `use_wsl=True` on Windows)
- **GTalign**: GPU-accelerated alignment
- **MMseqs2**: Sequence searching

### Important Notes

1. **Always use explicit paths** - The automatic path resolution has bugs
2. **Check empty DataFrames** - Many operations fail silently on empty data
3. **Enforce data types** - Use `apply_dtypes=True` when loading datasets
4. **Handle WSL paths** - Windows paths need conversion for Linux tools
5. **Test data exists** - Check `tests/test-data/` for test fixtures

The codebase is transitioning from a dual-path system (separate reference/user data) to a unified data directory. Some legacy code may still reference the old system.

## Common Issues & Solutions

### Path Resolution Errors
```python
# ❌ WRONG - Relies on buggy automatic path resolution
processor = CifBaseProcessor(name="test")

# ✅ CORRECT - Explicit path configuration
import os
data_dir = Path("/absolute/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
processor = CifBaseProcessor(name="test", data_root=str(data_dir.absolute()))
```

### Empty DataFrame Operations
```python
# ❌ WRONG - Fails silently on empty DataFrames
df_chain = df[df['chain_id'] == 'A']
ca_atoms = df_chain[df_chain['atom_name'] == 'CA']  # May be empty!

# ✅ CORRECT - Check for empty DataFrames
df_chain = df[df['chain_id'] == 'A'].copy()
if df_chain.empty:
    print(f"No chain A found")
    continue
```

### Coordinate Type Issues
```python
# ❌ WRONG - Coordinates may be strings
coords = df[['x', 'y', 'z']].values

# ✅ CORRECT - Ensure numeric types
for coord in ['x', 'y', 'z']:
    df[coord] = pd.to_numeric(df[coord], errors='coerce')
coords = df[['x', 'y', 'z']].values
```

### Windows/WSL Path Issues
```python
# When using external tools on Windows
from protos.processing.structure.foldmason import FoldMason
fm = FoldMason(use_wsl=True)  # Required on Windows

# Convert Windows paths for Linux tools
win_path = "C:\\Users\\user\\data"
wsl_path = "/mnt/c/Users/user/data"
```
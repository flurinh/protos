# Test Data Directory

This directory contains test data for the Protos framework using the simplified path system.

## Structure

- All data is contained in a single directory hierarchy
- No separation between 'reference' and 'user' data
- Automatic directory creation when processors are initialized

## Usage

```python
from protos.io.paths.path_config_simplified import ProtosPaths

# Initialize with this test directory
paths = ProtosPaths(data_root='tests/test-data')
```

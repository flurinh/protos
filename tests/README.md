# Protos Test Suite

This directory contains the test suite for the Protos package. The tests are organized to match the structure of the core package.

## Running Tests

To run all tests:

```bash
# From the protos directory
pytest
```

To run tests for a specific module:

```bash
# Run only the GRN module tests
pytest tests/test_processing/test_grn

# Run only the embedding tests
pytest tests/test_embedding
```

For more verbose output:

```bash
pytest -v
```

To get code coverage information:

```bash
pytest --cov=protos
```

## Test Organization

The test suite is organized to mirror the structure of the protos package:

- **test_processing/**: Tests for the processing submodules
  - **test_grn/**: Tests for GRN-related functionality
  - **test_structure/**: Tests for structure processing
  - **test_sequence/**: Tests for sequence alignment and analysis
  - ...
- **test_embedding/**: Tests for embedding functionality
- **test_io/**: Tests for file I/O utilities
- **test_visualization/**: Tests for visualization components
- **test_loaders/**: Tests for data loading utilities

## Writing New Tests

When adding new functionality to the Protos package, please also add corresponding tests. Tests should be placed in the appropriate directory matching the structure of the functionality being tested.

Each test file should:
1. Be named with a `test_` prefix, e.g. `test_grn_utils.py`
2. Import the necessary pytest fixtures and functions
3. Include clear docstrings explaining what each test case is verifying
4. Cover expected behavior as well as edge cases and error conditions

Example:

```python
"""
Tests for some_module functionality
"""

import pytest
from protos.some_module import some_function

def test_some_function():
    """Test that some_function works as expected."""
    result = some_function()
    assert result == expected_value
```
#!/usr/bin/env python
"""
Test runner script for the protos package.
"""

import os
import sys
import subprocess
import argparse

def run_tests(module=None, verbose=False, coverage=False):
    """Run pytest with the given options."""
    pytest_args = ["pytest"]
    
    # Add target module if specified
    if module:
        pytest_args.append(f"old_tests/{module}")
    
    # Add verbosity flag if requested
    if verbose:
        pytest_args.append("-v")
    
    # Add coverage flag if requested
    if coverage:
        pytest_args.append("--cov=protos")
    
    # Run the old_tests
    return subprocess.call(pytest_args)

def main():
    """Parse arguments and run old_tests."""
    parser = argparse.ArgumentParser(description="Run protos package old_tests")
    parser.add_argument(
        "-m", "--module", 
        help="Run old_tests for a specific module (e.g., 'test_processing/test_grn')"
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Enable verbose output"
    )
    parser.add_argument(
        "-c", "--coverage", 
        action="store_true", 
        help="Generate coverage report"
    )
    
    args = parser.parse_args()
    
    # Ensure we're in the protos directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    result = run_tests(
        module=args.module,
        verbose=args.verbose,
        coverage=args.coverage
    )
    
    sys.exit(result)

if __name__ == "__main__":
    main()
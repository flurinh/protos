"""
File utility functions for Protos.

This module provides simple utility functions for common file operations.
For more complex format handling, use the handlers in formats.py.
"""

import os
import json
import pickle
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List

def ensure_dir(directory: str) -> None:
    """Ensure that a directory exists; if it doesn't, create it."""
    Path(directory).mkdir(parents=True, exist_ok=True)
    
def save_json(data: Any, filepath: str, **kwargs) -> None:
    """
    Save data as a JSON file.
    
    Args:
        data: Data to save
        filepath: Path to save to
        **kwargs: Additional arguments for json.dump
    """
    ensure_dir(os.path.dirname(filepath))
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=kwargs.get('indent', 2))
        
def load_json(filepath: str) -> Any:
    """
    Load data from a JSON file.
    
    Args:
        filepath: Path to load from
        
    Returns:
        Loaded data
    """
    with open(filepath, 'r') as f:
        return json.load(f)
        
def save_pickle(data: Any, filepath: str, protocol: int = None) -> None:
    """
    Save data as a pickle file.
    
    Args:
        data: Data to save
        filepath: Path to save to
        protocol: Pickle protocol version
    """
    ensure_dir(os.path.dirname(filepath))
    with open(filepath, 'wb') as f:
        pickle.dump(data, f, protocol=protocol or pickle.HIGHEST_PROTOCOL)
        
def load_pickle(filepath: str) -> Any:
    """
    Load data from a pickle file.
    
    Args:
        filepath: Path to load from
        
    Returns:
        Loaded data
    """
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def validate_fasta_format(filepath):
    """Validate that a file is in FASTA format."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    # Check if file starts with a header line
    if not lines[0].startswith('>'):
        return False
        
    # Check that alternating lines are headers and sequences
    has_sequence = False
    for i, line in enumerate(lines):
        if i == 0:
            continue  # Skip first line, already checked
            
        line = line.strip()
        if not line:
            continue  # Skip empty lines
            
        if has_sequence:
            if line.startswith('>'):
                has_sequence = False
            else:
                has_sequence = True
        else:
            if line.startswith('>'):
                has_sequence = False
            else:
                has_sequence = True
                
    return True

def get_filenames(path):
    """Get filenames in directory without extension.
    
    Args:
        path (str): Directory path to search for files
        
    Returns:
        list: List of filenames without extensions
    """
    if not os.path.exists(path):
        return []
    return [os.path.splitext(f)[0] for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
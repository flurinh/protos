import os
import json
import pickle
import pandas as pd
from pathlib import Path

def ensure_dir(directory):
    """Ensure that a directory exists; if it doesn't, create it."""
    Path(directory).mkdir(parents=True, exist_ok=True)
    
def save_json(data, filepath):
    """Save data as a JSON file."""
    ensure_dir(os.path.dirname(filepath))
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
        
def load_json(filepath):
    """Load data from a JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)
        
def save_pickle(data, filepath):
    """Save data as a pickle file."""
    ensure_dir(os.path.dirname(filepath))
    with open(filepath, 'wb') as f:
        pickle.dump(data, f)
        
def load_pickle(filepath):
    """Load data from a pickle file."""
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
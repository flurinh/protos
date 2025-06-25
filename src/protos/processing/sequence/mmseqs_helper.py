"""
MMseqs2 helper functions for WSL compatibility
"""

import os
import subprocess
import platform
from pathlib import Path
from typing import Dict


def is_wsl():
    """Check if running in WSL."""
    return 'microsoft' in platform.uname().release.lower()

def is_windows():
    """Check if running on Windows."""
    return platform.system() == 'Windows'

def get_mmseqs_command():
    """Get the appropriate MMseqs2 command based on the environment."""
    
    # Check if MMSEQS_PATH is set
    mmseqs_path = '~/MMseqs2/build/bin/mmseqs'
    print("mmseqs path: {}".format(mmseqs_path))
    
    if mmseqs_path:
        # If running on Windows, use WSL
        if is_windows():
            return ['wsl', mmseqs_path]
        else:
            # On Linux/WSL, use path directly
            return [mmseqs_path]
    
    # Try to find mmseqs in PATH
    try:
        # Check if mmseqs is available directly
        result = subprocess.run(['which', 'mmseqs'], capture_output=True, text=True)
        if result.returncode == 0:
            return ['mmseqs']
    except:
        pass
    
    # Common installation paths
    common_paths = [
        '/usr/local/bin/mmseqs',
        '/usr/bin/mmseqs',
        '~/mmseqs2/bin/mmseqs',
        os.path.expanduser('~/MMseqs2/build/bin/mmseqs')
    ]
    
    for path in common_paths:
        expanded_path = os.path.expanduser(path)
        if os.path.exists(expanded_path):
            return [expanded_path]
    
    # If nothing found, return None
    return None

def run_mmseqs_safe(query_seq: str, ref_seqs: Dict[str, str], temp_folder: str = 'temp'):
    """
    Run MMseqs2 safely, returning None if not available.
    """
    mmseqs_cmd = get_mmseqs_command()
    
    if not mmseqs_cmd:
        print("MMseqs2 not found. Using fallback method.")
        return None
    
    try:
        # Import the original function
        from protos.processing.sequence.seq_alignment import mmseqs2_align
        
        # Temporarily set the command in environment
        old_path = os.getenv("MMSEQS_PATH")
        os.environ["MMSEQS_PATH"] = mmseqs_cmd[-1]  # Last element is the actual path
        
        try:
            # If we need to use wsl, modify subprocess calls in mmseqs2_align
            # For now, just try to run it
            result = mmseqs2_align(query_seq, ref_seqs, temp_folder)
            return result
        finally:
            # Restore old path
            if old_path:
                os.environ["MMSEQS_PATH"] = old_path
            else:
                del os.environ["MMSEQS_PATH"]
                
    except Exception as e:
        print(f"MMseqs2 error: {e}")
        return None
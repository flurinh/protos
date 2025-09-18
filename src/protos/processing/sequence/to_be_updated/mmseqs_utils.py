"""Utilities for MMseqs2 integration with automatic path detection."""

import os
import subprocess
import platform
from typing import Optional, Tuple
import re


def detect_mmseqs2() -> Tuple[Optional[str], bool]:
    """
    Detect MMseqs2 installation and determine if WSL is needed.

    Returns:
        Tuple of (mmseqs_path, use_wsl)
        - mmseqs_path: Path to mmseqs binary or None if not found
        - use_wsl: True if WSL should be used, False otherwise
    """
    is_windows = platform.system() == 'Windows' or os.name == 'nt'

    if is_windows:
        # Case 1 & 2: Windows
        # Check common WSL paths
        wsl_paths = [
            '~/MMseqs2/build/bin/mmseqs',
            '~/mmseqs2/bin/mmseqs',
            '/usr/local/bin/mmseqs',
            '/usr/bin/mmseqs'
        ]

        # First try 'which' command
        try:
            result = subprocess.run(['wsl', 'which', 'mmseqs'],
                                    capture_output=True, text=True, shell=True)
            if result.returncode == 0:
                # Case 1: Windows, found it
                path = result.stdout.strip()
                return path, True
        except:
            pass

        # Try known paths
        for path in wsl_paths:
            try:
                cmd = f'wsl bash -c "test -f {path} && echo exists"'
                result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
                if result.stdout.strip() == 'exists':
                    # Case 1: Windows, found it
                    return path, True
            except:
                pass

        # Case 2: Windows, cannot find it
        print("ERROR: MMseqs2 not found in WSL. Please install MMseqs2 in WSL:")
        print("  wsl sudo apt install mmseqs2")
        print("  or build from source in WSL")
        return None, True

    else:
        # Case 3 & 4: Linux/WSL
        # Try which command first
        try:
            result = subprocess.run(['which', 'mmseqs'],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                # Case 3: WSL/Linux, found it
                path = result.stdout.strip()
                return path, False
        except:
            pass

        # Check common paths
        linux_paths = [
            '/usr/local/bin/mmseqs',
            '/usr/bin/mmseqs',
            os.path.expanduser('~/MMseqs2/build/bin/mmseqs'),
            os.path.expanduser('~/mmseqs2/bin/mmseqs')
        ]

        for path in linux_paths:
            if os.path.exists(path):
                # Case 3: WSL/Linux, found it
                return path, False

        # Case 4: WSL/Linux, cannot find it
        print("ERROR: MMseqs2 not found. Please install MMseqs2:")
        print("  sudo apt install mmseqs2")
        print("  or build from source")
        return None, False


def get_mmseqs_command(mmseqs_path: str, use_wsl: bool) -> list:
    """
    Get the command prefix for running MMseqs2.
    
    Args:
        mmseqs_path: Path to mmseqs binary
        use_wsl: Whether to use WSL
        
    Returns:
        List of command components
    """
    if use_wsl:
        return ['wsl', mmseqs_path]
    else:
        return [mmseqs_path]


def windows_to_wsl_path(windows_path: str) -> str:
    """
    Convert Windows path to WSL path.
    
    Examples:
        C:\\Users\\user\\file.txt -> /mnt/c/Users/user/file.txt
        C:/Users/user/file.txt -> /mnt/c/Users/user/file.txt
    
    Args:
        windows_path: Windows-style path
        
    Returns:
        WSL-compatible path
    """
    # Handle empty or None paths
    if not windows_path:
        return windows_path
    
    # Replace backslashes with forward slashes
    path = windows_path.replace('\\', '/')
    
    # Convert drive letter to WSL mount point
    # Match patterns like C:/ or c:/
    drive_match = re.match(r'^([A-Za-z]):[/\\]', path)
    if drive_match:
        drive_letter = drive_match.group(1).lower()
        # Replace C:/ with /mnt/c/
        path = re.sub(r'^[A-Za-z]:[/\\]', f'/mnt/{drive_letter}/', path)
    
    return path


def ensure_mmseqs2() -> Tuple[Optional[str], bool]:
    """
    Ensure MMseqs2 is available, with helpful error messages.
    
    Returns:
        Tuple of (mmseqs_path, use_wsl)
        
    Raises:
        RuntimeError: If MMseqs2 is not found
    """
    mmseqs_path, use_wsl = detect_mmseqs2()
    
    if not mmseqs_path:
        is_windows = platform.system() == 'Windows'
        
        if is_windows:
            raise RuntimeError(
                "MMseqs2 not found!\n"
                "On Windows, you need to:\n"
                "1. Install WSL (Windows Subsystem for Linux)\n"
                "2. In WSL, install MMseqs2:\n"
                "   - git clone https://github.com/soedinglab/MMseqs2.git\n"
                "   - cd MMseqs2 && mkdir build && cd build\n"
                "   - cmake -DCMAKE_BUILD_TYPE=RELEASE ..\n"
                "   - make -j 4\n"
                "3. Set environment variables:\n"
                "   - set MMSEQS_PATH=~/MMseqs2/build/bin/mmseqs\n"
                "   - set USE_WSL_MMSEQS=true"
            )
        else:
            raise RuntimeError(
                "MMseqs2 not found!\n"
                "Please install MMseqs2:\n"
                "- Ubuntu/Debian: sudo apt install mmseqs2\n"
                "- Or build from source:\n"
                "  git clone https://github.com/soedinglab/MMseqs2.git\n"
                "  cd MMseqs2 && mkdir build && cd build\n"
                "  cmake -DCMAKE_BUILD_TYPE=RELEASE ..\n"
                "  make && sudo make install"
            )
    
    return mmseqs_path, use_wsl
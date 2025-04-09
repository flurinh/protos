"""
Standardized GRN utility functions.

This module provides standardized utility functions for working with GRN (Generic Residue
Number) notations in protein structures. These functions follow the standard schema
definitions and provide consistent interfaces for GRN operations.
"""

import re
import warnings
from typing import Dict, List, Tuple, Union, Optional, Any, Pattern

from protos.processing.schema.schema_definitions import (
    GRN_PATTERNS,
    GRN_GAP_SYMBOL,
    GRN_UNKNOWN_SYMBOL
)


def parse_grn_str2float(grn: str) -> float:
    """
    Convert a GRN string to its float representation for numerical operations.
    
    Handles multiple formats:
    - 'n.XX' for N-terminal: converts to negative values
    - 'c.XX' for C-terminal: adds to 100
    - 'TxYY' for transmembrane regions with 'x' notation
    - 'T.YY' for transmembrane regions with dot notation
    
    Args:
        grn: GRN string to convert (e.g., '1x50', 'n.10', 'c.5', '2.45')
        
    Returns:
        Float representation of the GRN position, or 0.0 for invalid strings
        
    Examples:
        >>> parse_grn_str2float('1x50')
        1.5
        >>> parse_grn_str2float('n.10')
        -0.1
        >>> parse_grn_str2float('c.5')
        100.05
        >>> parse_grn_str2float('2.45')
        2.045i a
    """
    try:
        # N-terminal region
        if 'n.' in grn:
            # Check if the format matches expected pattern
            n_pattern = re.compile(GRN_PATTERNS['n_term'])
            if not n_pattern.match(grn):
                warnings.warn(f"Invalid N-terminal GRN format: {grn}")
                return 0.0
                
            # Parse position number
            position = int(grn.split('n.')[1])
            return -0.01 * position
            
        # C-terminal region
        elif 'c.' in grn:
            # Check if the format matches expected pattern
            c_pattern = re.compile(GRN_PATTERNS['c_term'])
            if not c_pattern.match(grn):
                warnings.warn(f"Invalid C-terminal GRN format: {grn}")
                return 0.0
                
            # Parse position number
            position = int(grn.split('c.')[1])
            return 100.0 + 0.01 * position
            
        # Transmembrane region with 'x' notation
        elif 'x' in grn:
            # Check if the format matches expected pattern
            standard_pattern = re.compile(GRN_PATTERNS['standard'])
            if not standard_pattern.match(grn):
                warnings.warn(f"Invalid standard GRN format: {grn}")
                return 0.0
                
            # Parse helix and position
            helix_str, position_str = grn.split('x')
            helix = int(helix_str)
            position = int(position_str)
            return helix + position / 100.0
            
        # Transmembrane region with dot notation
        elif '.' in grn and 'n.' not in grn and 'c.' not in grn:
            # Check if the format matches expected pattern
            loop_pattern = re.compile(GRN_PATTERNS['loop'])
            if not loop_pattern.match(grn):
                warnings.warn(f"Invalid loop GRN format: {grn}")
                return 0.0
                
            # Parse helix and position
            helix_str, position_str = grn.split('.')
            helix = int(helix_str)
            position = int(position_str)
            return helix + position / 1000.0
            
        # Invalid or unrecognized format
        else:
            warnings.warn(f"Unrecognized GRN format: {grn}")
            return 0.0
            
    except (ValueError, IndexError) as e:
        warnings.warn(f"Error parsing GRN string '{grn}': {e}")
        return 0.0


def parse_grn_float2str(grn_float: float) -> str:
    """
    Convert a GRN float representation to its standardized string format.
    
    Handles multiple formats:
    - Negative values: convert to 'n.XX' for N-terminal
    - Values 100+: convert to 'c.XX' for C-terminal
    - Values between 0 and 100: convert to standard GRN notation
    
    Args:
        grn_float: Float representation of GRN (e.g., 1.5, -0.1, 100.05, 2.045)
        
    Returns:
        Standardized GRN string representation
        
    Examples:
        >>> parse_grn_float2str(1.5)
        '1x50'
        >>> parse_grn_float2str(-0.1)
        'n.10'
        >>> parse_grn_float2str(100.05)
        'c.5'
        >>> parse_grn_float2str(2.045)
        '2.45'
    """
    # Round to 3 decimal places to avoid floating point issues
    grn_float = round(grn_float, 3)
    
    # N-terminal region (negative values)
    if grn_float < 0:
        # Convert to n.XX format
        position = int(abs(grn_float) * 100)
        return f"n.{position}"
    
    # C-terminal region (100+)
    elif grn_float >= 100:
        # Handle special case for values between 100 and 101
        if 100 <= grn_float < 101:
            # For values like 100.05, convert to c.5 format
            position = int(round((grn_float - 100) * 100))
            return f"c.{position}"
        # For values 101+, subtract 100 for the position
        else:
            position = int(grn_float - 100)
            return f"c.{position}"
    
    # Standard transmembrane region or loop
    else:
        # Check if it's a loop position (low decimal)
        decimal_part = grn_float - int(grn_float)
        
        # Use standard notation (e.g., 1x50) for normal positions
        if decimal_part >= 0.1:
            # Split into helix and position parts
            helix = int(grn_float)
            position = int(round(decimal_part * 100))
            
            # Format with proper zero padding for single-digit positions
            if position < 10:
                return f"{helix}x0{position}"
            else:
                return f"{helix}x{position}"
        # Use loop notation (e.g., 2.45) for loop positions
        else:
            # Split into helix and position parts
            helix = int(grn_float)
            position = int(round(decimal_part * 1000))
            
            # Format with proper zero padding for single-digit positions
            if position < 10:
                return f"{helix}.0{position}"
            else:
                return f"{helix}.{position}"


def normalize_grn_format(grn: str) -> str:
    """
    Normalize a GRN string to the standardized format.
    
    Converts legacy formats to the new standard:
    - '12x05' -> '12.005' (loop with x notation)
    - '12.5' -> '12.005' (loop without zero padding)
    - '1.2' -> '1x20' (standard GRN with dot instead of x)
    
    Args:
        grn: GRN string to normalize
        
    Returns:
        Normalized GRN string
    
    Examples:
        >>> normalize_grn_format('12x05')
        '12.005'
        >>> normalize_grn_format('12.5')
        '12.005'
        >>> normalize_grn_format('1x50')
        '1x50'
        >>> normalize_grn_format('1.50')
        '1x50'
    """
    # Check if already in standard format
    for pattern_name, pattern_str in GRN_PATTERNS.items():
        if re.match(pattern_str, grn):
            return grn
    
    # Legacy loop format with x (e.g., '12x05')
    loop_x_pattern = re.compile(r'^([1-8])([1-8])x(\d+)$')
    match = loop_x_pattern.match(grn)
    if match:
        helix_pair = match.group(1) + match.group(2)
        distance = int(match.group(3))
        return f"{helix_pair}.{distance:03d}"
    
    # Legacy loop format without zero padding (e.g., '12.5')
    loop_no_padding_pattern = re.compile(r'^([1-8])([1-8])\.(\d+)$')
    match = loop_no_padding_pattern.match(grn)
    if match and len(match.group(3)) < 3:
        helix_pair = match.group(1) + match.group(2)
        distance = int(match.group(3))
        return f"{helix_pair}.{distance:03d}"
    
    # Standard GRN with dot instead of x (e.g., '1.50')
    std_dot_pattern = re.compile(r'^([1-8])\.(\d+)$')
    match = std_dot_pattern.match(grn)
    if match and len(match.group(1)) == 1:
        helix = match.group(1)
        position = int(match.group(2))
        return f"{helix}x{position:02d}"
    
    # Return as is if can't normalize
    return grn


def validate_grn_string(grn: str) -> Tuple[bool, str]:
    """
    Validate a GRN string against standard patterns and return validation status.
    
    Args:
        grn: GRN string to validate
        
    Returns:
        Tuple of (is_valid, message) where:
          - is_valid: Boolean indicating if the GRN string is valid
          - message: Validation message (error message if invalid, success message if valid)
    
    Examples:
        >>> validate_grn_string('1x50')
        (True, 'Valid standard GRN format')
        >>> validate_grn_string('1x')
        (False, 'Invalid GRN format: does not match any recognized pattern')
    """
    # If empty or None, it's invalid
    if not grn:
        return False, "Empty or None GRN string"
    
    # Check against all defined patterns
    for pattern_name, pattern_str in GRN_PATTERNS.items():
        pattern = re.compile(pattern_str)
        if pattern.match(grn):
            # Now validate against more specific rules
            
            # N-terminal rules
            if pattern_name == 'n_term':
                if grn[2] == '0':  # n.01 not allowed (leading zero)
                    return False, f"Invalid N-terminal GRN format: leading zero not allowed in {grn}"
                return True, "Valid N-terminal GRN format"
                
            # C-terminal rules
            elif pattern_name == 'c_term':
                if grn[2] == '0':  # c.01 not allowed (leading zero)
                    return False, f"Invalid C-terminal GRN format: leading zero not allowed in {grn}"
                return True, "Valid C-terminal GRN format"
                
            # Standard GRN rules
            elif pattern_name == 'standard':
                # Additional validation for standard format if needed
                helix_str, position_str = grn.split('x')
                try:
                    helix = int(helix_str)
                    position = int(position_str)
                    
                    # Check helix range (typically 1-8 for GPCRs)
                    if not (1 <= helix <= 8):
                        return False, f"Invalid helix number: {helix} (expected 1-8)"
                        
                    # Check position (typically 1-99)
                    if not (1 <= position <= 99):
                        return False, f"Invalid position number: {position} (expected 1-99)"
                except ValueError:
                    return False, f"Non-numeric values in GRN: {grn}"
                
                return True, "Valid standard GRN format"
                
            # Loop GRN rules
            elif pattern_name == 'loop':
                # Additional validation for loop format if needed
                helix_str, position_str = grn.split('.')
                try:
                    helix = int(helix_str)
                    position = int(position_str)
                    
                    # Check helix range (typically 1-8 for GPCRs)
                    if not (1 <= helix <= 8):
                        return False, f"Invalid helix number: {helix} (expected 1-8)"
                        
                    # Check position (typically 1-999)
                    if not (1 <= position <= 999):
                        return False, f"Invalid position number: {position} (expected 1-999)"
                except ValueError:
                    return False, f"Non-numeric values in GRN: {grn}"
                
                return True, "Valid loop GRN format"
    
    # If we got here, no pattern matched
    return False, "Invalid GRN format: does not match any recognized pattern"


def sort_grns(grn_list: List[Union[str, float]]) -> List[Union[str, float]]:
    """
    Sort a list of GRNs in standard order (N-terminal, TM helices, loops, C-terminal).
    
    This function accepts either string or float representations of GRNs and
    returns a sorted list of the same type.
    
    Args:
        grn_list: List of GRNs to sort (either strings or floats)
        
    Returns:
        Sorted list of GRNs in the same format as the input
        
    Examples:
        >>> sort_grns(['3x50', '1x50', 'n.10', 'c.5'])
        ['n.10', '1x50', '3x50', 'c.5']
        >>> sort_grns([3.5, 1.5, -0.1, 100.05])
        [-0.1, 1.5, 3.5, 100.05]
    """
    # Determine input type (string or float)
    input_type = str if all(isinstance(g, str) for g in grn_list) else float
    
    # Convert all to float representation for sorting
    if input_type == str:
        # Convert strings to floats
        float_values = [parse_grn_str2float(g) for g in grn_list]
    else:
        # Already floats
        float_values = grn_list
    
    # Create a mapping to preserve original values
    value_map = dict(zip(float_values, grn_list))
    
    # Sort the float values
    sorted_floats = sorted(float_values)
    
    # Convert back to original format
    return [value_map[f] for f in sorted_floats]
    

def get_grn_region(grn: Union[str, float]) -> str:
    """
    Determine the region of a protein structure that a GRN belongs to.
    
    Args:
        grn: GRN in either string or float representation
        
    Returns:
        Region name: 'N-terminal', 'TM1' through 'TM8', 'Loop', or 'C-terminal'
        
    Examples:
        >>> get_grn_region('1x50')
        'TM1'
        >>> get_grn_region('n.10')
        'N-terminal'
        >>> get_grn_region('c.5')
        'C-terminal'
        >>> get_grn_region('2.45')
        'Loop'
    """
    # Convert to float if string
    if isinstance(grn, str):
        grn_float = parse_grn_str2float(grn)
    else:
        grn_float = grn
    
    # N-terminal region
    if grn_float < 0:
        return 'N-terminal'
    
    # C-terminal region
    elif grn_float >= 100:
        return 'C-terminal'
    
    # TM or Loop region
    else:
        # Get the helix number and decimal part
        helix = int(grn_float)
        decimal = grn_float - helix
        
        # If decimal is small, it's a loop position
        if decimal < 0.1:
            return 'Loop'
        # Otherwise it's a TM helix
        else:
            return f'TM{helix}'
"""
Updated GRN utility functions with corrected loop formatting.

This module provides standardized utility functions for working with GRN (Generic Residue
Number) notations in protein structures. These functions follow the standard schema
definitions and provide consistent interfaces for GRN operations.

The key improvement in this update is the correct handling of loop regions:
- Format: <closer helix><further helix>.<distance>
- Example: 12.003 (loop between helix 1-2, closer to 1, distance 3)
"""

import re
import warnings
import logging
from typing import Dict, List, Tuple, Union, Optional, Any, Pattern

from protos.processing.schema.schema_definitions import (
    GRN_PATTERNS,
    GRN_FORMAT_DOCS,
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
    - 'AB.CCC' for loop regions (between helix A and B, closer to A, distance CCC)
    
    Args:
        grn: GRN string to convert (e.g., '1x50', 'n.10', 'c.5', '12.003')
        
    Returns:
        Float representation of the GRN position, or 0.0 for invalid strings
        
    Examples:
        >>> parse_grn_str2float('1x50')
        1.5
        >>> parse_grn_str2float('n.10')
        -0.1
        >>> parse_grn_str2float('c.5')
        100.05
        >>> parse_grn_str2float('12.003')
        12.003
    """
    try:
        # N-terminal region
        if 'n.' in grn:
            # Parse position number
            position = int(grn.split('n.')[1])
            return -0.01 * position
            
        # C-terminal region
        elif 'c.' in grn:
            # Parse position number
            position = int(grn.split('c.')[1])
            return 100.0 + 0.01 * position
            
        # Loop region with format AB.CCC
        elif '.' in grn and len(grn.split('.')[1]) == 3:
            # Parse helix pair and distance
            helix_pair = grn.split('.')[0]
            distance = int(grn.split('.')[1]) / 1000.0
            
            if len(helix_pair) == 2:
                # Extract closer and further helix
                closer_helix = int(helix_pair[0])
                further_helix = int(helix_pair[1])
                
                # Create float representation:
                # Integer part: 10*smaller_helix + larger_helix
                # Decimal part: normalized distance (0-1)
                helix_min = min(closer_helix, further_helix)
                helix_max = max(closer_helix, further_helix)
                
                # Calculate the float representation
                return float(f"{helix_min}{helix_max}") + distance
            else:
                raise ValueError(f"Invalid loop format: {grn}, expected format: AB.CCC")
                
        # Standard GRN format (TM regions)
        elif 'x' in grn:
            # Parse helix and position
            helix_str, position_str = grn.split('x')
            helix = int(helix_str)
            position = int(position_str)
            return helix + position / 100.0
            
        # Invalid or unrecognized format
        else:
            raise ValueError(f"Unrecognized GRN format: {grn}")
            
    except (ValueError, IndexError) as e:
        # Log the error with the GRN string for debugging
        logging.error(f"Error parsing GRN string '{grn}': {e}")
        return 0.0


def parse_grn_float2str(grn_float: float) -> str:
    """
    Convert a GRN float representation to its string format.
    
    Handles:
    - Standard: 1.50 -> '1x50'
    - N-terminal: -0.10 -> 'n.10'
    - C-terminal: 100.05 -> 'c.5'
    - Loop: 12.003 -> '12.003' (loop between helix 1-2, closer to 1, distance 3)
    
    Args:
        grn_float: Float representation of GRN (e.g., 1.5, -0.1, 100.05, 12.003)
        
    Returns:
        Standardized GRN string representation
        
    Examples:
        >>> parse_grn_float2str(1.5)
        '1x50'
        >>> parse_grn_float2str(-0.1)
        'n.10'
        >>> parse_grn_float2str(100.05)
        'c.5'
        >>> parse_grn_float2str(12.003)
        '12.003'
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
        # For values like 100.05, convert to c.5 format
        position = int(round((grn_float - 100) * 100))
        return f"c.{position}"
    
    # Loop region (values between 10 and 100)
    elif grn_float >= 10:
        # Extract the parts
        int_part = int(grn_float)
        decimal_part = round((grn_float - int_part) * 1000)
        
        # Get helix numbers
        helix1 = int(int_part // 10)  # First digit
        helix2 = int(int_part % 10)   # Second digit
        
        # Format with proper zero padding for the distance
        return f"{helix1}{helix2}.{decimal_part:03d}"
    
    # Standard transmembrane region
    else:
        # Split into helix and position parts
        helix = int(grn_float)
        position = int(round((grn_float - helix) * 100))
        
        # Format with proper zero padding
        return f"{helix}x{position:02d}"


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
        >>> validate_grn_string('12.003')
        (True, 'Valid loop GRN format')
        >>> validate_grn_string('9x50')
        (False, 'Invalid helix number: 9 (expected 1-8)')
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
                # Parse the components
                helix_pair = grn.split('.')[0]
                distance_str = grn.split('.')[1]
                
                try:
                    # Validate helices
                    closer_helix = int(helix_pair[0])
                    further_helix = int(helix_pair[1])
                    
                    # Validate helix values
                    if not (1 <= closer_helix <= 8):
                        return False, f"Invalid closer helix: {closer_helix} (expected 1-8)"
                    if not (1 <= further_helix <= 8):
                        return False, f"Invalid further helix: {further_helix} (expected 1-8)"
                        
                    # Validate distance value
                    distance = int(distance_str)
                    if not (0 <= distance <= 999):
                        return False, f"Invalid distance: {distance} (expected 0-999)"
                        
                    return True, "Valid loop GRN format"
                except ValueError:
                    return False, f"Non-numeric values in loop GRN: {grn}"
    
    # Try to normalize and validate again
    normalized = normalize_grn_format(grn)
    if normalized != grn:
        is_valid, message = validate_grn_string(normalized)
        if is_valid:
            return True, f"Valid GRN after normalization: {grn} -> {normalized}"
    
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
        >>> sort_grns(['3x50', '1x50', 'n.10', 'c.5', '12.003'])
        ['n.10', '1x50', '3x50', '12.003', 'c.5']
        >>> sort_grns([3.5, 1.5, -0.1, 100.05, 12.003])
        [-0.1, 1.5, 3.5, 12.003, 100.05]
    """
    # Handle empty list
    if not grn_list:
        return []
        
    # Determine input type (string or float)
    input_type = type(grn_list[0])
    
    # Convert all to float representation for sorting
    if input_type == str:
        # Convert strings to floats
        try:
            float_values = [parse_grn_str2float(g) for g in grn_list]
            # Create a mapping to preserve original values
            value_map = dict(zip(float_values, grn_list))
        except Exception as e:
            # If conversion fails, just sort alphabetically
            logging.warning(f"Error converting GRNs to float for sorting: {e}")
            return sorted(grn_list)
    else:
        # Already floats
        float_values = grn_list
        value_map = dict(zip(float_values, float_values))
    
    # Sort the float values
    sorted_floats = sorted(float_values)
    
    # Convert back to original type
    if input_type == str:
        return [value_map[f] for f in sorted_floats]
    else:
        return sorted_floats
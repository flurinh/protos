import re
import os
import json
import pandas as pd
import logging
from decimal import Decimal, getcontext
from typing import List, Tuple, Dict, Optional, Union, Any
from pathlib import Path

ERROR_FLOAT = -999.999  # Unique float to indicate parsing error

# Configure logger
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
#   GRN Pattern Definitions (from schema_definitions.py)
# ---------------------------------------------------------------------------

# GRN format specifications
GRN_PATTERNS = {
    'standard': r'^(\d+)x(\d+)$',  # e.g., 1x50
    'standard_dot': r'^([0-9])\.(\d+)$',  # e.g., 1.50, 0.00, 9.99 (dot notation)
    'n_term': r'^n\.(\d+)$',       # e.g., n.10
    'c_term': r'^c\.(\d+)$',       # e.g., c.5
    'loop': r'^([1-8])([1-8])\.(\d+)$'  # e.g., 12.003, 65.011, also 12.47
}

# Documentation for GRN formats
GRN_FORMAT_DOCS = {
    'standard': "Standard GRN format: <helix>x<position> (e.g., 1x50)",
    'standard_dot': "Standard GRN format with dot notation: <helix>.<position> (e.g., 1.50)",
    'n_term': "N-terminal format: n.<position> (e.g., n.10)",
    'c_term': "C-terminal format: c.<position> (e.g., c.5)",
    'loop': """Loop region format: <closer helix><further helix>.<distance> where:
            - First digit: Closer helix (1-8)
            - Second digit: Further helix (1-8)
            - Three-digit decimal: Distance from closer helix (001-999)
            Examples: 12.003 (between helix 1-2, closer to 1, distance 3)
                     65.011 (between helix 5-6, closer to 6, distance 11)"""
}

# Symbol definitions
GRN_GAP_SYMBOL = '-'
GRN_UNKNOWN_SYMBOL = 'X'

# ---------------------------------------------------------------------------
#   1. String to Float Parser (Updated version with improved handling)
# ---------------------------------------------------------------------------
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
            
        # Check for dot notation (need to distinguish between loops and standard dot notation)
        elif '.' in grn:
            helix_part = grn.split('.')[0]
            position_part = grn.split('.')[1]
            
            # Loop region with format AB.CCC (2 digits before dot, 3 after)
            if len(helix_part) == 2 and len(position_part) == 3:
                # Parse helix pair and distance
                distance = int(position_part) / 1000.0
                
                # Extract closer and further helix
                closer_helix = int(helix_part[0])
                further_helix = int(helix_part[1])
                
                # Create float representation:
                # Integer part: 10*smaller_helix + larger_helix
                # Decimal part: normalized distance (0-1)
                helix_min = min(closer_helix, further_helix)
                helix_max = max(closer_helix, further_helix)
                
                # Calculate the float representation
                return float(f"{helix_min}{helix_max}") + distance
            
            # Standard dot notation (1 digit before dot, 1-2 after) e.g., 1.50
            elif len(helix_part) == 1:
                helix = int(helix_part)
                position = int(position_part)
                return helix + position / 100.0
            
            else:
                raise ValueError(f"Invalid dot notation format: {grn}")
                
        # Standard GRN format with x notation (TM regions)
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
        logger.error(f"Error parsing GRN string '{grn}': {e}")
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
        '1.50'
        >>> normalize_grn_format('1.50')
        '1.50'
        >>> normalize_grn_format('1.5')
        '1.50'
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
    
    # Legacy loop format without zero padding (e.g., '12.5', '12.47')
    loop_no_padding_pattern = re.compile(r'^([0-9])([0-9])\.(\d+)$')
    match = loop_no_padding_pattern.match(grn)
    if match:
        helix_pair = match.group(1) + match.group(2)
        distance_str = match.group(3)
        
        # If already 3 digits, keep as is
        if len(distance_str) == 3:
            return grn
        # Otherwise pad to 3 digits
        else:
            distance = int(distance_str)
            return f"{helix_pair}.{distance:03d}"
    
    # Standard GRN with x notation (e.g., '1x50')
    std_x_pattern = re.compile(r'^([1-8])x(\d+)$')
    match = std_x_pattern.match(grn)
    if match:
        helix = match.group(1)
        position = int(match.group(2))
        return f"{helix}.{position:02d}"
    
    # Standard GRN with dot notation - normalize to 2-digit format (e.g., '1.5' -> '1.50')
    std_dot_pattern = re.compile(r'^([1-8])\.(\d+)$')
    match = std_dot_pattern.match(grn)
    if match and len(match.group(1)) == 1:
        helix = match.group(1)
        position = int(match.group(2))
        # Normalize to 2-digit format (1.5 -> 1.50, 1.50 -> 1.50)
        return f"{helix}.{position:02d}"
    
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
                
            # Standard dot notation rules
            elif pattern_name == 'standard_dot':
                # Additional validation for standard dot format
                helix_str, position_str = grn.split('.')
                try:
                    helix = int(helix_str)
                    position = int(position_str)
                    
                    # Check helix range
                    if not (0 <= helix <= 9):
                        return False, f"Invalid helix number: {helix} (expected 0-9)"
                        
                    # Check position
                    if not (0 <= position <= 99):
                        return False, f"Invalid position number: {position} (expected 0-99)"
                except ValueError:
                    return False, f"Non-numeric values in GRN: {grn}"
                
                return True, "Valid standard dot notation GRN format"
                
            # Loop region rules  
            elif pattern_name == 'loop':
                # Additional validation for loop format
                match = pattern.match(grn)
                helix1, helix2, distance = match.groups()
                try:
                    helix1_int = int(helix1)
                    helix2_int = int(helix2)
                    distance_int = int(distance)
                    
                    # Check helix range
                    if not (1 <= helix1_int <= 8) or not (1 <= helix2_int <= 8):
                        return False, f"Invalid helix numbers in loop: {helix1}, {helix2} (expected 1-8)"
                        
                    # Helices should be adjacent or at least make sense
                    if abs(helix1_int - helix2_int) > 1 and not (helix1_int == 1 and helix2_int == 8):
                        return False, f"Non-adjacent helices in loop: {helix1}, {helix2}"
                        
                    # Check distance range (typically 1-999)
                    if not (1 <= distance_int <= 999):
                        return False, f"Invalid distance in loop: {distance} (expected 1-999)"
                except ValueError:
                    return False, f"Non-numeric values in loop GRN: {grn}"
                
                return True, "Valid loop GRN format"
    
    # If no pattern matched
    return False, f"GRN string '{grn}' does not match any known pattern"


# ---------------------------------------------------------------------------
#   2. Legacy check function (kept for compatibility)
# ---------------------------------------------------------------------------
def check_str_grn_valid(grn_str):
    float_val = parse_grn_str2float(grn_str)
    return float_val != ERROR_FLOAT


# ---------------------------------------------------------------------------
#   3. Helper Functions
# ---------------------------------------------------------------------------
def get_prev_next_tm(loop_grn_float):
    """
    For a given loop GRN float, extract the previous and next TM helices.
    
    Loop format: AB.CCC where A and B are helix numbers
    """
    if loop_grn_float < 10:  # Not a loop
        return None, None
        
    int_part = int(loop_grn_float)
    prev_tm = int_part // 10
    next_tm = int_part % 10
    
    return prev_tm, next_tm


# ---------------------------------------------------------------------------
#   4. Sorting Functions
# ---------------------------------------------------------------------------
def sort_grns(grn_floats: List[float]) -> List[float]:
    """
    Sort a list of GRN floats in the standard order.
    
    Order:
    1. N-terminal (negative values)
    2. TM regions (1-8)
    3. Loops (10-99)
    4. C-terminal (100+)
    
    Within each category, sort by value.
    """
    # Separate into categories
    n_term = [g for g in grn_floats if g < 0]
    tm_regions = [g for g in grn_floats if 0 <= g < 10]
    loops = [g for g in grn_floats if 10 <= g < 100]
    c_term = [g for g in grn_floats if g >= 100]
    
    # Sort each category
    n_term.sort(reverse=True)  # Most negative last (closer to TM1)
    tm_regions.sort()
    loops.sort()
    c_term.sort()
    
    # Combine in order
    return n_term + tm_regions + loops + c_term


def sort_grns_str(grn_strs: List[str]) -> List[str]:
    """Sort a list of GRN strings."""
    # Convert to floats, sort, convert back
    grn_floats = [parse_grn_str2float(g) for g in grn_strs]
    sorted_floats = sort_grns(grn_floats)
    return [parse_grn_float2str(f) for f in sorted_floats]


# ---------------------------------------------------------------------------
#   5. GRN Interval and Configuration Management
# ---------------------------------------------------------------------------
def init_grn_intervals(protein_family, min_seq_id=0.3, max_gaps=20, max_e_value=10e-5, config_dir=None):
    """Initialize GRN intervals from configuration files."""
    if config_dir is None:
        # Use default config directory
        config_dir = Path(__file__).parent / 'configs'
    else:
        config_dir = Path(config_dir)
        
    config_file = config_dir / f'{protein_family}.json'
    
    if not config_file.exists():
        logger.warning(f"Config file not found: {config_file}")
        return {}
        
    with open(config_file, 'r') as f:
        config_data = json.load(f)
        
    intervals = {}
    for grn_range in config_data.get('grn_ranges', []):
        start_grn = grn_range['start']
        end_grn = grn_range['end']
        intervals[f"{start_grn}-{end_grn}"] = (
            parse_grn_str2float(start_grn),
            parse_grn_str2float(end_grn)
        )
    
    return intervals


class GRNConfigManager:
    """Manages GRN configurations for different protein families."""
    
    def __init__(self, config_dir=None):
        if config_dir is None:
            self.config_dir = Path(__file__).parent / 'configs'
        else:
            self.config_dir = Path(config_dir)
            
        self.configs = {}
        self._load_configs()
    
    def _load_configs(self):
        """Load all available configuration files."""
        if not self.config_dir.exists():
            logger.warning(f"Config directory not found: {self.config_dir}")
            return
            
        for config_file in self.config_dir.glob('*.json'):
            family_name = config_file.stem
            try:
                with open(config_file, 'r') as f:
                    self.configs[family_name] = json.load(f)
            except Exception as e:
                logger.error(f"Error loading config {config_file}: {e}")
    
    def get_intervals(self, protein_family):
        """Get GRN intervals for a specific protein family."""
        if protein_family not in self.configs:
            logger.warning(f"No config found for protein family: {protein_family}")
            return {}
            
        return init_grn_intervals(protein_family, config_dir=self.config_dir)


# ---------------------------------------------------------------------------
#   6. GRN Interval Calculation
# ---------------------------------------------------------------------------
def get_grn_interval(grn_target, intervals):
    """
    Get the GRN interval that contains the target GRN.
    
    Args:
        grn_target: Target GRN (string or float)
        intervals: Dictionary of interval_name -> (start_float, end_float)
        
    Returns:
        Interval name or None if not found
    """
    if isinstance(grn_target, str):
        grn_target_float = parse_grn_str2float(grn_target)
    else:
        grn_target_float = grn_target
        
    for interval_name, (start, end) in intervals.items():
        if start <= grn_target_float <= end:
            return interval_name
            
    return None


def init_std_grns(intervals):
    """Initialize standard GRNs from intervals."""
    std_grns = []
    for interval_name, (start, end) in intervals.items():
        # Add the start and end of each interval
        std_grns.append(parse_grn_float2str(start))
        std_grns.append(parse_grn_float2str(end))
    
    # Remove duplicates and sort
    std_grns = list(set(std_grns))
    return sort_grns_str(std_grns)


# ---------------------------------------------------------------------------
#   7. Visualization Support
# ---------------------------------------------------------------------------
def map_grn_to_color(grn, color_map=None):
    """Map a GRN to a color for visualization."""
    if color_map is None:
        # Default color scheme based on helix
        color_map = {
            1: '#FF0000',  # Red
            2: '#FF7F00',  # Orange
            3: '#FFFF00',  # Yellow
            4: '#00FF00',  # Green
            5: '#0000FF',  # Blue
            6: '#4B0082',  # Indigo
            7: '#9400D3',  # Violet
            8: '#FF1493',  # Deep Pink
        }
    
    grn_float = parse_grn_str2float(grn)
    
    # N-terminal: Gray
    if grn_float < 0:
        return '#808080'
    # C-terminal: Black
    elif grn_float >= 100:
        return '#000000'
    # Loops: Light gray
    elif grn_float >= 10:
        return '#C0C0C0'
    # TM regions: Use helix color
    else:
        helix = int(grn_float)
        return color_map.get(helix, '#FFFFFF')


# ---------------------------------------------------------------------------
#   8. Sequence Utilities
# ---------------------------------------------------------------------------
def get_seq(name, grn_table):
    """Get sequence from GRN table for a specific protein."""
    if name not in grn_table.index:
        return None
        
    row = grn_table.loc[name]
    # Filter out gaps and unknowns
    seq_parts = [val[0] for val in row.values if val not in [GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL, '-', 'X']]
    return ''.join(seq_parts)


def get_annot_seq(name, grn_table):
    """Get annotated sequence with residue numbers from GRN table."""
    if name not in grn_table.index:
        return None
        
    row = grn_table.loc[name]
    # Include residue position information
    annot_parts = []
    for grn, val in row.items():
        if val not in [GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL, '-', 'X']:
            annot_parts.append(f"{val}({grn})")
    
    return ''.join(annot_parts)


def remove_gaps_from_sequences(seq_dict):
    """Remove gaps from sequences in a dictionary."""
    cleaned_dict = {}
    for name, seq in seq_dict.items():
        cleaned_seq = seq.replace('-', '').replace('.', '')
        cleaned_dict[name] = cleaned_seq
    return cleaned_dict


# ---------------------------------------------------------------------------
#   9. Utility Functions
# ---------------------------------------------------------------------------
def flatten(lst):
    """Flatten a nested list."""
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result


def get_tm_residues(grn_list):
    """Filter TM residues from a GRN list."""
    tm_residues = []
    for grn in grn_list:
        grn_float = parse_grn_str2float(grn)
        # TM regions are between 1 and 8 (not loops)
        if 1 <= grn_float < 10:
            tm_residues.append(grn)
    return tm_residues


# ---------------------------------------------------------------------------
#   10. Testing Function
# ---------------------------------------------------------------------------
def run_all_tests():
    """Run all test cases for GRN utilities."""
    test_cases = [
        # Standard formats
        ('1x50', 1.50),
        ('7x53', 7.53),
        ('3x25', 3.25),
        
        # N-terminal
        ('n.10', -0.10),
        ('n.5', -0.05),
        ('n.25', -0.25),
        
        # C-terminal
        ('c.5', 100.05),
        ('c.10', 100.10),
        ('c.25', 100.25),
        
        # Loops
        ('12.003', 12.003),
        ('23.015', 23.015),
        ('67.100', 67.100),
    ]
    
    print("Testing parse_grn_str2float and parse_grn_float2str...")
    for grn_str, expected_float in test_cases:
        # Test string to float
        result_float = parse_grn_str2float(grn_str)
        assert abs(result_float - expected_float) < 0.0001, f"Failed: {grn_str} -> {result_float} (expected {expected_float})"
        
        # Test float to string
        result_str = parse_grn_float2str(expected_float)
        # Normalize both for comparison
        normalized_input = normalize_grn_format(grn_str)
        normalized_result = normalize_grn_format(result_str)
        assert normalized_input == normalized_result, f"Failed: {expected_float} -> {result_str} (expected {grn_str})"
    
    print("All tests passed!")


if __name__ == "__main__":
    run_all_tests()
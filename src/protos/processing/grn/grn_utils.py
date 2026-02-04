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
    'standard_insertion': r'^([0-9])\.(\d+)\.(\d+)$',
    'n_term': r'^n\.(\d+)$',       # e.g., n.10
    'c_term': r'^c\.(\d+)$',       # e.g., c.5
    'loop': r'^([1-8])([1-8])\.(\d+)$'  # e.g., 12.003, 65.011, also 12.047
}

# Documentation for GRN formats
GRN_FORMAT_DOCS = {
    'standard': "Standard GRN format: <helix>x<position> (e.g., 1x50)",
    'standard_dot': "Standard GRN format with dot notation: <helix>.<position> (e.g., 1.50)",
    'standard_insertion': "Standard GRN insertion format: <helix>x<position> (e.g., 1x521)",
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
    - 'n.XX' for N-terminal: converts to negative values (n.10 -> -10.0, n.1 -> -1.0)
    - 'c.XX' for C-terminal: adds to 100
    - 'TxYY' or 'T.YY' for transmembrane regions
    - 'AB.CCC' for loop regions where:
      - A is the closer helix (the residue is closer to this helix)
      - B is the further helix
      - CCC is the distance from the closer helix
      Example: 67.001 = closer to helix 6, between 6 and 7, distance 1
               76.001 = closer to helix 7, between 6 and 7, distance 1
    
    Args:
        grn: GRN string to convert
        
    Returns:
        Float representation of the GRN position, or 0.0 for invalid strings
        
    Examples:
        >>> parse_grn_str2float('1x50')
        1.5
        >>> parse_grn_str2float('n.10')
        -10.0
        >>> parse_grn_str2float('n.1')
        -1.0
        >>> parse_grn_str2float('c.5')
        100.05
        >>> parse_grn_str2float('67.001')  # Closer to 6, between 6-7
        67.001
        >>> parse_grn_str2float('76.001')  # Closer to 7, between 6-7
        76.001
        >>> parse_grn_str2float('76.01')  # Closer to 7, between 6-7
        76.010
        >>> parse_grn_str2float('76.1')  # Closer to 7, between 6-7
        76.100
        >>> parse_grn_str2float('76.100')  # Closer to 7, between 6-7
        76.100
    """
    try:
        # N-terminal region
        if grn.startswith('n.'):
            # Parse position number - keep as negative for sorting
            position = int(grn.split('n.')[1])
            return -1.0 * position
            
        # C-terminal region
        elif grn.startswith('c.'):
            # Parse position number
            position = int(grn.split('c.')[1])
            return 100.0 + position

        else:
            helix_part = grn.split('.')[0]
            position_part = grn.split('.')[1]
            
            # Loop region: 2 digits before dot, typically 3 after (e.g., 67.001, 76.001)
            if len(helix_part) == 2 and helix_part.isdigit():
                # Extract closer and further helix
                closer_helix = int(helix_part[0])
                further_helix = int(helix_part[1])
                
                # Parse distance, ensuring it's treated as a decimal
                position_part.ljust(3, '0')
                distance = int(position_part) / 1000.0
                
                # For sorting, we need to preserve the exact notation
                # 67.001 should sort differently from 76.001
                # Return the literal float value
                return round(float(grn), 3)
            
            # Standard TM position e.g., 1.50, 7.50
            elif len(helix_part) == 1 and helix_part.isdigit():
                helix = int(helix_part)
                position = int(position_part)
                return round(helix + position / 100.0, 3)
            
            # Otherwise treat as a literal float for special cases
            else:
                return round(float(grn), 3)

    except (ValueError, IndexError) as e:
        # Log the error with the GRN string for debugging
        logger.error(f"Error parsing GRN string '{grn}': {e}")
        return 0.0


def parse_grn_float2str(grn_float: float, notation_type: str = 'dot') -> str:
    """
    Convert a GRN float representation to its string format.
    X notation is deprecated.

    Handles:
    - Standard dot: 1.50 -> '1.50' (default)
    - Standard x: 1.50 -> '1x50' (when notation_type='x')
    - N-terminal: -10.0 -> 'n.10', -1.0 -> 'n.1'
    - C-terminal: 100.05 -> 'c.5'
    - Loop: 67.001 -> '67.001' (closer to helix 6, between 6-7)
            76.001 -> '76.001' (closer to helix 7, between 6-7)
    
    Args:
        grn_float: Float representation of GRN
        notation_type: 'dot' or 'x' for standard notation (default: 'dot')
        
    Returns:
        Standardized GRN string representation
    """
    # Round to 3 decimal places to avoid floating point issues
    grn_float = round(grn_float, 3)
    
    # N-terminal region (negative values)
    if grn_float < 0:
        # Convert to n.XX format
        position = int(abs(grn_float))
        return f"n.{position}"
    
    # C-terminal region (100+)
    elif grn_float >= 100:
        # For values like 101.0, convert to c.1 format
        position = int(round(grn_float - 100))
        return f"c.{position}"
    
    # Loop region (values between 10 and 100)
    elif grn_float >= 10:
        # For loop regions, preserve the exact format
        # This handles cases like 67.001, 76.001
        int_part = int(grn_float)
        decimal_part = grn_float - int_part
        
        # Convert decimal part to 3-digit distance
        distance = int(round(decimal_part * 1000))
        
        # Format with proper zero padding for the distance
        return f"{int_part}.{distance:03d}"
    
    # Standard transmembrane region (0-10)
    else:
        # Split into helix and position parts
        helix = int(grn_float)
        position = int(round((grn_float - helix) * 100))
        
        # Format based on notation type
        return f"{helix}.{position:02d}"


def normalize_grn_format(grn: str) -> str:
    """
    Normalize a GRN string to the standardized format.

    Converts legacy formats to the new standard:
    - '12x05' -> '12.050' (loop with x notation)
    - '12.5' -> '12.500' (loop without zero padding)
    - '1x50' -> '1.50' (standard GRN with x notation)
    - '1.2' -> '1.20' (standard GRN with dot)
    - '1.511' -> '1.511' (insertion with additional zero padding)

    Args:
        grn: GRN string to normalize

    Returns:
        Normalized GRN string

    Examples:
        >>> normalize_grn_format('12.5')
        '12.500'
        >>> normalize_grn_format('12x05')
        '12.050'
        >>> normalize_grn_format('12x005')
        '12.005'
        >>> normalize_grn_format('1x50')
        '1.50'
        >>> normalize_grn_format('1x5')
        '1.50'
        >>> normalize_grn_format('1.50')
        '1.50'
        >>> normalize_grn_format('1.5')
        '1.50'
        >>> normalize_grn_format('1.521') # insertion notation
        '1.521'
    """
    # First, replace 'x' with '.' to standardize notation
    grn = str(grn).replace('x', '.')

    # N-terminal (no normalization needed)
    if grn.startswith('n.'):
        return grn

    # C-terminal (no normalization needed)
    elif grn.startswith('c.'):
        return grn

    # Loop helix identifiers (inter-helix loops)
    loop_helices = {'12', '23', '34', '45', '54', '56', '67', '76', '78'}

    # Loop format with dot notation (e.g., '12.5' -> '12.500', '12.05' -> '12.050')
    loop_dot_pattern = re.compile(r'^([1-8])([1-8])\.(\d+)$')
    match = loop_dot_pattern.match(grn)
    if match:
        helix_pair = match.group(1) + match.group(2)
        distance_str = match.group(3)

        # Loop positions MUST have exactly 3 digits
        if len(distance_str) == 3:
            return grn
        elif len(distance_str) == 1:
            # Single digit: multiply by 100 (e.g., '5' -> '500')
            distance = int(distance_str) * 100
        elif len(distance_str) == 2:
            # Two digits: multiply by 10 (e.g., '01' -> '010', '10' -> '100')
            distance = int(distance_str) * 10
        else:
            # More than 3 digits - keep as-is (unusual case)
            return grn

        return f"{helix_pair}.{distance:03d}"

    # Standard GRN with dot notation - normalize to 2-digit format (e.g., '1.5' -> '1.50')
    std_dot_pattern = re.compile(r'^([0-9])\.(\d+)$')
    match = std_dot_pattern.match(grn)
    if match and len(match.group(1)) == 1:
        helix = match.group(1)
        position_str = match.group(2)
        position = int(position_str)

        # Handle single digit positions that should be x0 (e.g., '5' -> '50')
        if len(position_str) == 1:
            position = position * 10

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
    1. N-terminal (negative values, sorted descending: n.10, n.9, ..., n.1)
    2. Helix 1 and positions (1.xx)
    3. Loop 1-2 (12.xxx, 21.xxx)
    4. Helix 2 and positions (2.xx)
    5. Loop 2-3 (23.xxx, 32.xxx)
    ... and so on
    8. C-terminal (100+, sorted ascending: c.1, c.2, ...)
    
    Loop positions are placed between their respective helices.
    For example: 67.001, 76.001 come between helix 6 and helix 7.
    """
    # Separate into categories
    n_term = [g for g in grn_floats if g < 0]
    c_term = [g for g in grn_floats if g >= 100]
    
    # Sort N-terminal and C-terminal
    n_term.sort()  # Most negative first (-10 before -1)
    c_term.sort()  # c.1 before c.2
    
    # Process helices and loops together
    result = n_term.copy()
    
    # Process each helix and its following loop
    for helix_num in range(1, 9):  # Helices 1-8
        # Add positions from this helix (e.g., 1.39, 1.50, 1.57)
        helix_positions = [g for g in grn_floats if helix_num <= g < helix_num + 1]
        helix_positions.sort()
        result.extend(helix_positions)
        
        # Add loop positions after this helix
        if helix_num < 8:  # No loop after helix 8
            next_helix = helix_num + 1
            
            # Loop positions closer to current helix (e.g., 12.xxx for loop 1-2)
            loop_closer_current = [g for g in grn_floats 
                                 if float(f"{helix_num}{next_helix}") <= g < float(f"{helix_num}{next_helix}") + 1]
            
            # Loop positions closer to next helix (e.g., 21.xxx for loop 1-2)
            loop_closer_next = [g for g in grn_floats 
                              if float(f"{next_helix}{helix_num}") <= g < float(f"{next_helix}{helix_num}") + 1]
            
            # Sort loop positions
            loop_closer_current.sort()      # 12.001 before 12.002
            loop_closer_next.sort(reverse=True)  # 21.002 before 21.001 (reverse for closer to next)
            
            # Add loop positions: first those closer to current helix, then those closer to next
            result.extend(loop_closer_current)
            result.extend(loop_closer_next)
    
    # Add C-terminal
    result.extend(c_term)
    
    return result


def sort_grns_str(grn_strs: List[str], notation_type: str = None) -> List[str]:
    """Sort a list of GRN strings with proper helix-loop-helix ordering.
    
    Args:
        grn_strs: List of GRN strings to sort
        notation_type: If specified, convert output to this notation ('dot' or 'x').
                      If None, preserve original notation.
    """
    # Convert to floats for sorting
    grn_floats = [parse_grn_str2float(g) for g in grn_strs]
    
    # Sort using our custom sort function
    sorted_floats = sort_grns(grn_floats)
    
    # Create a mapping from float to original string for notation preservation
    float_to_orig = {parse_grn_str2float(g): g for g in grn_strs}
    
    # Convert back to strings
    if notation_type is None:
        # Preserve original notation
        sorted_strs = []
        for f in sorted_floats:
            orig = float_to_orig.get(f, '')
            if orig:
                sorted_strs.append(str(orig))
            else:
                # Fallback if original not found
                sorted_strs.append(parse_grn_float2str(f, notation_type='dot'))
        return sorted_strs
    else:
        # Convert all to specified notation
        return [parse_grn_float2str(f, notation_type=notation_type) for f in sorted_floats]


# ---------------------------------------------------------------------------
#   5. GRN Interval and Configuration Management
# ---------------------------------------------------------------------------
def get_grn_interval(start_grn: str, end_grn: str, grns_str: List[str] = None) -> List[str]:
    """Get GRN interval between start and end.
    
    Args:
        start_grn: Starting GRN position
        end_grn: Ending GRN position
        grns_str: Optional list of valid GRN strings to filter from.
                 If None, auto-generates GRNs with 0.01 stepsize.
        
    Returns:
        List of GRN strings within the interval [start_grn, end_grn]
    """
    start_float = parse_grn_str2float(start_grn)
    end_float = parse_grn_str2float(end_grn)
    
    if grns_str is not None:
        # Filter from provided list
        return [g for g in grns_str if start_float <= parse_grn_str2float(g) <= end_float]
    else:
        # Auto-generate with 0.01 stepsize
        result = []
        
        # Determine if we're in a TM region or loop region
        if '.' in start_grn and len(start_grn.split('.')[0]) == 1:
            # TM region (e.g., 1.28 to 1.64)
            helix = int(start_grn.split('.')[0])
            start_pos = int(round((start_float - helix) * 100))
            end_pos = int(round((end_float - helix) * 100))
            
            for pos in range(start_pos, end_pos + 1):
                grn = f"{helix}.{pos:02d}"
                result.append(grn)
        
        return result


def init_grn_intervals(grn_config_str):

    grns_str_str = []
    if grn_config_str:
        for region_name, (start_grn, end_grn) in grn_config_str.items():
            # Generate GRNs for this interval (auto-generate since we don't have a predefined list)
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_str.extend(region_grns)
    return grns_str_str





class GRNConfigManager:
    """Manages GRN configurations for different protein families."""
    
    def __init__(self, paths=None):
        from protos.io.paths import ProtosPaths
        
        # Use ProtosPaths for all path management
        self.paths = paths or ProtosPaths()
        
        # Get config directory through ProtosPaths
        self.config_dir = Path(self.paths.get_processor_path("grn")) / "configs"
        
        self.configs = {}
        self._load_configs()
    
    def _load_configs(self):
        """Load all available configuration files."""
        if not self.config_dir.exists():
            logger.warning(f"Config directory not found: {self.config_dir}")
            return
            
        # Check for main config.json with multiple families
        main_config = self.config_dir / 'config.json'
        print("main config", main_config)
        if main_config.exists():
            try:
                with open(main_config, 'r') as f:
                    all_configs = json.load(f)
                    # If it contains protein families as keys, load them
                    if isinstance(all_configs, dict):
                        for family_name, family_config in all_configs.items():
                            if isinstance(family_config, dict):
                                self.configs[family_name] = family_config

            except Exception as e:
                logger.error(f"Error loading main config {main_config}: {e}")
        print(self.configs)

    def get_intervals(self, protein_family):
        """Get GRN intervals for a specific protein family."""
        if protein_family not in self.configs:
            logger.warning(f"No config found for protein family: {protein_family}")
            return {}
            
        config = self.configs[protein_family]
        intervals = {}
        
        # Handle the structure from config.json (with standard/strict subdivisions)
        if 'standard' in config or 'strict' in config:
            # Use standard intervals by default
            interval_data = config.get('standard', config.get('strict', {}))
            for region_name, (start, end) in interval_data.items():
                intervals[f"{start}-{end}"] = (
                    parse_grn_str2float(start),
                    parse_grn_str2float(end)
                )
        # Handle the structure with grn_ranges
        elif 'grn_ranges' in config:
            for grn_range in config['grn_ranges']:
                start_grn = grn_range['start']
                end_grn = grn_range['end']
                intervals[f"{start_grn}-{end_grn}"] = (
                    parse_grn_str2float(start_grn),
                    parse_grn_str2float(end_grn)
                )
        
        return intervals
    
    def get_config(self, protein_family, strict=False):
        """Get configuration for a specific protein family.
        
        Args:
            protein_family: Name of the protein family
            strict: Whether to return strict or standard config
            
        Returns:
            Configuration dictionary
        """
        if protein_family not in self.configs:
            logger.warning(f"No config found for protein family: {protein_family}")
            return {}
            
        config = self.configs[protein_family]
        
        # If config has strict/standard subdivisions
        if strict and 'strict' in config:
            return config['strict']
        elif not strict and 'standard' in config:
            return config['standard']
        else:
            return config


# ---------------------------------------------------------------------------
#   6. GRN Interval Calculation
# ---------------------------------------------------------------------------
def get_interval_for_grn(grn_target, intervals):
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


def generate_all_grns_from_config(grn_config, grns_str=None):
    """Generate all GRN positions from a configuration dictionary.
    
    Args:
        grn_config: Dictionary with TM regions as keys and [start, end] as values
        grns_str: Optional list of valid GRN strings to filter from.
                 If None, auto-generates GRNs.
        
    Returns:
        List of all GRN strings within the defined intervals
    """
    all_grns = []
    
    for region_name, (start_grn, end_grn) in grn_config.items():
        # Use get_grn_interval to get GRNs for this region
        region_grns = get_grn_interval(start_grn, end_grn, grns_str)
        all_grns.extend(region_grns)
    
    # Remove duplicates and sort
    all_grns = list(set(all_grns))
    return sort_grns_str(all_grns)


# ---------------------------------------------------------------------------
#   7. Visualization Support
# ---------------------------------------------------------------------------
def map_grn_to_color(grn, color_map=None):
    """Map a GRN to a color for visualization."""
    if color_map is None:
        # Default color scheme using RGB values
        color_map = {
            'n_term': 'rgb(31, 119, 180)',     # N-terminal
            'c_term': 'rgb(255, 127, 14)',     # C-terminal  
            'tm': 'rgb(214, 39, 40)',          # Transmembrane regions
            'loop': 'rgb(44, 160, 44)',        # Loop regions
        }
    
    grn_float = parse_grn_str2float(grn)
    
    # N-terminal
    if grn_float < 0:
        return color_map.get('n_term', 'rgb(31, 119, 180)')
    # C-terminal
    elif grn_float >= 100:
        return color_map.get('c_term', 'rgb(255, 127, 14)')
    # Loops (between 10 and 100)
    elif grn_float >= 10:
        return color_map.get('loop', 'rgb(44, 160, 44)')
    # TM regions (less than 10)
    else:
        return color_map.get('tm', 'rgb(214, 39, 40)')


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

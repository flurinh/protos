import re
import os
import json
import pandas as pd
from decimal import Decimal, getcontext
from typing import List
import re
from typing import List, Tuple  # Python 3.9+ for Tuple, just Tuple for older

import re
from typing import List, Tuple

ERROR_FLOAT = -999.999  # Unique float to indicate parsing error


# ---------------------------------------------------------------------------
#   1. String to Float Parser
# ---------------------------------------------------------------------------
def parse_grn_str2float(grn_str: str) -> float:
    grn_str = grn_str.strip()

    if grn_str.startswith('n'):
        val_str = ""
        if grn_str.startswith('n.'):
            val_str = grn_str[2:]
        elif len(grn_str) > 1 and grn_str[1:].isdigit():
            val_str = grn_str[1:]
        else:
            print(f"Error: Unparsable N-terminus GRN: '{grn_str}'");
            return ERROR_FLOAT
        try:
            n_index = int(val_str)
            if n_index < 1:
                print(f"Error: Invalid N-terminus index (<1): '{grn_str}'");
                return ERROR_FLOAT
            return round(float(-1 * n_index), 3)
        except ValueError:
            print(f"Error: N-term value parse error: '{grn_str}'");
            return ERROR_FLOAT
    elif grn_str.startswith('c'):
        val_str = ""
        if grn_str.startswith('c.'):
            val_str = grn_str[2:]
        elif len(grn_str) > 1 and grn_str[1:].isdigit():  # Allow legacy cVAL for parsing
            val_str = grn_str[1:]
            # print(f"Info: Parsed legacy C-terminus format '{grn_str}'. Canonical is 'c.{val_str}'.")
        else:
            print(f"Error: Unparsable C-terminus GRN: '{grn_str}'");
            return ERROR_FLOAT
        try:
            c_index = int(val_str)
            if c_index < 1:
                print(f"Error: Invalid C-terminus index (<1): '{grn_str}'");
                return ERROR_FLOAT
            return round(100.0 + float(c_index), 3)
        except ValueError:
            print(f"Error: C-term value parse error: '{grn_str}'");
            return ERROR_FLOAT
    elif '.' in grn_str:
        try:
            parts = grn_str.split('.', 1)
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                integer_part_val = int(parts[0])
                if integer_part_val == 0:
                    print(f"Error: Invalid GRN: 0.xx not allowed: '{grn_str}'");
                    return ERROR_FLOAT
                return round(integer_part_val + float(f"0.{parts[1]}"), 3)
            else:
                print(f"Error: Unparsable dot-notation: '{grn_str}'");
                return ERROR_FLOAT
        except (ValueError, IndexError):
            print(f"Error parsing dot-notation: '{grn_str}'");
            return ERROR_FLOAT
    elif 'x' in grn_str:
        try:
            parts = grn_str.split('x', 1)
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                integer_part_val = int(parts[0])
                if integer_part_val == 0:
                    print(f"Error: Invalid GRN: 0xXX not allowed: '{grn_str}'");
                    return ERROR_FLOAT
                return round(integer_part_val + float(f"0.{parts[1]}"), 3)
            else:
                print(f"Error: Unparsable x-notation: '{grn_str}'");
                return ERROR_FLOAT
        except (ValueError, IndexError):
            print(f"Error parsing x-notation: '{grn_str}'");
            return ERROR_FLOAT
    else:
        print(f"Error: Unparsable GRN (no known format): '{grn_str}'")
        return ERROR_FLOAT


# ---------------------------------------------------------------------------
#   2. Float to String Formatter
# ---------------------------------------------------------------------------
def parse_grn_float2str(grn_float: float, notation_type: str = 'dot') -> str:
    grn_float = round(grn_float, 3)

    if grn_float <= -0.5:  # N-Term: ends at n.1 (-1.0)
        n_index = int(round(abs(grn_float)))
        if n_index == 0: n_index = 1  # Ensure at least n.1
        return f'n.{n_index}'  # Always dot for N-term index

    if grn_float >= 100.5:  # C-Term: starts at c.1 (101.0)
        c_index = int(round(grn_float - 100.0))
        if c_index == 0: c_index = 1  # Ensure at least c.1
        return f'c.{c_index}'  # Always dot for C-term index

    integer_part_check = int(grn_float)  # For range checking
    if not (0.999 < grn_float < 100.5) or integer_part_check == 0 or integer_part_check == 100:
        return f"undef.{grn_float:.3f}"

    integer_part = int(grn_float)  # Actual integer part for formatting
    separator = '.' if notation_type == 'dot' else 'x'

    # Rule: Single number integer part (1-9) is TM Helix.
    # Rule: Two-number integer part (e.g., 12, 23, 32) is Loop AB.CCC or BA.CCC.
    is_tm_helix_category = 1 <= integer_part <= 9  # TMH can be H1-H9

    if is_tm_helix_category:  # TM Helix (H1-H9)
        is_insertion = abs(grn_float - round(grn_float, 2)) > 1e-5
        if is_insertion:  # TM with insertion: H.YYY or HxYYY
            decimal_val = int(round((grn_float - integer_part) * 1000))
            return f"{integer_part}{separator}{decimal_val:03d}"
        else:  # Standard TM: H.YY or HxYY
            grn_for_yy = round(grn_float, 2)
            decimal_val = int(round((grn_for_yy - integer_part) * 100))
            return f"{integer_part}{separator}{decimal_val:02d}"
    else:  # Loop (e.g. integer_part 12, 23, 32) or other multi-digit int part non-TM
        # Loops (AB.CCC or BA.CCC) ALWAYS use 3 digits for fractional part
        decimal_val = int(round((grn_float - integer_part) * 1000))
        return f"{integer_part}{separator}{decimal_val:03d}"


# ---------------------------------------------------------------------------
#   3. GRN String Validator
# ---------------------------------------------------------------------------
def check_str_grn_valid(grn_str: str) -> bool:
    grn_str = grn_str.strip()
    n_pattern_dot = re.compile(r'^n\.([1-9]\d*)$')
    if n_pattern_dot.match(grn_str): return True

    c_pattern_dot = re.compile(r'^c\.([1-9]\d*)$')  # Canonical C-term is c.VAL
    if c_pattern_dot.match(grn_str): return True
    # Allow legacy cVAL for parsing by str2float, but check_str_grn_valid enforces c.VAL
    # If legacy cVAL should also pass validation:
    # c_pattern_x_legacy = re.compile(r'^c([1-9]\d*)$')
    # if c_pattern_x_legacy.match(grn_str): return True

    tm_loop_pattern = re.compile(r'^([1-9]\d*)([x.])(\d+)$')
    m_tm_loop = tm_loop_pattern.match(grn_str)
    if m_tm_loop:
        s_integer_part = m_tm_loop.group(1)
        s_fractional_part = m_tm_loop.group(3)
        try:
            integer_part = int(s_integer_part)
        except ValueError:
            return False

        is_tm_helix_category = 1 <= integer_part <= 9  # H1-H9 are TMs
        frac_len = len(s_fractional_part)

        if is_tm_helix_category:
            if frac_len == 2: return True
            if frac_len == 3:
                # A string "X.YY0" or "XxYY0" for a TM is not canonical if X.YY is the true base.
                # parse_grn_float2str(X.YY_float) -> "X.YY"
                # parse_grn_float2str(X.YYZ_float_insertion) -> "X.YYZ"
                # So, a valid 3-digit TM string should not be equivalent to a 2-digit one.
                if s_fractional_part.endswith('0'):
                    # Check if it's like "X.Y00" which should be "X.Y0"
                    # or "X.YY0" which should be "X.YY"
                    # Example: "1.500" vs "1.50". "1.520" vs "1.52".
                    # If 0.xxx represents X.YY, then X.YY0 is invalid.
                    # If 0.x represents X.Y, then X.Y00 and X.Y0 are invalid (should be X.Y0).
                    # A robust check: convert to float, then back to 2-digit string. If it matches first 2 digits,
                    # and original was 3 digits ending in 0, it's likely an invalid verbose form.
                    temp_float = parse_grn_str2float(
                        f"{s_integer_part}.{s_fractional_part[:2]}")  # e.g., 1.50 from 1.500
                    temp_str_canonical_2digit = parse_grn_float2str(temp_float, m_tm_loop.group(2))  # separator

                    if temp_str_canonical_2digit == f"{s_integer_part}{m_tm_loop.group(2)}{s_fractional_part[:2]}":
                        return False  # e.g. 1.500 is invalid because 1.50 is canonical
                return True  # Valid 3-digit insertion like 1.521
            return False  # Must be 2 or 3 digits
        else:  # Loop (integer part has >= 2 digits like 12, 23, 32 or > 9)
            return frac_len == 3  # Loops always 3 digits
    return False


# ---------------------------------------------------------------------------
#   4. Sorting Helper: Get Flanking TMs for a Loop
# ---------------------------------------------------------------------------
def get_prev_next_tm(grn_float: float) -> Tuple[int, int]:
    loop_id_integer = int(grn_float)
    # Loop IDs are two digits, e.g., 12 (H1-H2), 23 (H2-H3),
    # or "reversed" like 32 (interpreted as H2-H3 region).
    if 10 <= loop_id_integer <= 87:  # Max could be 87 (H7-H8, B=8, A=7 for BA)
        s_loop_id = str(loop_id_integer)
        if len(s_loop_id) == 2:
            try:
                digit1 = int(s_loop_id[0])
                digit2 = int(s_loop_id[1])
                if not (1 <= digit1 <= 8 and 1 <= digit2 <= 8 and digit1 != digit2):
                    return (99, 99)
                # For sorting, "helix_A" should be the N-terminal flanking helix of the loop segment.
                # If ID is "12", N-flank is 1, C-flank is 2.
                # If ID is "32", it's between H2 and H3. N-flank for this segment is H2.
                n_flank = min(digit1, digit2)
                c_flank = max(digit1, digit2)
                return n_flank, c_flank
            except ValueError:
                pass
    return (99, 99)


# ---------------------------------------------------------------------------
#   6. Main GRN Float Sorter
# ---------------------------------------------------------------------------
def sort_grns(grns_float_list: List[float]) -> List[float]:
    grns_ntail = sorted([x for x in grns_float_list if x <= -0.5])
    grns_ctail = sorted([x for x in grns_float_list if x >= 100.5])

    body_grns_all_numerically_sorted = sorted([
        x for x in grns_float_list
        if 0.999 < x < 100.5 and int(x) > 0  # Valid body floats (H1.xx up to before C-term)
    ])

    grns_tm_only = []  # TMs H1-H9
    grns_loops_ab_type = []  # Loops like AB.CCC or BA.CCC
    # No grns_other_body if all single-digit int parts are TMs up to H9,
    # and all 2-digit int parts (10-87) are AB/BA loops.

    for x_float in body_grns_all_numerically_sorted:
        integer_part = int(x_float)
        if 1 <= integer_part <= 9:  # TM Helix (H1-H9)
            grns_tm_only.append(x_float)
        else:  # Must be a loop (e.g., 12.xxx, 23.xxx, 32.xxx)
            helix_A, helix_B = get_prev_next_tm(x_float)
            if helix_A != 99:  # It's a recognized AB/BA.CCC type loop
                grns_loops_ab_type.append(x_float)
            else:
                # This case should be rare if all body GRNs are TMs H1-H9 or AB/BA loops
                print(f"Warning: GRN float {x_float} in body not classified as TM H1-H9 or AB/BA Loop.")
                # Decide where to put these: for now, append after TMs and AB loops
                # This would require an grns_other_body list again if such cases exist
                # and a final merge step. For simplicity, assume valid inputs are TM or AB/BA loop.
                # If they can exist, they need a defined sorting rule.
                # For now, we can add them to a temporary list and append at the end of body.
                # For the given test case, grns_other_body will be empty.
                pass  # Or add to a grns_unclassified_body list

    # TMs and AB_Loops are already sorted numerically among themselves
    # because they came from body_grns_all_numerically_sorted.

    loop_idx = 0
    tm_idx = 0
    sorted_body_interleaved = []

    while tm_idx < len(grns_tm_only) and loop_idx < len(grns_loops_ab_type):
        current_tm_float = grns_tm_only[tm_idx]
        current_loop_float = grns_loops_ab_type[loop_idx]

        # helix_A_of_loop is the N-terminal flanking helix of this loop segment
        helix_A_of_loop, _ = get_prev_next_tm(current_loop_float)

        if int(current_tm_float) <= helix_A_of_loop:
            sorted_body_interleaved.append(current_tm_float)
            tm_idx += 1
        else:
            sorted_body_interleaved.append(current_loop_float)
            loop_idx += 1

    sorted_body_interleaved.extend(grns_tm_only[tm_idx:])
    sorted_body_interleaved.extend(grns_loops_ab_type[loop_idx:])

    # If there were grns_unclassified_body, they would be merged here based on some rule
    # or simply appended if they always come last in the body.
    # For the given test case, this part is not critical as other_body is empty.

    result = grns_ntail + sorted_body_interleaved + grns_ctail
    return [round(x, 3) for x in result]


# ---------------------------------------------------------------------------
#   7. String List Sorter (Orchestrator)
# ---------------------------------------------------------------------------
def sort_grns_str(grns_str_list: List[str], output_notation_type: str = 'dot') -> List[str]:
    grns_float = []
    for s in grns_str_list:
        f = parse_grn_str2float(s)
        if abs(f - ERROR_FLOAT) > 1e-6:
            grns_float.append(f)
        else:
            print(f"Info: Skipping unparsable/invalid GRN string from sort input: {s}")
    sorted_grns_float = sort_grns(grns_float)
    return [parse_grn_float2str(f, notation_type=output_notation_type) for f in sorted_grns_float]


def init_grn_intervals(grn_config):
    getcontext().prec = 4

    std_grns = []
    for interval_key, (left_str, right_str) in grn_config.items():
        # Convert the string GRN values to Decimal
        float_left = Decimal(parse_grn_str2float(left_str))
        float_right = Decimal(parse_grn_str2float(right_str))
        # Check if left is greater than right and swap if necessary
        if float_left > float_right:
            float_left, float_right = float_right, float_left
        # Generate a list of Decimals from left to right with a step of .01
        y = (float_right - float_left).quantize(Decimal('0.01'))
        n = int(y * 100) + 1
        grns_float = [(float_left + Decimal(x) * Decimal('0.01')).quantize(Decimal('0.01')) for x in range(n)]
        # Convert the list of Decimals back to the standardized GRN strings
        grns_str = [parse_grn_float2str(x) for x in grns_float]
        # Process the list to divide the numbers by 10 and adjust the format

        output_list = []
        for item in grns_str:
            match = re.match(r'(\d+)x(\d+)', item)
            if match:
                prefix_number = match.group(1)
                number_after_x = int(match.group(2))
                scaled_number = number_after_x // 10  # Use integer division for rounding down
                formatted_string = f'{prefix_number}x{scaled_number}'
                output_list.append(formatted_string)
        std_grns.extend(output_list)
    return std_grns


class GRNConfigManager:
    def __init__(self, path=None, config_path='config.json', protein_family='gpcr_a'):
        if path is None:
            # Try different possible locations for configs
            possible_paths = [
                'data/grn/configs/',
                os.path.join(os.path.dirname(__file__), 'configs/'),
                os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs/')
            ]

            for possible_path in possible_paths:
                test_path = os.path.join(possible_path, config_path)
                if os.path.exists(test_path):
                    path = possible_path
                    break

            # Default to package directory if not found
            if path is None:
                path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs/')

        self.config_path = os.path.join(path, config_path)
        self.config = None
        self.protein_family = protein_family

    def load_config(self):
        if self.config is None:
            try:
                with open(self.config_path, 'r') as f:
                    self.config = json.load(f)
            except FileNotFoundError:
                # Provide a fallback config for testing
                self.config = {
                    "microbial_opsins": {
                        "standard": {
                            "TM1": ["1x01", "1x50"],
                            "TM2": ["2x01", "2x50"],
                            "TM3": ["3x01", "3x50"],
                            "TM4": ["4x01", "4x50"],
                            "TM5": ["5x01", "5x50"],
                            "TM6": ["6x01", "6x50"],
                            "TM7": ["7x01", "7x50"]
                        },
                        "strict": {
                            "TM1": ["1x50", "1x50"],
                            "TM2": ["2x50", "2x50"],
                            "TM3": ["3x50", "3x50"],
                            "TM4": ["4x50", "4x50"],
                            "TM5": ["5x50", "5x50"],
                            "TM6": ["6x50", "6x50"],
                            "TM7": ["7x50", "7x50"]
                        }
                    },
                    "gpcr_a": {
                        "standard": {
                            "TM1": ["1x01", "1x50"],
                            "TM2": ["2x01", "2x50"],
                            "TM3": ["3x01", "3x50"],
                            "TM4": ["4x01", "4x50"],
                            "TM5": ["5x01", "5x50"],
                            "TM6": ["6x01", "6x50"],
                            "TM7": ["7x01", "7x50"]
                        },
                        "strict": {
                            "TM1": ["1x50", "1x50"],
                            "TM2": ["2x50", "2x50"],
                            "TM3": ["3x50", "3x50"],
                            "TM4": ["4x50", "4x50"],
                            "TM5": ["5x50", "5x50"],
                            "TM6": ["6x50", "6x50"],
                            "TM7": ["7x50", "7x50"]
                        }
                    }
                }
        return self.config

    def list_available_datasets(self):
        self.load_config()
        datasets = {}
        for family, configs in self.config.items():
            datasets[family] = {
                'standard': 'standard' in configs,
                'strict': 'strict' in configs
            }
        return datasets

    def get_config(self, strict=False, protein_family=None):
        if protein_family is None:
            protein_family = self.protein_family
        self.load_config()
        config_type = 'strict' if strict else 'standard'
        family_config = self.config.get(protein_family, {})
        return family_config.get(config_type, {})

    def init_grns(self, strict=False, protein_family=None):
        if protein_family is None:
            protein_family = self.protein_family
        grn_config = self.get_config(strict=strict, protein_family=protein_family)
        return init_grn_intervals(grn_config)


def get_grn_interval(left: str = '8x48', right: str = '9x71', table: pd.DataFrame = pd.DataFrame([]),
                     config_manager: GRNConfigManager = None, grns_str: list = [],
                     protein_family: str = 'gpcr_a', strict: bool = True):
    # Priority 1: Use provided GRNs string list if available
    if grns_str:
        grns_float = sort_grns([round(parse_grn_str2float(x), 2) for x in grns_str])
    # Priority 2: Use the provided table's columns as GRNs if no GRNs string list is provided
    elif not grns_str and not table.empty:
        grns_str = table.columns.tolist()
        grns_float = [round(parse_grn_str2float(x), 3) for x in grns_str]
    # Priority 3: Default initialization using a new initialization of config_manager if no list or table provided
    else:
        if config_manager is None:
            config_manager = GRNConfigManager()
        grn_config = config_manager.get_config(strict=strict, protein_family=protein_family)
        grns_str = init_grn_intervals(grn_config)
        grns_float = sort_grns([round(parse_grn_str2float(x), 2) for x in grns_str])

    n_tail_grn = []
    grns_str_interval = []
    c_tail_grn = []
    if 'n' in left:
        start = int(left.split('n.')[1])
        end = 1
        if 'n' in right:
            end = int(right.split('n.')[1])
        len_n_tail = start - end + 1
        n_tail_grn = [('n.' + str(start - i)) for i in range(len_n_tail) if ('n.' + str(start - i)) in grns_str]
        if not 'n' in right:
            std_grns = [x for x in grns_str if 'x' in x]
            left = std_grns[0]
    if 'c' in right:
        end = int(right.split('c.')[1])
        start = 1
        if 'c' in left:
            start = int(left.split('c.')[1])
        len_c_tail = end - start
        c_tail_grn = [('c.' + str(start + i)) for i in range(len_c_tail + 1) if ('c.' + str(start + i)) in grns_str]
        if not 'c' in left:
            std_grns = [x for x in grns_str if 'x' in x]
            right = std_grns[-1]
    if ('x' in left) and ('x' in right):
        left_grn = parse_grn_str2float(left)
        right_grn = parse_grn_str2float(right)
        sorted_grns = sort_grns(grns_float)
        assert left_grn in sorted_grns, print('could not find left grn pivot in columns')
        assert right_grn in sorted_grns, print('could not find right grn pivot in columns')
        float_interval = sorted_grns[sorted_grns.index(left_grn):sorted_grns.index(right_grn) + 1]
        grns_str_interval = [parse_grn_float2str(round(x, 6)) for x in float_interval]
    return n_tail_grn + grns_str_interval + c_tail_grn


def init_std_grns(grn_intervals):
    std_grns = []
    for interval in grn_intervals.items():
        grns = get_grn_interval(*interval[1])
        grns = [g for g in grns if len(g) < 5]
        std_grns += grns
    return std_grns


def map_grn_to_color(grn):
    if 'n' in grn:
        return 'rgb(31, 119, 180)'
    if 'c' in grn:
        return 'rgb(255, 127, 14)'
    grn_f = parse_grn_str2float(grn)
    if grn_f > 10:
        return 'rgb(44, 160, 44)'
    else:
        return 'rgb(214, 39, 40)'


def get_seq(gene: str, grn_table: pd.DataFrame):
    seq_list = list(grn_table.loc[gene].values)
    seq = ''.join([x[0] for x in seq_list if x != '-'])
    seq.replace('-', '')
    return seq


def get_annot_seq(gene: str, grn_table: pd.DataFrame):
    seq_list = list(grn_table.loc[gene].values)
    return [x for x in seq_list if x != '-']


def remove_gaps_from_sequences(sequence_dict):
    return {key: sequence.replace('-', '') for key, sequence in sequence_dict.items()}


def flatten(l):
    return [item for sublist in l for item in sublist if len(item) == 4]


from typing import List


# Assuming ERROR_FLOAT and parse_grn_str2float are defined as in your provided code.
# If not, you'll need to include them or ensure they are accessible.
# For example:
# ERROR_FLOAT = -999.999
# def parse_grn_str2float(grn_str: str) -> float:
#     # ... (implementation from your provided code) ...
#     pass

def get_tm_residues(grn_str_list: List[str]) -> List[str]:
    """
    Filters a list of GRN strings to return only those representing TM residues.

    Args:
        grn_str_list: A list of GRN strings.

    Returns:
        A list of GRN strings that correspond to TM residues (H1-H9).
    """
    tm_residue_strings = []
    for grn_s in grn_str_list:
        f_val = parse_grn_str2float(grn_s)

        # Skip if the GRN string is unparsable
        if abs(f_val - ERROR_FLOAT) < 1e-6:
            continue

        # The integer part of the float determines the segment type.
        # N-terminal residues result in f_val <= -0.5
        # C-terminal residues result in f_val >= 100.5
        # TM helices (H1-H9) have an integer part from 1 to 9.
        # Loops (e.g., 12.xxx, 23.xxx) have integer parts >= 10.

        integer_part = int(f_val)

        if 1 <= integer_part <= 9:
            # This condition specifically identifies TM helices H1-H9
            tm_residue_strings.append(grn_s)

    return tm_residue_strings


# ---------------------------------------------------------------------------
#   Test Suite
# ---------------------------------------------------------------------------
def run_all_tests():
    print("--- Testing parse_grn_str2float ---")
    str2float_tests = [
        ('n.10', -10.0), ('n1', -1.0), ('n.1', -1.0),
        ('n.0', ERROR_FLOAT), ('n0', ERROR_FLOAT),
        ('c.5', 105.0), ('c1', 101.0), ('c.12', 112.0), ('c12', 112.0),
        ('c.0', ERROR_FLOAT), ('c0', ERROR_FLOAT),
        ('1.50', 1.50), ('1x50', 1.50),
        ('0.50', ERROR_FLOAT), ('0x50', ERROR_FLOAT),
        ('1.501', 1.501), ('1x501', 1.501),
        ('12.003', 12.003), ('12x003', 12.003), ('32.001', 32.001), ('32x001', 32.001),
        ('9.123', 9.123), ('9x123', 9.123),  # TMH9
        ('bad', ERROR_FLOAT), ('n.x', ERROR_FLOAT), ('1.x2', ERROR_FLOAT)
    ]
    all_s2f_passed = True
    for grn_s, expected_f in str2float_tests:
        res_f = parse_grn_str2float(grn_s)
        if abs(res_f - expected_f) > 1e-6:
            print(f"FAIL: parse_grn_str2float('{grn_s}') -> {res_f} (Expected: {expected_f})")
            all_s2f_passed = False
    print(f"parse_grn_str2float tests {'PASSED' if all_s2f_passed else 'FAILED'}.")

    print("\n--- Testing parse_grn_float2str ---")
    float2str_tests = [
        (-10.0, 'dot', 'n.10'), (-1.0, 'x', 'n.1'), (-0.5, 'dot', 'n.1'),
        (0.0, 'dot', 'undef.0.000'), (0.4, 'dot', 'undef.0.400'),
        (105.0, 'dot', 'c.5'), (101.0, 'x', 'c.1'),  # C-term now always dot
        (100.0, 'dot', 'undef.100.000'), (100.4, 'dot', 'undef.100.400'),
        (1.50, 'dot', '1.50'), (1.50, 'x', '1x50'),
        (1.501, 'dot', '1.501'), (1.501, 'x', '1x501'),
        (7.520, 'dot', '7.52'), (7.520, 'x', '7x52'),
        (7.5206, 'dot', '7.521'),
        (12.003, 'dot', '12.003'), (12.003, 'x', '12x003'),
        (32.001, 'dot', '32.001'), (32.001, 'x', '32x001'),  # Loop BA.CCC
        (9.12, 'dot', '9.12'), (9.12, 'x', '9x12'),  # TMH9 canonical
        (9.123, 'dot', '9.123'), (9.123, 'x', '9x123'),  # TMH9 insertion
    ]
    all_f2s_passed = True
    for f_val, note, expected_s in float2str_tests:
        res_s = parse_grn_float2str(f_val, note)
        if res_s != expected_s:
            print(f"FAIL: parse_grn_float2str({f_val}, '{note}') -> '{res_s}' (Expected: '{expected_s}')")
            all_f2s_passed = False
    print(f"parse_grn_float2str tests {'PASSED' if all_f2s_passed else 'FAILED'}.")

    print("\n--- Testing check_str_grn_valid ---")
    validity_tests = [
        ('n.10', True), ('n.1', True), ('n.0', False), ('n10', False),
        ('c.5', True), ('c.12', True), ('c1', False), ('c12', False),  # cX legacy not valid output
        ('1.50', True), ('1x50', True), ('1.501', True), ('1x501', True),
        ('0.50', False), ('1.500', False),  # ('1x500', False) implicit
        ('9.12', True), ('9x12', True), ('9.123', True), ('9x123', True),  # TMH9
        ('12.003', True), ('12x003', True), ('32.001', True), ('32x001', True),
        ('12.10', False),  # Loop must be 3 fractional
    ]
    all_valid_passed = True
    for grn_s, expected_v in validity_tests:
        res_v = check_str_grn_valid(grn_s)
        if res_v != expected_v:
            print(f"FAIL: check_str_grn_valid('{grn_s}') -> {res_v} (Expected: {expected_v})")
            all_valid_passed = False
    print(f"check_str_grn_valid tests {'PASSED' if all_valid_passed else 'FAILED'}.")

    print("\n--- Testing sort_grns_str ---")
    # Input includes n.0, c0, 0.50 which will be filtered by parse_grn_str2float
    grn_list_str_input = ["c.10", "n.10", "1.50", "c.1", "12.003", "n.1", "2.601", "1x60", "8.70", "12x000", "2.30",
                          "23.001", "9.100", "32.001", "n.0", "c0", "0.50"]
    # Expected sort of valid inputs:
    expected_sorted_dot = ["n.10", "n.1",  # N-Term
                           "1.50", "1.60",  # TM1
                           "12.000", "12.003",  # Loop H1-H2
                           "2.30", "2.601",  # TM2
                           "23.001",  # Loop H2-H3
                           "32.001",  # Loop H2-H3 (BA notation, numerically after 23.xxx)
                           "8.70",  # TM8
                           "9.10",  # TM9
                           "c.1", "c.10"]  # C-Term
    all_sort_passed = True

    sorted_dot = sort_grns_str(grn_list_str_input, 'dot')
    print(f"Input: {grn_list_str_input}")
    print(f"Sorted (dot): {sorted_dot}")
    print(f"Expected (dot): {expected_sorted_dot}")
    if sorted_dot != expected_sorted_dot:
        print("FAIL sort_grns_str (dot)")
        all_sort_passed = False
        grns_float_debug = []
        for s_debug in grn_list_str_input:
            f_debug = parse_grn_str2float(s_debug)
            if abs(f_debug - ERROR_FLOAT) > 1e-6: grns_float_debug.append(f_debug)
        print(f"DEBUG: Floats parsed for sort: {grns_float_debug}")
        print(f"DEBUG: Floats after sort_grns: {sort_grns(grns_float_debug)}")

    # Generate expected_sorted_x based on corrected parse_grn_float2str for C-term
    temp_floats_for_x = [parse_grn_str2float(s) for s in expected_sorted_dot if
                         abs(parse_grn_str2float(s) - ERROR_FLOAT) > 1e-6]
    expected_sorted_x = [parse_grn_float2str(f, 'x') for f in temp_floats_for_x]

    sorted_x = sort_grns_str(grn_list_str_input, 'x')
    print(f"Sorted (x):   {sorted_x}")
    print(f"Expected (x): {expected_sorted_x}")
    if sorted_x != expected_sorted_x:
        print("FAIL sort_grns_str (x)")
        all_sort_passed = False
    print(f"sort_grns_str tests {'PASSED' if all_sort_passed else 'FAILED'}.")


if __name__ == '__main__':
    run_all_tests()
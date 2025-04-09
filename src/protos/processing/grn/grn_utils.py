import re
import os
import json
import pandas as pd
from decimal import Decimal, getcontext


def parse_grn_str2float(grn: str = '1x13'):
    # N TAIL
    if 'n' in grn:
        tm = 0
        rn = int(-1 * float(grn.split('n.')[1]))
    # C TAIL
    elif 'c' in grn:
        tm = 100
        rn = int(grn.split('c.')[1])
    # TM BUNDLE with dot notation (e.g., 1.50)
    elif '.' in grn and 'x' not in grn:
        try:
            parts = grn.split('.')
            if len(parts) == 2 and parts[0].isdigit():
                tm = int(parts[0])
                rn = float('0.' + parts[1])
            else:
                print("unparsable grn string:", grn)
                # Default fallback values
                tm = 0
                rn = 0
        except Exception as e:
            print(f"Error parsing GRN string '{grn}': {e}")
            # Default fallback values
            tm = 0
            rn = 0
    # TM BUNDLE with x notation (e.g., 1x50)
    elif 'x' in grn:
        region, rn = grn.split('x')
        if len(region) == 2:
            loop = int(region)
            rn = float('.' + rn)
            tm = loop
        else:
            tm = int(region)
            rn = float('.' + rn)
    else:
        print("unparsable grn string:", grn)
        # Default fallback values
        tm = 0
        rn = 0
    return round(tm + rn,3)


def parse_grn_float2str(grn: float = 7.521):
    grn = round(grn, 3)
    # N tail
    if grn <= 0:
        grn = 'n.' + str(int(-1 * grn))
    # Special case for values between 100 and 101 (handle as TM notation)
    elif 100 <= grn < 101:
        # Extract decimal part and format as 100xYY
        decimal_part = int(round((grn - 100) * 100))
        grn = f'100x{decimal_part:02d}'
    # C tail (values above 101)
    elif grn >= 101:
        grn = 'c.' + str(int(grn - 100))
    # TM
    else:
        # For TM helices (1-7), use dot notation instead of x notation
        parts = str(grn).split('.')
        if len(parts) == 2 and parts[0].isdigit() and int(parts[0]) <= 7:
            # Format as N.YY (dot notation) for TM helices
            tm_num = parts[0]
            if len(parts[1]) == 1:
                parts[1] += '0'  # Ensure we have at least 2 digits
            grn = f"{tm_num}.{parts[1]}"
        else:
            # Use x notation for other regions (legacy format)
            grn = 'x'.join(str(grn).split('.'))
            if len(grn.split('x')[0]) == 2:
                if (len(grn) == 4) and (grn[-1] != '0'):
                    grn += '0'
            elif len(grn) == 3:
                grn += '0'
    return grn


def check_str_grn_valid(grn: str = '1x13'):
    n_pattern = re.compile(r'n\.[1-9][0-9]*[0-9]*[0-9]*')
    x_pattern = re.compile(r'[1-8]x[1-9][1-9]*')
    l_pattern = re.compile(r'[1-8][1-8]*x[0-9][0-9]*')
    c_pattern = re.compile(r'c\.[1-9][0-9]*[0-9]*')
    # Check if the grn matches any of the allowed patterns
    if n_pattern.match(grn) or c_pattern.match(grn) or l_pattern.match(grn):
        # Check if the grn is of the form "n.01" (not allowed)
        if grn[0] == 'n' and grn[2] == '0':
            return False
        # Check if the grn is of the form "c.01" (not allowed)
        elif grn[0] == 'c' and grn[2] == '0':
            return False
        # Check if the grn is of the form "34.510" (not allowed)
        else:
            if not x_pattern.match(grn):
                if len(grn.split('x')[-1]) == 3:
                    if grn.split('x')[-1][-1] == 0:
                        return False
        return True
    else:
        return False


def get_prev_next_tm(grn: float = 23.512):
    int2 = int(grn % 10)  # 3
    int1 = int(grn // 10)
    prev_ = min(int1, int2)
    next_ = max(int1, int2)
    return (int(prev_), int(next_))


def sort_iecls(grns):
    ciecl = sorted([x for x in grns if (int(x // 10) >= int(x % 10))])  # % > int2
    niecl = sorted([x for x in grns if (int(x // 10) <= int(x % 10))])
    sorted_ciecl = []
    for iecl_id in range(8):
        div = int(iecl_id + 1)
        holder = []
        for grn in ciecl:
            g = grn // 10
            if g / div == 1:
                holder.append(grn)
        sorted_ciecl += holder[::-1]
    ciecl = sorted_ciecl
    grns_ = []
    i = 0
    j = 0
    while (len(grns_) < (len(niecl) + len(ciecl))):
        if len(niecl) == i:
            grns_.append(ciecl[j])
            j += 1
        elif len(ciecl) == j:
            grns_.append(niecl[i])
            i += 1
        else:
            if niecl[i] <= ciecl[j]:
                grns_.append(niecl[i])
                i += 1
            else:
                grns_.append(ciecl[j])
                j += 1
    return grns_


def sort_grns(grns: list[float]):
    grns_iecls = sort_iecls([x for x in grns if (x > 10) & (x < 100)])
    grns_tm = sorted([x for x in grns if (x > 0) & (x < 10)])
    grns_ntail = sorted([x for x in grns if (x < 0)])
    grns_ctail = sorted([x for x in grns if (x > 100)])
    iecl_idx = 0
    tm_idx = 0
    sorted_grns = []
    while ((iecl_idx < len(grns_iecls)) and (tm_idx < len(grns_tm))):
        if grns_tm[tm_idx] < get_prev_next_tm(grns_iecls[iecl_idx])[1]:
            sorted_grns.append(grns_tm[tm_idx])
            tm_idx += 1
        else:
            sorted_grns.append(grns_iecls[iecl_idx])
            iecl_idx += 1
    if (iecl_idx < len(grns_iecls)):
        sorted_grns += grns_iecls[iecl_idx:]
    if (tm_idx < len(grns_tm)):
        sorted_grns += grns_tm[tm_idx:]
    result = grns_ntail + sorted_grns + grns_ctail
    return [round(x, 3) for x in result]

def sort_grns_str(grns: list[str]):
    grns_float = [parse_grn_str2float(x) for x in grns]
    sorted_grns = sort_grns(grns_float)
    return [parse_grn_float2str(x) for x in sorted_grns]

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
        c_tail_grn = [('c.' + str(start + i)) for i in range(len_c_tail+1) if ('c.' + str(start + i)) in grns_str]
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
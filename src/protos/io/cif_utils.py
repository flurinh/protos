"""
CIF file utilities for Protos.

This module provides functions for parsing and generating mmCIF/PDBx files.
It allows conversion between pandas DataFrames and CIF format.
"""

import os
import re
import random
import logging
import pandas as pd
import numpy as np
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Tuple, Union, Optional, Any, Set

# Configure logger
logger = logging.getLogger(__name__)

# Define standard amino acid codes for reference
STD_AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]

# Define CIF column mapping for standardization
CIF_COLUMN_MAPPING = {
    # mmCIF standard to internal representation
    'group_PDB': 'group',
    'id': 'atom_id',
    'type_symbol': 'element',
    'label_atom_id': 'atom_name',
    'label_alt_id': 'alt_id',
    'label_comp_id': 'res_name',
    'label_asym_id': 'label_chain_id',
    'label_entity_id': 'entity_id',
    'label_seq_id': 'label_seq_id',
    'pdbx_PDB_ins_code': 'insertion',
    'Cartn_x': 'x',
    'Cartn_y': 'y',
    'Cartn_z': 'z',
    'occupancy': 'occupancy',
    'B_iso_or_equiv': 'b_factor',
    'pdbx_formal_charge': 'charge',
    'auth_seq_id': 'auth_seq_id',
    'auth_comp_id': 'auth_comp_id',
    'auth_asym_id': 'auth_chain_id',
    'auth_atom_id': 'auth_atom_name',
    'pdbx_PDB_model_num': 'model_num'
}

# Reverse mapping for writing
REVERSE_COLUMN_MAPPING = {v: k for k, v in CIF_COLUMN_MAPPING.items()}

# Define required and optional columns for a valid CIF DataFrame
REQUIRED_COLUMNS = [
    'group', 'atom_id', 'atom_name', 'res_name', 'auth_chain_id',
    'auth_seq_id', 'x', 'y', 'z'
]

OPTIONAL_COLUMNS = [
    'element', 'alt_id', 'label_chain_id', 'entity_id', 'label_seq_id',
    'insertion', 'occupancy', 'b_factor', 'charge', 'auth_comp_id',
    'auth_atom_name', 'model_num', 'pdb_id', 'res_name3l', 'res_name1l',
    'res_id', 'gen_seq_id'
]

# Define standard fixed-width CIF format positions
CIF_PIVOTS = {
    "group": 0,  # ATOM/HETATM
    "atom_id": 7,  # e.g., 1869
    "type_symbol": 12,  # e.g., O
    "atom_name": 15,  # e.g., O5 or "C1'"
    "alt_id": 21,  # e.g., .
    "res_name": 23,  # e.g., BOG
    "chain_id": 27,  # e.g., E
    "entity_id": 29,  # e.g., 4
    "label_seq_id": 31,  # e.g., .
    "ins_code": 35,  # e.g., ?
    "x_coord": 37,  # e.g., -24.610
    "y_coord": 45,  # e.g., -14.819
    "z_coord": 52,  # e.g., 10.118
    "occupancy": 60,  # e.g., 0.93
    "b_factor": 65,  # e.g., 54.64
    "formal_charge": 72,  # e.g., ?
    "auth_seq_id": 74,  # e.g., 1265
    "auth_comp_id": 79,  # e.g., BOG
    "auth_asym_id": 83,  # e.g., A
    "auth_atom_id": 85,  # e.g., O5 or "C1'"
    "model_num": 91  # e.g., 1
}


def cif_to_df(cif_content: Union[str, bytes], pdb_id: Optional[str] = None) -> pd.DataFrame:
    """
    Parse CIF format content into a pandas DataFrame.

    Args:
        cif_content: String or bytes content of a CIF file
        pdb_id: Optional PDB ID to assign (extracted from filename if None)

    Returns:
        DataFrame with standardized column names and data types

    Raises:
        ValueError: If CIF parsing fails
    """
    # Initialize data dictionary
    data = defaultdict(list)

    # Extract data from CIF content
    try:
        # Ensure content is string, not bytes
        if isinstance(cif_content, bytes):
            cif_content = cif_content.decode('utf-8')

        # Find the loop_ section containing atom_site data
        # Try to match a specific pattern for atom_site loop section
        loop_match = re.search(r'loop_\s+((?:_atom_site\.[^\n]+\s+)+)(.*?)(?=\s*#\s*$|\s*loop_|\s*data_|\Z)',
                               cif_content, re.DOTALL)

        # If that fails, try a more general pattern looking for ATOM/HETATM
        if not loop_match:
            logger.debug("Trying alternative loop pattern match")
            loop_match = re.search(
                r'loop_\s+((?:_atom_site\.[^\n]+\s+)+)((?:ATOM|HETATM).*?)(?=\s*#|\s*loop_|\s*data_|\Z)',
                cif_content, re.DOTALL)

        if not loop_match:
            raise ValueError("No atom_site loop found in CIF content")

        # Extract column definitions and data section
        column_defs = loop_match.group(1).strip().split('\n')
        data_section = loop_match.group(2).strip()

        # Clean up column names
        columns = []
        for col_def in column_defs:
            col_name = col_def.strip().replace('_atom_site.', '')
            columns.append(col_name)

        # Process data lines
        lines = data_section.split('\n')
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Parse line using simplified method
            values = _parse_fixed_width_line(line, columns)

            # Map values to standard column names
            for i, col in enumerate(columns):
                if i < len(values):
                    mapped_col = CIF_COLUMN_MAPPING.get(col, col)
                    value = values[i]
                    # Convert '.' and '?' to None
                    if value in ['.', '?']:
                        value = None
                    data[mapped_col].append(value)

            # Add PDB ID if provided or extract from data_block
            if pdb_id is not None:
                data['pdb_id'].append(pdb_id)
            elif 'pdb_id' not in data:
                # Try to extract PDB ID from data_block
                data_block_match = re.search(r'data_(\w+)', cif_content)
                if data_block_match:
                    extracted_pdb_id = data_block_match.group(1).lower()
                    data['pdb_id'].append(extracted_pdb_id)
                else:
                    data['pdb_id'].append('UNKNOWN')

            # Track information needed for res_id
            res_name = None
            auth_seq_id = None
            chain_id = None

            for i, col in enumerate(columns):
                if col == 'label_comp_id' and i < len(values):
                    res_name = values[i]
                elif col == 'auth_seq_id' and i < len(values):
                    auth_seq_id = values[i]
                elif col == 'auth_asym_id' and i < len(values):
                    chain_id = values[i]

            # Create res_id (resname_authseqid_chain)
            if res_name and auth_seq_id and chain_id:
                res_id = f"{res_name}_{auth_seq_id}_{chain_id}"
                data['res_id'].append(res_id)
            else:
                data['res_id'].append(None)

            # Add residue name in different formats
            if res_name:
                data['res_name3l'].append(res_name)

                # Convert 3-letter code to 1-letter code if possible
                aa_map = {
                    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                }
                data['res_name1l'].append(aa_map.get(res_name, 'X'))
            else:
                data['res_name3l'].append(None)
                data['res_name1l'].append('X')

            # Add gen_seq_id (sequential numbering)
            data['gen_seq_id'].append(len(data['gen_seq_id']) + 1)

        # Create DataFrame
        df = pd.DataFrame(data)

        # Convert numeric columns
        numeric_cols = ['atom_id', 'auth_seq_id', 'label_seq_id', 'x', 'y', 'z',
                        'occupancy', 'b_factor', 'model_num', 'gen_seq_id']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # Make sure we have the key columns from standard mapping
        if 'atom_name' not in df.columns and 'label_atom_id' in df.columns:
            df['atom_name'] = df['label_atom_id']

        if 'res_name' not in df.columns and 'label_comp_id' in df.columns:
            df['res_name'] = df['label_comp_id']

        # Fill any missing columns with reasonable defaults
        if 'group' not in df.columns:
            df['group'] = 'ATOM'

        # Make sure res_id is properly constructed if not already present
        if 'res_id' not in df.columns:
            # Generate res_id from component columns
            df['res_id'] = df.apply(
                lambda row: f"{row['res_name']}_{row['auth_seq_id']}_{row['auth_chain_id']}"
                if pd.notna(row['res_name']) and pd.notna(row['auth_seq_id']) and pd.notna(row['auth_chain_id'])
                else None,
                axis=1
            )

        return df

    except Exception as e:
        logger.error(f"Error parsing CIF content: {e}")
        raise ValueError(f"Failed to parse CIF content: {e}")

def read_cif_file(file_path: str) -> pd.DataFrame:
    """
    Read a CIF file and convert to a DataFrame.
    
    Args:
        file_path: Path to the CIF file
        
    Returns:
        DataFrame with standardized column names
    
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    # Extract PDB ID from filename
    pdb_id = os.path.basename(file_path).split('.')[0].lower()
    if '_v' in pdb_id:  # Handle versioned filenames
        pdb_id = pdb_id.split('_v')[0]
    
    # Read the file and parse content
    with open(file_path, 'r') as f:
        cif_content = f.read()
        
    return cif_to_df(cif_content, pdb_id)


def df_to_cif(df: pd.DataFrame, pdb_id: Optional[str] = None) -> str:
    """
    Convert a DataFrame to CIF format string using a fixed width line.

    Args:
        df: DataFrame with atomic structure data
        pdb_id: Optional PDB ID to use (extracted from df if not provided)

    Returns:
        String with CIF format content

    Raises:
        ValueError: If required columns are missing
    """
    # Validate required columns
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Use provided PDB ID or extract from DataFrame
    if pdb_id is None:
        if 'pdb_id' in df.columns:
            pdb_id = df['pdb_id'].iloc[0]
        else:
            pdb_id = "UNKNOWN"

    # Generate default values for optional columns
    df = df.copy()  # Create a copy to avoid modifying the original

    # Set defaults for optional columns
    if 'b_factor' not in df.columns and 'B_iso_or_equiv' not in df.columns:
        df['b_factor'] = 30.0

    if 'occupancy' not in df.columns:
        df['occupancy'] = 1.00

    if 'element' not in df.columns:
        # Try to infer element from atom name
        df['element'] = df['atom_name'].apply(
            lambda x: x[0] if isinstance(x, str) and len(x) > 0 else 'C'
        )

    if 'model_num' not in df.columns:
        df['model_num'] = 1

    # Ensure label_seq_id is present
    if 'label_seq_id' not in df.columns:
        if 'auth_seq_id' in df.columns:
            # Get unique residues
            unique_residues = df[['auth_seq_id', 'res_name']].drop_duplicates()

            # Create a sequential numbering map
            label_seq_map = {}
            for i, (_, row) in enumerate(unique_residues.iterrows()):
                label_seq_map[(row['auth_seq_id'], row['res_name'])] = i + 1

            # Apply mapping to create label_seq_id
            df['label_seq_id'] = df.apply(
                lambda row: label_seq_map.get((row['auth_seq_id'], row['res_name']), 1),
                axis=1
            )

    # Generate CIF content
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cif_content = f"data_{pdb_id}\n"
    cif_content += f"#\n# Created on {timestamp}\n#\n"
    cif_content += f"_entry.id {pdb_id}\n"

    # Basic structure metadata
    cif_content += "_struct.title 'Generated Structure'\n"

    # Write atom_site category with proper loop_ format
    cif_content += "\nloop_\n"
    cif_content += "_atom_site.group_PDB \n"
    cif_content += "_atom_site.id \n"
    cif_content += "_atom_site.type_symbol \n"
    cif_content += "_atom_site.label_atom_id \n"
    cif_content += "_atom_site.label_alt_id \n"
    cif_content += "_atom_site.label_comp_id \n"
    cif_content += "_atom_site.label_asym_id \n"
    cif_content += "_atom_site.label_entity_id \n"
    cif_content += "_atom_site.label_seq_id \n"
    cif_content += "_atom_site.pdbx_PDB_ins_code \n"
    cif_content += "_atom_site.Cartn_x \n"
    cif_content += "_atom_site.Cartn_y \n"
    cif_content += "_atom_site.Cartn_z \n"
    cif_content += "_atom_site.occupancy \n"
    cif_content += "_atom_site.B_iso_or_equiv \n"
    cif_content += "_atom_site.pdbx_formal_charge \n"
    cif_content += "_atom_site.auth_seq_id \n"
    cif_content += "_atom_site.auth_comp_id \n"
    cif_content += "_atom_site.auth_asym_id \n"
    cif_content += "_atom_site.auth_atom_id \n"
    cif_content += "_atom_site.pdbx_PDB_model_num \n"

    # Define standard amino acid codes for reference
    std_amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    ]

    # Helper function to place a field value at a specific position in the line
    def place_at_position(line_chars, position, value):
        # Convert to list for easier character replacement
        for i, char in enumerate(value):
            if position + i < len(line_chars):
                line_chars[position + i] = char
        return line_chars

    # Variable to store the first row's field values for debugging
    first_row_values = None

    # Write atom records using fixed width lines
    for i, (_, row) in enumerate(df.iterrows()):
        try:
            # Initialize a line with 92 spaces
            line_chars = [' '] * 92

            # Extract and safely convert values
            group = row.get('group', 'ATOM')

            try:
                atom_id = str(int(row['atom_id'])) if pd.notna(row['atom_id']) else "1"
            except (ValueError, TypeError):
                atom_id = "1"

            element = str(row.get('element', 'C')) if pd.notna(row.get('element')) else 'C'
            atom_name = str(row['atom_name']) if pd.notna(row['atom_name']) else 'X'
            alt_id = str(row.get('alt_id', '.')) if pd.notna(row.get('alt_id', '.')) else '.'
            res_name = str(row['res_name']) if pd.notna(row['res_name']) else 'UNK'

            # Determine if this is a protein residue
            is_protein = res_name in std_amino_acids

            auth_chain_id = str(row['auth_chain_id']) if pd.notna(row['auth_chain_id']) else 'A'
            display_chain_id = auth_chain_id if is_protein else "R"
            entity_id = "1" if is_protein else "4"

            try:
                label_seq_id = str(int(row['label_seq_id'])) if pd.notna(row['label_seq_id']) else "1"
            except (ValueError, TypeError):
                label_seq_id = "1"

            try:
                auth_seq_id = str(int(row['auth_seq_id'])) if pd.notna(row['auth_seq_id']) else "1"
            except (ValueError, TypeError):
                auth_seq_id = "1"

            # Format coordinates with 3 decimal places
            x_coord = f"{float(row['x']):.3f}" if pd.notna(row['x']) else "0.000"
            y_coord = f"{float(row['y']):.3f}" if pd.notna(row['y']) else "0.000"
            z_coord = f"{float(row['z']):.3f}" if pd.notna(row['z']) else "0.000"

            # Handle occupancy and b_factor
            if res_name == "HOH":
                occupancy = "1.00"
            else:
                occupancy = "0.93"

            if 'b_factor' in row and pd.notna(row['b_factor']):
                b_factor = f"{float(row['b_factor']):.2f}"
            elif 'B_iso_or_equiv' in row and pd.notna(row['B_iso_or_equiv']):
                b_factor = f"{float(row['B_iso_or_equiv']):.2f}"
            else:
                b_factor = "30.00"

            try:
                model_num = str(int(row.get('model_num', 1))) if pd.notna(row.get('model_num', 1)) else "1"
            except (ValueError, TypeError):
                model_num = "1"

            # Save field values from the first row for debugging
            if i == 0:
                first_row_values = {
                    "group": group,
                    "atom_id": atom_id,
                    "type_symbol": element,
                    "atom_name": atom_name,
                    "alt_id": alt_id,
                    "res_name": res_name,
                    "chain_id": display_chain_id,
                    "entity_id": entity_id,
                    "label_seq_id": label_seq_id,
                    "ins_code": "?",
                    "x_coord": x_coord,
                    "y_coord": y_coord,
                    "z_coord": z_coord,
                    "occupancy": occupancy,
                    "b_factor": b_factor,
                    "formal_charge": "?",
                    "auth_seq_id": auth_seq_id,
                    "auth_comp_id": res_name,
                    "auth_asym_id": display_chain_id,
                    "auth_atom_id": atom_name,
                    "model_num": model_num
                }

                # Print field values by position
                print("Field values for first residue by position:\n")
                print("{:<5} {:<15} {:<15} {:<5}".format("Pos", "Field", "Value", "End"))
                print("=" * 45)

                for field, position in sorted(CIF_PIVOTS.items(), key=lambda x: x[1]):
                    value = first_row_values.get(field, "?")
                    end_position = position + len(value) - 1
                    print("{:<5} {:<15} {:<15} {:<5}".format(
                        position, field, value, end_position
                    ))

            # Place each field value at its exact position
            line_chars = place_at_position(line_chars, CIF_PIVOTS["group"], group)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["atom_id"], atom_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["type_symbol"], element)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["atom_name"], atom_name)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["alt_id"], alt_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["res_name"], res_name)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["chain_id"], display_chain_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["entity_id"], entity_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["label_seq_id"], label_seq_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["ins_code"], "?")
            line_chars = place_at_position(line_chars, CIF_PIVOTS["x_coord"], x_coord)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["y_coord"], y_coord)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["z_coord"], z_coord)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["occupancy"], occupancy)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["b_factor"], b_factor)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["formal_charge"], "?")
            line_chars = place_at_position(line_chars, CIF_PIVOTS["auth_seq_id"], auth_seq_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["auth_comp_id"], res_name)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["auth_asym_id"], display_chain_id)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["auth_atom_id"], atom_name)
            line_chars = place_at_position(line_chars, CIF_PIVOTS["model_num"], model_num)

            # Convert character list back to string
            line = ''.join(line_chars)
            cif_content += line + "\n"

        except Exception as e:
            logger.warning(f"Error formatting atom record: {e}")
            continue

    # Proper termination of the file
    cif_content += "#\n"

    return cif_content


def write_cif_file(file_path: str, df: pd.DataFrame,
                   versioned: bool = False, force_overwrite: bool = False) -> str:
    """
    Write a DataFrame to a CIF file.

    Args:
        file_path: Path to save the CIF file
        df: DataFrame with atomic structure data
        versioned: Whether to add version numbering (_v1, _v2, etc.)
        force_overwrite: Whether to allow overwriting existing files

    Returns:
        Path to the written file

    Raises:
        ValueError: If required data is missing
        FileExistsError: If file exists and force_overwrite=False
    """
    import shutil
    from pathlib import Path

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)

    # Convert to Path object for reliable path handling
    file_path = Path(file_path)

    # Handle versioning if requested
    final_path = file_path
    if versioned:
        # Get base path without extension
        base_path = file_path.parent / file_path.stem
        extension = file_path.suffix

        # Find next available version number
        version = 1
        while True:
            versioned_path = f"{base_path}_v{version}{extension}"
            if not os.path.exists(versioned_path):
                final_path = Path(versioned_path)
                break
            version += 1
    else:
        # Check if file exists and handle force_overwrite
        if file_path.exists() and not force_overwrite:
            raise FileExistsError(
                f"File {file_path} already exists. Use versioned=True or force_overwrite=True."
            )

        # Create backup if overwriting
        if file_path.exists() and force_overwrite:
            backup_path = f"{file_path}.backup"
            shutil.copy2(file_path, backup_path)
            logger.info(f"Created backup at {backup_path}")

    # Extract PDB ID from file path if not in DataFrame
    pdb_id = None
    if 'pdb_id' not in df.columns:
        pdb_id = file_path.stem.lower()
        if '_v' in pdb_id:  # Handle versioned filenames
            pdb_id = pdb_id.split('_v')[0]

    # Generate CIF content
    cif_content = df_to_cif(df, pdb_id)

    # Write to file - convert Path to string for compatibility
    with open(str(final_path), 'w') as f:
        f.write(cif_content)

    # Return the exact path of the written file as a string
    return str(final_path)


def validate_cif_data(df: pd.DataFrame) -> dict:
    """
    Validate structure data before writing to CIF.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        Dictionary with validation results and issues
    """
    validation_results = {
        "atom_count": len(df),
        "issues": [],
        "is_valid": True
    }
    
    # Check for missing required columns
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        validation_results["issues"].append({
            "type": "missing_required_columns",
            "details": f"Missing required columns: {missing_cols}",
            "severity": "error"
        })
        validation_results["is_valid"] = False
    
    # Check for missing coordinates
    missing_coords = df[df[['x', 'y', 'z']].isna().any(axis=1)]
    if not missing_coords.empty:
        validation_results["issues"].append({
            "type": "missing_coordinates",
            "count": len(missing_coords),
            "affected_rows": missing_coords.index.tolist(),
            "details": "Some atoms have missing x, y, or z coordinates",
            "severity": "error"
        })
        validation_results["is_valid"] = False
        
    # Check for missing atom IDs
    missing_atom_ids = df[df['atom_id'].isna()]
    if not missing_atom_ids.empty:
        validation_results["issues"].append({
            "type": "missing_atom_ids",
            "count": len(missing_atom_ids),
            "affected_rows": missing_atom_ids.index.tolist(),
            "details": "Some atoms have missing atom IDs",
            "severity": "warning"
        })
        
    # Check for missing chain IDs
    missing_chain_ids = df[df['auth_chain_id'].isna()]
    if not missing_chain_ids.empty:
        validation_results["issues"].append({
            "type": "missing_chain_ids",
            "count": len(missing_chain_ids),
            "affected_rows": missing_chain_ids.index.tolist(),
            "details": "Some atoms have missing chain IDs",
            "severity": "warning"
        })
    
    # Check for valid occupancy values
    if 'occupancy' in df.columns:
        invalid_occupancy = df[(df['occupancy'] < 0) | (df['occupancy'] > 1)]
        if not invalid_occupancy.empty:
            validation_results["issues"].append({
                "type": "invalid_occupancy",
                "count": len(invalid_occupancy),
                "affected_rows": invalid_occupancy.index.tolist(),
                "details": "Occupancy values should be between 0 and 1",
                "severity": "warning"
            })
    
    # Check for residue continuity
    if 'auth_seq_id' in df.columns and 'auth_chain_id' in df.columns:
        for chain in df['auth_chain_id'].unique():
            chain_df = df[df['auth_chain_id'] == chain]
            seq_ids = sorted(chain_df['auth_seq_id'].unique())
            if len(seq_ids) > 1:
                gaps = []
                for i in range(len(seq_ids) - 1):
                    if seq_ids[i+1] - seq_ids[i] > 1:
                        gaps.append((seq_ids[i], seq_ids[i+1]))
                
                if gaps:
                    validation_results["issues"].append({
                        "type": "residue_gaps",
                        "chain": chain,
                        "gaps": gaps,
                        "details": f"Chain {chain} has residue numbering gaps",
                        "severity": "info"
                    })
    
    return validation_results

def fix_cif_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fix common issues in CIF data.
    
    Args:
        df: DataFrame with atomic structure data
    
    Returns:
        Corrected DataFrame
    """
    # Create a copy to avoid modifying the original
    fixed_df = df.copy()
    
    # Generate sequential atom IDs if missing
    if 'atom_id' not in fixed_df.columns or fixed_df['atom_id'].isna().any():
        fixed_df['atom_id'] = range(1, len(fixed_df) + 1)
    
    # Set default chain ID if missing
    if 'auth_chain_id' not in fixed_df.columns or fixed_df['auth_chain_id'].isna().any():
        fixed_df['auth_chain_id'] = fixed_df['auth_chain_id'].fillna('A')
    
    # Set default group type if missing
    if 'group' not in fixed_df.columns or fixed_df['group'].isna().any():
        fixed_df['group'] = fixed_df['group'].fillna('ATOM')
    
    # Fix coordinates to ensure they're numeric
    for col in ['x', 'y', 'z']:
        if col in fixed_df.columns:
            fixed_df[col] = pd.to_numeric(fixed_df[col], errors='coerce')
            # Replace NaN values with zeros
            fixed_df[col].fillna(0.0, inplace=True)
    
    # Set default occupancy if missing or invalid
    if 'occupancy' not in fixed_df.columns:
        fixed_df['occupancy'] = 1.0
    else:
        # Ensure occupancy is numeric and within valid range
        fixed_df['occupancy'] = pd.to_numeric(fixed_df['occupancy'], errors='coerce')
        fixed_df['occupancy'].fillna(1.0, inplace=True)
        fixed_df.loc[fixed_df['occupancy'] < 0, 'occupancy'] = 0.0
        fixed_df.loc[fixed_df['occupancy'] > 1, 'occupancy'] = 1.0
    
    # Set default model number if missing
    if 'model_num' not in fixed_df.columns or fixed_df['model_num'].isna().any():
        fixed_df['model_num'] = 1
    
    # Ensure res_name is populated
    if 'res_name' not in fixed_df.columns and 'res_name3l' in fixed_df.columns:
        fixed_df['res_name'] = fixed_df['res_name3l']
    elif 'res_name' in fixed_df.columns and fixed_df['res_name'].isna().any():
        fixed_df['res_name'].fillna('UNK', inplace=True)
    
    # Generate res_id if missing
    if 'res_id' not in fixed_df.columns or fixed_df['res_id'].isna().any():
        fixed_df['res_id'] = fixed_df.apply(
            lambda row: f"{row['res_name']}_{row['auth_seq_id']}_{row['auth_chain_id']}" 
            if pd.notna(row['res_name']) and pd.notna(row['auth_seq_id']) and pd.notna(row['auth_chain_id']) 
            else None, 
            axis=1
        )
    
    # Generate res_name3l if missing
    if 'res_name3l' not in fixed_df.columns and 'res_name' in fixed_df.columns:
        fixed_df['res_name3l'] = fixed_df['res_name']
    
    # Generate res_name1l if missing
    if 'res_name1l' not in fixed_df.columns and 'res_name3l' in fixed_df.columns:
        aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
            'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
            'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        fixed_df['res_name1l'] = fixed_df['res_name3l'].map(lambda x: aa_map.get(x, 'X') if pd.notna(x) else 'X')
    
    return fixed_df

def merge_structures(dfs: List[pd.DataFrame], 
                    chain_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Merge multiple structure DataFrames into a single structure.
    
    Args:
        dfs: List of DataFrames to merge
        chain_mapping: Optional mapping to reassign chain IDs
    
    Returns:
        DataFrame with merged structure data
    """
    if not dfs:
        raise ValueError("No DataFrames provided for merging")
    
    # Create copies to avoid modifying originals
    copies = [df.copy() for df in dfs]
    
    # Assign new chain IDs if mapping provided
    if chain_mapping:
        for i, df in enumerate(copies):
            if 'auth_chain_id' in df.columns:
                # Apply chain mapping or keep original
                df['auth_chain_id'] = df['auth_chain_id'].map(
                    lambda x: chain_mapping.get(x, x) if pd.notna(x) else 'X'
                )
                
                # Update label_chain_id if present
                if 'label_chain_id' in df.columns:
                    df['label_chain_id'] = df['auth_chain_id']
                
                # Update res_id to reflect new chain
                if 'res_id' in df.columns:
                    df['res_id'] = df.apply(
                        lambda row: f"{row['res_name']}_{row['auth_seq_id']}_{row['auth_chain_id']}" 
                        if pd.notna(row['res_name']) and pd.notna(row['auth_seq_id']) and pd.notna(row['auth_chain_id']) 
                        else None, 
                        axis=1
                    )
    
    # Concatenate DataFrames
    merged_df = pd.concat(copies, ignore_index=True)
    
    # Reassign atom IDs sequentially
    merged_df['atom_id'] = range(1, len(merged_df) + 1)
    
    # Regenerate gen_seq_id sequentially
    if 'gen_seq_id' in merged_df.columns:
        merged_df['gen_seq_id'] = range(1, len(merged_df) + 1)
    
    return merged_df

def extract_bioassembly(df: pd.DataFrame, 
                       assembly_id: int = 1, 
                       include_chains: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Extract biological assembly from a structure DataFrame.
    
    Args:
        df: DataFrame with atomic structure data
        assembly_id: ID of the biological assembly to extract
        include_chains: Optional list of chains to include
    
    Returns:
        DataFrame with biological assembly structure
    """
    # Create a copy to avoid modifying the original
    result_df = df.copy()
    
    # Filter by specified chains if provided
    if include_chains:
        result_df = result_df[result_df['auth_chain_id'].isin(include_chains)]
    
    # Filter to keep only protein, ligand, and important heteroatoms
    if 'group' in result_df.columns:
        # Keep all ATOM records
        is_atom = result_df['group'] == 'ATOM'
        
        # Keep HETATM records that are not water
        is_important_hetatm = (result_df['group'] == 'HETATM') & (result_df['res_name'] != 'HOH')
        
        # Combine filters
        result_df = result_df[is_atom | is_important_hetatm]
    
    return result_df


def _parse_fixed_width_line(line: str, columns: list) -> list:
    """
    Parse a CIF line into values using simple whitespace splitting.

    Args:
        line: The line to parse
        columns: List of column names (for reference)

    Returns:
        List of values extracted from the line
    """
    import re

    # Split by one or more whitespace characters
    values = re.split(r'\s+', line.strip())

    # Ensure we have enough values to match columns
    if len(values) < len(columns):
        values.extend([None] * (len(columns) - len(values)))
    elif len(values) > len(columns):
        # If we have more values than columns, truncate
        values = values[:len(columns)]

    # Convert '.' and '?' to None
    values = [None if v in ['.', '?'] else v for v in values]

    return values

def _add_field(current_line: str, position: int, value: str) -> str:
    """
    Add a field at a specific position in a line.
    
    Args:
        current_line: The current line being built
        position: The position to place the field
        value: The value to add
        
    Returns:
        Updated line with the new field added
    """
    # Pad current line with spaces if needed
    while len(current_line) < position:
        current_line += " "
    return current_line + value

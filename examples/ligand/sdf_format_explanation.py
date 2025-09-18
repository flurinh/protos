"""
SDF Format Explanation and Parser

This script explains the SDF file format and demonstrates parsing the coordinate block.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

def explain_sdf_format():
    """Explain the SDF file format structure."""
    
    print("=== SDF (Structure Data Format) File Structure ===\n")
    
    print("1. HEADER BLOCK (3 lines):")
    print("   Line 1: Molecule name")
    print("   Line 2: Program/user info (often blank)")
    print("   Line 3: Comment (often blank)")
    print()
    
    print("2. CONNECTION TABLE BLOCK:")
    print("   Line 4: Counts line")
    print("   Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv")
    print("   - aaa: number of atoms")
    print("   - bbb: number of bonds")
    print("   - Other fields for properties, version, etc.")
    print()
    
    print("3. ATOM BLOCK (one line per atom):")
    print("   Columns:")
    print("   1-10:   x coordinate (10.4 format)")
    print("   11-20:  y coordinate (10.4 format)")
    print("   21-30:  z coordinate (10.4 format)")
    print("   31-33:  atom symbol (left justified)")
    print("   34-35:  mass difference")
    print("   36-38:  charge")
    print("   39-41:  atom stereo parity")
    print("   42-44:  hydrogen count")
    print("   45-47:  stereo care box")
    print("   48-50:  valence")
    print("   51-63:  not used")
    print("   64-66:  atom-atom mapping number")
    print()
    
    print("4. BOND BLOCK (one line per bond):")
    print("   Columns:")
    print("   1-3:    first atom number")
    print("   4-6:    second atom number")
    print("   7-9:    bond type (1=single, 2=double, 3=triple, 4=aromatic)")
    print("   10-12:  stereo configuration")
    print("   13-15:  not used")
    print("   16-18:  bond topology")
    print("   19-21:  reacting center status")
    print()
    
    print("5. PROPERTIES BLOCK:")
    print("   > <PropertyName>")
    print("   Property value")
    print("   (blank line)")
    print()
    
    print("6. END:")
    print("   $$$$")


def parse_sdf_atom_line(line):
    """
    Parse a single atom line from SDF format.
    
    Args:
        line: Atom line from SDF file
        
    Returns:
        Dictionary with parsed atom data
    """
    # SDF atom line format (V2000):
    # Columns 1-30: coordinates (3 x F10.4)
    # Columns 31-33: atom symbol
    # Remaining columns: various atom properties
    
    atom_data = {
        'x': float(line[0:10].strip()),
        'y': float(line[10:20].strip()),
        'z': float(line[20:30].strip()),
        'element': line[31:34].strip(),
        'mass_diff': line[34:36].strip(),
        'charge': line[36:39].strip(),
        'stereo_parity': line[39:42].strip(),
        'hydrogen_count': line[42:45].strip(),
        'stereo_care': line[45:48].strip(),
        'valence': line[48:51].strip(),
        'atom_mapping': line[63:66].strip() if len(line) >= 66 else ''
    }
    
    return atom_data


def demonstrate_sdf_parsing():
    """Demonstrate parsing an actual SDF file."""
    
    print("\n\n=== SDF Parsing Example ===\n")
    
    # Create a simple example SDF content
    example_sdf = """Aspirin
  ChEMBL  08052024 2D

 13 13  0  0  0  0  0  0  0  0999 V2000
    2.5000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7500    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2500    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    2.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7500    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5000    2.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2500    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  2  0  0  0  0
  2 12  1  0  0  0  0
 12 13  1  0  0  0  0
M  END
> <CHEMBL_ID>
CHEMBL25

> <MW>
180.16

$$$$"""
    
    print("Example SDF content for Aspirin:")
    print("-" * 60)
    
    lines = example_sdf.strip().split('\n')
    
    # Parse header
    print(f"Molecule name: {lines[0]}")
    print(f"Program info: {lines[1]}")
    print(f"Comment: {lines[2]}")
    print()
    
    # Parse counts line
    counts_line = lines[3]
    num_atoms = int(counts_line[0:3])
    num_bonds = int(counts_line[3:6])
    print(f"Counts: {num_atoms} atoms, {num_bonds} bonds")
    print()
    
    # Parse atom block
    print("Atom Block (parsed):")
    print("Index | X      | Y      | Z      | Element | Charge")
    print("-" * 60)
    
    for i in range(num_atoms):
        atom_line = lines[4 + i]
        atom = parse_sdf_atom_line(atom_line)
        print(f"{i+1:5d} | {atom['x']:6.4f} | {atom['y']:6.4f} | {atom['z']:6.4f} | {atom['element']:7s} | {atom['charge']:6s}")
    
    print()
    
    # Show bond block
    print("Bond Block:")
    print("Atom1 | Atom2 | Type | Stereo")
    print("-" * 40)
    
    for i in range(num_bonds):
        bond_line = lines[4 + num_atoms + i]
        atom1 = int(bond_line[0:3])
        atom2 = int(bond_line[3:6])
        bond_type = int(bond_line[6:9])
        stereo = bond_line[9:12].strip()
        
        bond_name = {1: "Single", 2: "Double", 3: "Triple", 4: "Aromatic"}.get(bond_type, "Unknown")
        print(f"{atom1:5d} | {atom2:5d} | {bond_name:6s} | {stereo:6s}")
    
    print()
    
    # Show properties
    print("Properties:")
    in_property = False
    prop_name = ""
    for line in lines[4 + num_atoms + num_bonds:]:
        if line.startswith("> <"):
            prop_name = line[3:-1]
            in_property = True
        elif in_property and line.strip():
            print(f"  {prop_name}: {line}")
            in_property = False


def show_rdkit_access():
    """Show how to access SDF data using RDKit."""
    
    print("\n\n=== Accessing SDF Data with RDKit ===\n")
    
    print("Using RDKit to read SDF files:")
    print("```python")
    print("from rdkit import Chem")
    print()
    print("# Read SDF file")
    print("suppl = Chem.SDMolSupplier('molecules.sdf')")
    print()
    print("for mol in suppl:")
    print("    if mol is None:")
    print("        continue")
    print("    ")
    print("    # Get conformer (3D coordinates)")
    print("    conf = mol.GetConformer()")
    print("    ")
    print("    # Access atom coordinates")
    print("    for atom_idx in range(mol.GetNumAtoms()):")
    print("        atom = mol.GetAtomWithIdx(atom_idx)")
    print("        pos = conf.GetAtomPosition(atom_idx)")
    print("        ")
    print("        print(f'Atom {atom_idx}: {atom.GetSymbol()}')")
    print("        print(f'  X: {pos.x:.4f}')")
    print("        print(f'  Y: {pos.y:.4f}')")
    print("        print(f'  Z: {pos.z:.4f}')")
    print("    ")
    print("    # Access properties")
    print("    for prop_name in mol.GetPropNames():")
    print("        print(f'{prop_name}: {mol.GetProp(prop_name)}')")
    print("```")


def show_dataframe_structure():
    """Show the DataFrame structure when SDF is converted."""
    
    print("\n\n=== SDF to DataFrame Structure ===\n")
    
    print("When converted to DataFrame using sdf_to_dataframe():")
    print()
    print("Standard columns:")
    print("- smiles: SMILES string representation")
    print("- index: Molecule index in SDF file")
    print("- mol: RDKit molecule object (optional)")
    print()
    print("Property columns (from SDF properties):")
    print("- Any properties stored in the SDF file")
    print("- Examples: CHEMBL_ID, MW, LogP, Activity_nM, etc.")
    print()
    print("Note: The 3D coordinates are NOT in the DataFrame.")
    print("They remain in the RDKit molecule object.")
    print()
    print("To access coordinates from DataFrame:")
    print("```python")
    print("df = sdf_to_dataframe('molecules.sdf', include_mol=True)")
    print()
    print("for idx, row in df.iterrows():")
    print("    mol = row['mol']")
    print("    conf = mol.GetConformer()")
    print("    ")
    print("    # Get all atom positions")
    print("    positions = []")
    print("    for i in range(mol.GetNumAtoms()):")
    print("        pos = conf.GetAtomPosition(i)")
    print("        positions.append([pos.x, pos.y, pos.z])")
    print("```")


def main():
    """Run all explanations."""
    explain_sdf_format()
    demonstrate_sdf_parsing()
    show_rdkit_access()
    show_dataframe_structure()
    
    print("\n✅ SDF format explanation completed!")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""Utilities for RFdiffusion2 input preparation."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Tuple


def cif_to_pdb_with_ori(
    cif_path: Path,
    pdb_path: Path,
    ori_coords: Optional[Tuple[float, float, float]] = None,
    chains_to_keep: Optional[List[str]] = None,
) -> Path:
    """Convert mmCIF file to PDB format with ORI HETATM for centering.

    Args:
        cif_path: Path to input CIF file
        pdb_path: Path to output PDB file
        ori_coords: Optional (x, y, z) coordinates for ORI. If None, uses center of mass.
        chains_to_keep: Optional list of chain IDs to include. If None, keeps all.

    Returns:
        Path to output PDB file
    """
    from Bio.PDB import MMCIFParser

    # Parse CIF
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(cif_path))

    # Collect atoms and calculate center of mass
    coords_all = []
    atoms_data = []  # (record_type, atom_name, alt_loc, res_name, chain, res_num, icode, x, y, z, occ, bfac, element)

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            if chains_to_keep and chain_id not in chains_to_keep:
                continue

            for residue in chain:
                res_id = residue.get_id()
                hetflag = res_id[0]
                res_num = res_id[1]
                icode = res_id[2] if res_id[2] != ' ' else ''
                res_name = residue.get_resname()

                # Determine record type
                record_type = "HETATM" if hetflag != ' ' else "ATOM"

                for atom in residue:
                    coord = atom.get_coord()
                    coords_all.append(coord)

                    # Get element from atom name if not set
                    element = atom.element if atom.element else atom.name[0]
                    if len(element) == 1:
                        element = f" {element}"

                    atoms_data.append({
                        'record': record_type,
                        'name': atom.name,
                        'alt_loc': atom.get_altloc(),
                        'res_name': res_name,
                        'chain': chain_id,
                        'res_num': res_num,
                        'icode': icode,
                        'x': coord[0],
                        'y': coord[1],
                        'z': coord[2],
                        'occ': atom.get_occupancy(),
                        'bfac': atom.get_bfactor(),
                        'element': element,
                    })

    # Calculate center of mass if not specified
    if ori_coords is None:
        coords_all = np.array(coords_all)
        ori_coords = tuple(coords_all.mean(axis=0))

    # Write PDB with proper formatting
    with open(pdb_path, 'w') as f:
        # Write ORI first
        ori_x, ori_y, ori_z = ori_coords
        f.write(f"HETATM    1  O   ORI Z   1    {ori_x:8.3f}{ori_y:8.3f}{ori_z:8.3f}  1.00  0.00           O\n")

        # Write all atoms with proper PDB format
        atom_num = 2
        for atom in atoms_data:
            # Format atom name (4 chars, left-justified if starts with element)
            name = atom['name']
            if len(name) < 4:
                if len(atom['element'].strip()) == 1:
                    name = f" {name:<3}"
                else:
                    name = f"{name:<4}"

            # PDB format:
            # ATOM   serial name altLoc resName chain resSeq iCode    x       y       z     occ   bfac          element
            # 1-6    7-11   13-16 17    18-20   22    23-26  27      31-38   39-46   47-54  55-60 61-66         77-78
            line = (
                f"{atom['record']:<6}"
                f"{atom_num:>5} "
                f"{name:4}"
                f"{atom['alt_loc']:1}"
                f"{atom['res_name']:>3} "
                f"{atom['chain']:1}"
                f"{atom['res_num']:>4}"
                f"{atom['icode']:1}   "
                f"{atom['x']:>8.3f}"
                f"{atom['y']:>8.3f}"
                f"{atom['z']:>8.3f}"
                f"{atom['occ']:>6.2f}"
                f"{atom['bfac']:>6.2f}"
                f"          "
                f"{atom['element']:>2}\n"
            )
            f.write(line)
            atom_num += 1

        f.write("END\n")

    return pdb_path


def get_triad_positions_from_cif(cif_path: Path) -> Dict[str, Dict]:
    """Extract catalytic triad positions from a rhodozyme placement CIF.

    Returns dict mapping residue type to {chain, resnum, atoms}.
    """
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(cif_path))

    # Look for catalytic residues (SER/CYS, HIS, ASP/GLU)
    # In the grafted triad, these will be non-standard residues or in chain B
    triad = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                resid = residue.get_id()[1]
                chain_id = chain.get_id()

                # Check for catalytic residue types
                if resname in ['SER', 'CYS'] and 'OG' in [a.name for a in residue.get_atoms()]:
                    if 'nucleophile' not in triad:
                        triad['nucleophile'] = {
                            'resname': resname,
                            'chain': chain_id,
                            'resnum': resid,
                            'key_atoms': ['OG', 'CB', 'CA'] if resname == 'SER' else ['SG', 'CB', 'CA']
                        }
                elif resname == 'HIS':
                    if 'histidine' not in triad:
                        triad['histidine'] = {
                            'resname': resname,
                            'chain': chain_id,
                            'resnum': resid,
                            'key_atoms': ['NE2', 'ND1', 'CB', 'CA']
                        }
                elif resname in ['ASP', 'GLU']:
                    if 'acidic' not in triad:
                        triad['acidic'] = {
                            'resname': resname,
                            'chain': chain_id,
                            'resnum': resid,
                            'key_atoms': ['OD1', 'OD2', 'CB', 'CA'] if resname == 'ASP' else ['OE1', 'OE2', 'CB', 'CA']
                        }

    return triad


def prepare_rfd2_input(
    placement_cif: Path,
    output_dir: Path,
    triad_positions: Optional[Dict[str, int]] = None,
    ligand_name: str = "RET",
) -> Dict:
    """Prepare input files and config for RFdiffusion2.

    Args:
        placement_cif: Path to placement CIF from Boltz/design
        output_dir: Directory for prepared files
        triad_positions: Optional manual triad positions {chain_resnum: residue_type}
        ligand_name: Name of ligand (default RET for retinal)

    Returns:
        Dict with 'pdb_path', 'contigs', 'contig_atoms', 'ligand'
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert CIF to PDB with ORI
    pdb_path = output_dir / f"{placement_cif.stem}.pdb"
    cif_to_pdb_with_ori(placement_cif, pdb_path)

    # Get structure info
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(placement_cif))

    # Count residues per chain
    chain_lengths = {}
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            residues = [r for r in chain if r.id[0] == ' ']  # Standard residues only
            if residues:
                chain_lengths[chain_id] = len(residues)

    # Build contigs for full structure (partial diffusion mode)
    # Format: "A1-N" for each chain
    contigs_parts = []
    for chain_id, length in sorted(chain_lengths.items()):
        contigs_parts.append(f"{chain_id}1-{length}")
    contigs = "/0 ".join(contigs_parts)  # Join chains with /0 for separate chains

    # Build contig_atoms for triad constraints if provided
    contig_atoms = None
    if triad_positions:
        contig_atoms = {}
        for res_key, atoms in triad_positions.items():
            contig_atoms[res_key] = atoms

    return {
        'pdb_path': pdb_path,
        'contigs': contigs,
        'contig_atoms': contig_atoms,
        'ligand': ligand_name,
        'chain_lengths': chain_lengths,
    }


if __name__ == "__main__":
    # Test with a placement
    import sys

    if len(sys.argv) > 1:
        cif_path = Path(sys.argv[1])
    else:
        cif_path = Path("/data/fast/projects/protos/data/models/boltzgen/20260202_rhodozyme/placement_00/placement_00.cif")

    if not cif_path.exists():
        print(f"CIF file not found: {cif_path}")
        sys.exit(1)

    output_dir = Path("/tmp/rfd2_test_prep")
    result = prepare_rfd2_input(cif_path, output_dir)

    print(f"Prepared input:")
    print(f"  PDB: {result['pdb_path']}")
    print(f"  Contigs: {result['contigs']}")
    print(f"  Chain lengths: {result['chain_lengths']}")
    print(f"  Ligand: {result['ligand']}")

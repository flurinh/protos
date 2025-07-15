#!/usr/bin/env python3
"""
GRN-Structure Integration Example
=================================

This script demonstrates how to integrate GRN (Generic Residue Numbering)
annotations with 3D protein structures for structural analysis.

Workflow:
1. Download/load microbial opsin structures from PDB
2. Extract sequences from structure data
3. Assign GRN numbers to sequences
4. Map GRN annotations back to structure residues
5. Analyze structure-GRN relationships (distances, interactions)

Key Features:
- Integrates StructureProcessor with GRNProcessor
- Maps sequence positions to 3D coordinates
- Calculates distances between key GRN positions
- Analyzes structural features by GRN position

Usage:
    python examples/grn_structure_integration.py
    
Output:
    - GRN-structure analysis saved to: grn_structure_analysis.csv
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.structure.structure_processor import StructureProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.sequence import (
    init_aligner, align_blosum62, format_alignment, mmseqs2_align2
)
from protos.processing.grn.grn_utils import (
    get_seq, sort_grns_str, GRNConfigManager
)
from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment, expand_annotation
)
from protos.loaders.download_structures import download_protein_structures

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GRNStructureIntegration:
    """Demonstrates GRN-Structure integration with real data."""
    
    def __init__(self):
        # Set up ProtosPaths to use the project's data directory
        from protos.io.paths.path_config import ProtosPaths
        
        # Use the project's data directory (non-deprecated approach)
        data_dir = Path(__file__).parent.parent / "data"  # Go up to project root
        paths = ProtosPaths(data_root=str(data_dir))
        
        # Initialize processors with zero configuration
        self.struct_proc = StructureProcessor(name="grn_struct_demo", paths=paths)
        self.grn_proc = GRNProcessor(name="grn_struct_demo", paths=paths)
        self.seq_proc = SequenceProcessor(name="grn_struct_demo", paths=paths)
        
        # Known microbial opsin structures
        self.opsin_structures = {
            '1M0K': 'Bacteriorhodopsin (dark state)',
            '1M0L': 'Bacteriorhodopsin (K intermediate)',
            '2NTU': 'Sensory rhodopsin II',
            '3DDL': 'Xanthorhodopsin',
            '4HYJ': 'Acetabularia rhodopsin',
            '5B6W': 'Channelrhodopsin chimera C1C2',
            '6EIG': 'Channelrhodopsin 2',
            '6CSM': 'Krokinobacter rhodopsin 2',
        }
        
    def download_structures(self):
        """Download microbial opsin structures."""
        logger.info("Downloading microbial opsin structures...")
        pdb_ids = list(self.opsin_structures.keys())
        
        try:
            downloaded = download_protein_structures(pdb_ids)
            logger.info(f"Downloaded {len(downloaded)} structures")
            return downloaded
        except Exception as e:
            logger.warning(f"Download failed: {e}")
            # Check which structures are already available
            available = []
            for pdb_id in pdb_ids:
                if self.struct_proc.entity_exists(pdb_id.lower()):
                    available.append(pdb_id.lower())
            logger.info(f"Found {len(available)} structures already available")
            return available
    
    def extract_sequences_from_structures(self, pdb_ids):
        """Extract sequences from structure data."""
        logger.info("\nExtracting sequences from structures...")
        sequences = {}
        
        for pdb_id in pdb_ids:
            try:
                # Load structure
                structure = self.struct_proc.load_entity(pdb_id)
                if structure is None or structure.empty:
                    continue
                
                # Get unique chains
                chains = structure['auth_chain_id'].unique()
                
                for chain in chains:
                    # Filter to chain
                    chain_data = structure[structure['auth_chain_id'] == chain]
                    
                    # Get CA atoms only for sequence
                    ca_atoms = chain_data[chain_data['atom_name'] == 'CA']
                    
                    if ca_atoms.empty:
                        continue
                    
                    # Sort by residue number
                    ca_atoms = ca_atoms.sort_values('auth_seq_id')
                    
                    # Convert 3-letter to 1-letter code
                    aa_map = {
                        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
                        'MSE': 'M', 'UNK': 'X'
                    }
                    
                    # Build sequence
                    sequence = []
                    for _, residue in ca_atoms.iterrows():
                        res_name = residue.get('res_name3l', residue.get('auth_comp_id', 'UNK')).upper()
                        aa = aa_map.get(res_name, 'X')
                        sequence.append(aa)
                    
                    seq_str = ''.join(sequence)
                    if len(seq_str) > 100:  # Only keep reasonable length sequences
                        seq_id = f"{pdb_id}_{chain}"
                        sequences[seq_id] = seq_str
                        logger.info(f"  {seq_id}: {len(seq_str)} residues")
                        
            except Exception as e:
                logger.warning(f"Failed to extract sequence from {pdb_id}: {e}")
        
        return sequences
    
    def create_grn_reference(self):
        """Create or load GRN reference table."""
        logger.info("\nPreparing GRN reference table...")
        
        # Try to load existing reference
        try:
            self.grn_proc.load_grn_table("ref/mo_ref")
            logger.info(f"Loaded reference table with {len(self.grn_proc.data)} proteins")
            return self.grn_proc.data
        except:
            # Create minimal reference for bacteriorhodopsin
            logger.info("Creating bacteriorhodopsin reference...")
            
            # Key positions for bacteriorhodopsin
            br_grns = {
                '1.46': 'W', '1.50': 'A', '1.53': 'T',
                '2.46': 'V', '2.50': 'Y', '2.53': 'A',
                '3.46': 'L', '3.50': 'D', '3.53': 'L',
                '4.46': 'L', '4.50': 'F', '4.53': 'V',
                '5.46': 'I', '5.50': 'G', '5.53': 'A',
                '6.46': 'G', '6.48': 'W', '6.50': 'P', '6.53': 'Y',
                '7.46': 'F', '7.50': 'K', '7.53': 'S',  # K216 is Schiff base
            }
            
            # Add surrounding positions
            for tm in range(1, 8):
                for offset in [-4, -3, -2, -1, 1, 2, 3, 4]:
                    pos = 50 + offset
                    grn = f"{tm}.{pos}"
                    if grn not in br_grns and 36 <= pos <= 62:
                        br_grns[grn] = 'X'
            
            df = pd.DataFrame([br_grns], index=['BR_reference'])
            df = df.fillna('-')
            
            # Sort columns
            sorted_cols = sort_grns_str(df.columns.tolist())
            df = df[sorted_cols]
            
            self.grn_proc.data = df
            self.grn_proc.ids = df.index.tolist()
            self.grn_proc.grns = df.columns.tolist()
            
            return df
    
    def assign_grns_to_sequences(self, sequences):
        """Assign GRN numbers to extracted sequences."""
        logger.info("\nAssigning GRNs to sequences...")
        
        # Get reference sequences
        ref_sequences = {}
        for protein_id in self.grn_proc.ids:
            seq = get_seq(protein_id, self.grn_proc.data)
            if seq and len(seq) > 50:
                ref_sequences[protein_id] = seq
        
        aligner = init_aligner()
        grn_assignments = {}
        
        for seq_id, query_seq in sequences.items():
            logger.info(f"\nProcessing {seq_id}...")
            
            # Find best reference
            best_ref = list(ref_sequences.keys())[0]
            best_score = -1000
            
            try:
                # Use MMseqs2 if available
                hits = mmseqs2_align2(
                    query_seqs={seq_id: query_seq},
                    ref_seqs=ref_sequences
                )
                if hits is not None and not hits.empty:
                    best_ref = hits.iloc[0]['target_id']
                    logger.info(f"  Best match: {best_ref}")
            except:
                # Fallback to pairwise alignment
                for ref_id, ref_seq in ref_sequences.items():
                    alignment = align_blosum62(query_seq, ref_seq, aligner)
                    if alignment.score > best_score:
                        best_score = alignment.score
                        best_ref = ref_id
            
            # Perform alignment with best reference
            ref_seq = ref_sequences[best_ref]
            alignment = align_blosum62(query_seq, ref_seq, aligner)
            formatted = format_alignment(alignment)
            
            # Transfer GRN annotations
            ref_row = self.grn_proc.data.loc[best_ref]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
            seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
            
            # Initialize GRN row from alignment
            new_row = init_row_from_alignment(formatted, seq_pos2grn)
            
            # Convert to position mapping
            grn_mapping = {}
            for grn, (res, pos) in new_row.items():
                if res != '-':
                    grn_mapping[pos] = grn
            
            grn_assignments[seq_id] = grn_mapping
            logger.info(f"  Assigned {len(grn_mapping)} GRN positions")
            
            # Show key positions
            key_positions = ['1.50', '2.50', '3.50', '7.50']
            for grn, (res, pos) in new_row.items():
                if grn in key_positions and res != '-':
                    logger.info(f"    {grn}: {res}{pos}")
        
        return grn_assignments
    
    def analyze_grn_structure_relationships(self, pdb_ids, grn_assignments):
        """Analyze relationships between GRN positions and structure."""
        logger.info("\nAnalyzing GRN-structure relationships...")
        
        results = []
        
        for pdb_id in pdb_ids:
            # Load structure
            structure = self.struct_proc.load_entity(pdb_id)
            if structure is None or structure.empty:
                continue
            
            # Get chains with GRN assignments
            for chain in structure['auth_chain_id'].unique():
                seq_id = f"{pdb_id}_{chain}"
                if seq_id not in grn_assignments:
                    continue
                
                grn_mapping = grn_assignments[seq_id]
                chain_data = structure[structure['auth_chain_id'] == chain]
                
                logger.info(f"\n{seq_id}:")
                
                # Analyze key GRN positions
                key_grns = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
                
                for grn in key_grns:
                    # Find residue with this GRN
                    res_pos = None
                    for pos, assigned_grn in grn_mapping.items():
                        if assigned_grn == grn:
                            res_pos = pos
                            break
                    
                    if res_pos is None:
                        continue
                    
                    # Get residue data
                    residue_data = chain_data[chain_data['auth_seq_id'] == res_pos]
                    
                    if not residue_data.empty:
                        res_name = residue_data.iloc[0].get('res_name3l', 'UNK')
                        
                        # Calculate properties
                        ca_atom = residue_data[residue_data['atom_name'] == 'CA']
                        if not ca_atom.empty:
                            # Get coordinates
                            x = ca_atom.iloc[0]['x_coord']
                            y = ca_atom.iloc[0]['y_coord']
                            z = ca_atom.iloc[0]['z_coord']
                            
                            # Store result
                            results.append({
                                'pdb_chain': seq_id,
                                'grn': grn,
                                'res_pos': res_pos,
                                'res_name': res_name,
                                'x': x,
                                'y': y,
                                'z': z
                            })
                            
                            logger.info(f"  {grn}: {res_name}{res_pos} at ({x:.1f}, {y:.1f}, {z:.1f})")
        
        return pd.DataFrame(results)
    
    def calculate_grn_distances(self, grn_structure_df):
        """Calculate distances between key GRN positions."""
        logger.info("\nCalculating distances between key GRN positions...")
        
        # Group by structure
        for pdb_chain, group in grn_structure_df.groupby('pdb_chain'):
            logger.info(f"\n{pdb_chain}:")
            
            # Calculate distances between conserved positions
            positions = group.set_index('grn')
            
            # Key pairs to analyze
            key_pairs = [
                ('3.50', '7.50'),  # D85-K216 in BR (proton transfer)
                ('1.50', '2.50'),  # TM1-TM2 interaction
                ('3.50', '6.50'),  # TM3-TM6 interaction
            ]
            
            for grn1, grn2 in key_pairs:
                if grn1 in positions.index and grn2 in positions.index:
                    pos1 = positions.loc[grn1]
                    pos2 = positions.loc[grn2]
                    
                    # Calculate distance
                    dist = np.sqrt(
                        (pos1['x'] - pos2['x'])**2 +
                        (pos1['y'] - pos2['y'])**2 +
                        (pos1['z'] - pos2['z'])**2
                    )
                    
                    logger.info(f"  {grn1} ({pos1['res_name']}{pos1['res_pos']}) - "
                               f"{grn2} ({pos2['res_name']}{pos2['res_pos']}): {dist:.1f} Å")
    
    def run_demonstration(self):
        """Run the complete GRN-structure integration demonstration."""
        logger.info("="*80)
        logger.info("GRN-STRUCTURE INTEGRATION DEMONSTRATION")
        logger.info("="*80)
        
        # 1. Download structures
        pdb_ids = self.download_structures()
        if not pdb_ids:
            logger.error("No structures available")
            return
        
        # 2. Extract sequences
        sequences = self.extract_sequences_from_structures(pdb_ids[:3])  # Use first 3
        if not sequences:
            logger.error("No sequences extracted")
            return
        
        # 3. Create/load GRN reference
        ref_table = self.create_grn_reference()
        
        # 4. Assign GRNs to sequences
        grn_assignments = self.assign_grns_to_sequences(sequences)
        
        # 5. Analyze GRN-structure relationships
        grn_structure_df = self.analyze_grn_structure_relationships(pdb_ids[:3], grn_assignments)
        
        # 6. Calculate distances
        if not grn_structure_df.empty:
            self.calculate_grn_distances(grn_structure_df)
        
        # 7. Save results
        if not grn_structure_df.empty:
            output_file = "grn_structure_analysis.csv"
            grn_structure_df.to_csv(output_file, index=False)
            logger.info(f"\nSaved analysis to: {output_file}")
        
        logger.info("\n" + "="*80)
        logger.info("DEMONSTRATION COMPLETE")
        logger.info("="*80)


def main():
    """Run the GRN-structure integration demonstration."""
    demo = GRNStructureIntegration()
    demo.run_demonstration()


if __name__ == "__main__":
    main()
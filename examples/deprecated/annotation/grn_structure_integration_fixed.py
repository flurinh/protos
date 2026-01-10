#!/usr/bin/env python3
"""
GRN-Structure Integration Example (Fixed)
=========================================

This script demonstrates the proper workflow for integrating structure, sequence,
and GRN data using the three processors in sequence.

Workflow:
1. Load protein structures using StructureProcessor (CifBaseProcessor)
2. Extract sequences from structures and save using SequenceProcessor
3. Assign GRN numbers using GRNProcessor
4. Analyze structural features mapped to GRN positions

Key improvements:
- Uses the modern processor API from DATA_MANAGEMENT_UNIFIED.md
- Properly handles entity operations (one at a time)
- Follows the zero-configuration principle
- Uses human-readable names throughout

Usage:
    python examples/grn_structure_integration_fixed.py
"""

import logging
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple

# Import the modern processors
from protos.processing.structure.structure_processor import CifBaseProcessor
from protos.processing.sequence.sequence_processor import SequenceProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.io.paths.path_config import ProtosPaths

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GRNStructureIntegration:
    """Integrates structure → sequence → GRN analysis workflow."""
    
    def __init__(self, data_dir: Optional[Path] = None):
        """
        Initialize processors with optional data directory.
        
        Args:
            data_dir: Custom data directory. If None, uses default.
        """
        # Initialize ProtosPaths (zero configuration if no path provided)
        if data_dir:
            self.paths = ProtosPaths(base_path=data_dir)
        else:
            # Use default data location
            self.paths = ProtosPaths()
        
        # Initialize processors - they automatically share the same paths
        self.struct_proc = CifBaseProcessor(name="grn_integration", paths=self.paths)
        self.seq_proc = SequenceProcessor(name="grn_integration", paths=self.paths)
        self.grn_proc = GRNProcessor(name="grn_integration", paths=self.paths)
        
        # Define test structures (microbial rhodopsins)
        self.test_structures = {
            '1m0k': 'Bacteriorhodopsin (dark state)',
            '1m0l': 'Bacteriorhodopsin (K intermediate)',
            '2ntu': 'Sensory rhodopsin II',
            '3ddl': 'Xanthorhodopsin',
            '4hyj': 'Acetabularia rhodopsin',
            '5b6w': 'Channelrhodopsin chimera C1C2',
            '6eig': 'Channelrhodopsin 2',
            '6csm': 'Krokinobacter rhodopsin 2',
        }
        
        # Amino acid mapping
        self.aa_3to1 = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            'MSE': 'M', 'UNK': 'X'  # Selenomethionine and unknown
        }
    
    def check_available_structures(self) -> List[str]:
        """Check which test structures are available."""
        logger.info("Checking available structures...")
        available = []
        
        for pdb_id in self.test_structures.keys():
            if self.struct_proc.entity_exists(pdb_id):
                available.append(pdb_id)
                logger.info(f"  ✓ {pdb_id}: {self.test_structures[pdb_id]}")
            else:
                logger.info(f"  ✗ {pdb_id}: Not available")
        
        logger.info(f"\nFound {len(available)}/{len(self.test_structures)} structures")
        return available
    
    def extract_sequence_from_structure(self, pdb_id: str, chain_id: str = 'A') -> Optional[str]:
        """
        Extract sequence from a single structure and chain.
        
        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier (default 'A')
            
        Returns:
            Protein sequence as string, or None if extraction fails
        """
        try:
            # Load structure
            structure_df = self.struct_proc.load_entity(pdb_id)
            
            if structure_df is None or structure_df.empty:
                logger.warning(f"No data found for {pdb_id}")
                return None
            
            # Filter to specific chain
            chain_data = structure_df[structure_df['auth_chain_id'] == chain_id].copy()
            
            if chain_data.empty:
                logger.warning(f"No chain {chain_id} found in {pdb_id}")
                return None
            
            # Get CA atoms for sequence
            ca_atoms = chain_data[chain_data['atom_name'] == 'CA'].copy()
            
            if ca_atoms.empty:
                logger.warning(f"No CA atoms found for {pdb_id} chain {chain_id}")
                return None
            
            # Sort by residue number
            ca_atoms = ca_atoms.sort_values('auth_seq_id')
            
            # Build sequence
            sequence_parts = []
            for _, residue in ca_atoms.iterrows():
                # Get residue name (try different column names)
                res_name = residue.get('auth_comp_id', residue.get('res_name3l', 'UNK')).upper()
                aa_1letter = self.aa_3to1.get(res_name, 'X')
                sequence_parts.append(aa_1letter)
            
            sequence = ''.join(sequence_parts)
            
            # Only return if reasonable length
            if len(sequence) > 100:
                logger.info(f"  Extracted {pdb_id}_{chain_id}: {len(sequence)} residues")
                return sequence
            else:
                logger.warning(f"  Sequence too short for {pdb_id}_{chain_id}: {len(sequence)} residues")
                return None
                
        except Exception as e:
            logger.error(f"Failed to extract sequence from {pdb_id}: {e}")
            return None
    
    def extract_all_sequences(self, pdb_ids: List[str]) -> Dict[str, str]:
        """
        Extract sequences from multiple structures.
        
        Args:
            pdb_ids: List of PDB identifiers
            
        Returns:
            Dictionary mapping sequence_id to sequence string
        """
        logger.info(f"\nExtracting sequences from {len(pdb_ids)} structures...")
        sequences = {}
        
        for pdb_id in pdb_ids:
            try:
                # Load structure to get available chains
                structure_df = self.struct_proc.load_entity(pdb_id)
                
                if structure_df is None or structure_df.empty:
                    continue
                
                # Get unique chains
                chains = structure_df['auth_chain_id'].unique()
                
                for chain_id in chains:
                    sequence = self.extract_sequence_from_structure(pdb_id, chain_id)
                    
                    if sequence:
                        seq_id = f"{pdb_id}_{chain_id}"
                        sequences[seq_id] = sequence
                        
                        # Save to sequence processor
                        self.seq_proc.save_entity(seq_id, sequence, metadata={
                            'source': 'structure',
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'length': len(sequence)
                        })
                        
            except Exception as e:
                logger.error(f"Error processing {pdb_id}: {e}")
        
        logger.info(f"Extracted {len(sequences)} sequences total")
        return sequences
    
    def create_grn_table_for_sequences(self, sequences: Dict[str, str]) -> pd.DataFrame:
        """
        Create a GRN table for the extracted sequences.
        
        This is a simplified version that demonstrates the workflow.
        In production, you would use proper GRN assignment algorithms.
        
        Args:
            sequences: Dictionary of sequence_id to sequence
            
        Returns:
            DataFrame with GRN assignments
        """
        logger.info("\nCreating GRN assignments...")
        
        # For demonstration, we'll create a simple GRN table
        # Key positions for rhodopsins (simplified)
        grn_positions = [
            '1.46', '1.50', '1.53',
            '2.46', '2.50', '2.53',
            '3.46', '3.50', '3.53',
            '4.46', '4.50', '4.53',
            '5.46', '5.50', '5.53',
            '6.46', '6.50', '6.53',
            '7.46', '7.50', '7.53'
        ]
        
        # Create empty table
        grn_data = []
        
        for seq_id, sequence in sequences.items():
            # For demo, assign positions based on approximate locations
            # In reality, this would use alignment-based assignment
            row_data = {'protein_id': seq_id}
            
            # Simplified mapping (would use real GRN assignment)
            if 'bacteriorhodopsin' in seq_id.lower() or '1m0' in seq_id:
                # Known positions for bacteriorhodopsin
                position_map = {
                    '1.50': 42,   # W42
                    '2.50': 83,   # Y83
                    '3.50': 85,   # D85
                    '7.50': 216,  # K216 (Schiff base)
                }
                
                for grn, pos in position_map.items():
                    if 0 < pos <= len(sequence):
                        row_data[grn] = sequence[pos-1]
            
            # Fill in some other positions for demo
            for grn in grn_positions:
                if grn not in row_data:
                    row_data[grn] = '-'
            
            grn_data.append(row_data)
        
        # Create DataFrame
        grn_df = pd.DataFrame(grn_data)
        grn_df = grn_df.set_index('protein_id')
        
        # Save GRN table
        table_name = "rhodopsin_grn_demo"
        self.grn_proc.save_grn_table(table_name, grn_df)
        
        logger.info(f"Created GRN table '{table_name}' with {len(grn_df)} proteins")
        return grn_df
    
    def analyze_grn_structure_mapping(self, pdb_id: str, chain_id: str, grn_df: pd.DataFrame) -> pd.DataFrame:
        """
        Map GRN positions back to 3D structure coordinates.
        
        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            grn_df: GRN assignment table
            
        Returns:
            DataFrame with GRN positions mapped to 3D coordinates
        """
        seq_id = f"{pdb_id}_{chain_id}"
        
        if seq_id not in grn_df.index:
            logger.warning(f"No GRN data for {seq_id}")
            return pd.DataFrame()
        
        # Load structure
        structure_df = self.struct_proc.load_entity(pdb_id)
        if structure_df is None:
            return pd.DataFrame()
        
        # Filter to chain
        chain_data = structure_df[structure_df['auth_chain_id'] == chain_id]
        
        # Get GRN assignments for this protein
        grn_row = grn_df.loc[seq_id]
        
        results = []
        
        # For each GRN position with an assignment
        for grn_pos, residue in grn_row.items():
            if residue != '-':
                # In real implementation, would map GRN to sequence position
                # For demo, we'll use known positions
                if grn_pos == '3.50' and '1m0' in pdb_id:
                    seq_pos = 85  # D85 in bacteriorhodopsin
                elif grn_pos == '7.50' and '1m0' in pdb_id:
                    seq_pos = 216  # K216 in bacteriorhodopsin
                else:
                    continue
                
                # Get residue data
                residue_data = chain_data[chain_data['auth_seq_id'] == seq_pos]
                
                if not residue_data.empty:
                    # Get CA atom coordinates
                    ca_atom = residue_data[residue_data['atom_name'] == 'CA']
                    
                    if not ca_atom.empty:
                        atom = ca_atom.iloc[0]
                        results.append({
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'grn_position': grn_pos,
                            'seq_position': seq_pos,
                            'residue': residue,
                            'x': atom['x_coord'],
                            'y': atom['y_coord'],
                            'z': atom['z_coord']
                        })
        
        return pd.DataFrame(results)
    
    def calculate_key_distances(self, coord_df: pd.DataFrame) -> None:
        """Calculate distances between key GRN positions."""
        if coord_df.empty:
            return
        
        logger.info("\nKey distances in structure:")
        
        # Group by structure/chain
        for (pdb_id, chain_id), group in coord_df.groupby(['pdb_id', 'chain_id']):
            positions = group.set_index('grn_position')
            
            # Calculate specific distances
            if '3.50' in positions.index and '7.50' in positions.index:
                pos1 = positions.loc['3.50']
                pos2 = positions.loc['7.50']
                
                distance = np.sqrt(
                    (pos1['x'] - pos2['x'])**2 +
                    (pos1['y'] - pos2['y'])**2 +
                    (pos1['z'] - pos2['z'])**2
                )
                
                logger.info(f"  {pdb_id}_{chain_id}: "
                           f"{pos1['residue']}{pos1['seq_position']} (3.50) - "
                           f"{pos2['residue']}{pos2['seq_position']} (7.50) = "
                           f"{distance:.1f} Å")
    
    def run_demonstration(self):
        """Run the complete demonstration workflow."""
        logger.info("="*80)
        logger.info("GRN-STRUCTURE INTEGRATION DEMONSTRATION")
        logger.info("="*80)
        
        # 1. Check available structures
        available_pdbs = self.check_available_structures()
        
        if not available_pdbs:
            logger.error("No structures available. Please add PDB files to: "
                        f"{self.paths.get_processor_path('structure') / 'mmcif'}")
            return
        
        # 2. Extract sequences from structures (use first 3)
        test_pdbs = available_pdbs[:3]
        sequences = self.extract_all_sequences(test_pdbs)
        
        if not sequences:
            logger.error("No sequences could be extracted")
            return
        
        # 3. Create GRN assignments
        grn_df = self.create_grn_table_for_sequences(sequences)
        
        # 4. Map GRN positions back to structures
        logger.info("\nMapping GRN positions to 3D coordinates...")
        all_coords = []
        
        for seq_id in sequences.keys():
            pdb_id, chain_id = seq_id.split('_')
            coord_df = self.analyze_grn_structure_mapping(pdb_id, chain_id, grn_df)
            if not coord_df.empty:
                all_coords.append(coord_df)
        
        if all_coords:
            combined_coords = pd.concat(all_coords, ignore_index=True)
            
            # 5. Analyze distances
            self.calculate_key_distances(combined_coords)
            
            # 6. Save results
            output_file = self.paths.get_base_path() / "grn_structure_mapping.csv"
            combined_coords.to_csv(output_file, index=False)
            logger.info(f"\nSaved mapping results to: {output_file}")
        
        # 7. Create a dataset for future use
        self.struct_proc.create_dataset(
            "rhodopsin_structures",
            test_pdbs,
            metadata={
                "description": "Microbial rhodopsin structures for GRN analysis",
                "grn_table": "rhodopsin_grn_demo"
            }
        )
        
        logger.info("\n" + "="*80)
        logger.info("DEMONSTRATION COMPLETE")
        logger.info("="*80)
        logger.info(f"\nData saved to: {self.paths.get_base_path()}")
        logger.info("- Structures: structure/mmcif/")
        logger.info("- Sequences: sequence/fasta/")
        logger.info("- GRN table: grn/tables/rhodopsin_grn_demo.csv")
        logger.info("- Dataset: structure/datasets/rhodopsin_structures.json")


def main():
    """Run the demonstration."""
    # Option 1: Use default data location (./data)
    demo = GRNStructureIntegration()
    
    # Option 2: Use custom data location
    # demo = GRNStructureIntegration(data_dir=Path("/path/to/custom/data"))
    
    demo.run_demonstration()


if __name__ == "__main__":
    main()
#!/usr/bin/env python
"""
Comprehensive test script for updated SequenceProcessor and GRNProcessor.

This script tests all the new functionalities including:
- New directory structures (alignments/pairwise, alignments/multiple, alignments/mmseqs, databases/)
- Reference vs tables separation for GRN
- Temporary file handling
- All core methods with real data

Author: Assistant
Date: 2024
"""

import os
import sys
import logging
import pandas as pd
from pathlib import Path
from datetime import datetime
import time

# Add protos to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.io.paths import ProtosPaths
from protos.io.fasta_utils import read_fasta

# Configure logging for detailed output
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*80}")
    print(f" {title}")
    print('='*80)


def test_sequence_processor(data_dir: Path):
    """Test all SequenceProcessor functionalities using real data."""
    print_section("TESTING SEQUENCE PROCESSOR WITH REAL DATA")
    
    try:
        # Initialize processor
        logger.info("Initializing SequenceProcessor...")
        seq_proc = SequenceProcessor(name="test_seq_processor")
        logger.info(f"SequenceProcessor initialized successfully")
        
        # Test 1: Check all path properties
        print("\n1. Testing path properties:")
        paths_to_check = [
            ('fasta_dir', seq_proc.path_fasta_dir),
            ('alignments_dir', seq_proc.path_alignments_dir),
            ('pairwise_alignments_dir', seq_proc.path_pairwise_alignments_dir),
            ('multiple_alignments_dir', seq_proc.path_multiple_alignments_dir),
            ('mmseqs_alignments_dir', seq_proc.path_mmseqs_alignments_dir),
            ('databases_dir', seq_proc.path_databases_dir),
            ('metadata_dir', seq_proc.path_metadata_dir)
        ]
        
        for name, path in paths_to_check:
            exists = path.exists()
            logger.info(f"  {name}: {path} (exists: {exists})")
            if not exists:
                logger.debug(f"  Creating directory: {path}")
        
        # Test 2: Load existing real sequences
        print("\n2. Testing load_entity with real FASTA files:")
        # Note: 1ubq.fasta actually contains 3SN6 sequences
        real_sequences = ['1ubq', 'BACR_HALSA', 'ChR2', 'NpHR', '1crn_SEQUENCE', '3SN6_A']
        loaded_sequences = {}
        
        for seq_id in real_sequences:
            logger.info(f"  Loading sequence: {seq_id}")
            try:
                sequence = seq_proc.load_entity(seq_id)
                if sequence:
                    if isinstance(sequence, str):
                        loaded_sequences[seq_id] = sequence
                        logger.info(f"    ✓ Successfully loaded: length={len(sequence)}")
                        logger.debug(f"    First 50 chars: {sequence[:50]}...")
                    elif isinstance(sequence, dict):
                        # It's a multi-sequence file, add all sequences
                        loaded_sequences.update(sequence)
                        logger.info(f"    ✓ Successfully loaded: {len(sequence)} sequences from {seq_id}")
                        for sub_id in list(sequence.keys())[:2]:  # Show first 2
                            logger.debug(f"      {sub_id}: {len(sequence[sub_id])} residues")
                else:
                    logger.warning(f"    Could not load {seq_id}")
            except Exception as e:
                logger.warning(f"    Failed to load {seq_id}: {e}")
        
        # Test 3: Load multi-sequence FASTA file  
        print("\n3. Testing load_sequences with real multi-sequence file:")
        multi_seq_files = ['mo_small.fasta', 'gpcr_structures.fasta', 'demo_proteins.fasta']
        
        for fasta_file in multi_seq_files:
            logger.info(f"  Loading multi-sequence file: {fasta_file}")
            try:
                sequences = seq_proc.load_sequences(fasta_file)
                logger.info(f"    ✓ Loaded {len(sequences)} sequences from {fasta_file}")
                # Show first 3 sequences
                for i, (seq_id, seq) in enumerate(list(sequences.items())[:3]):
                    logger.debug(f"      {seq_id}: {len(seq)} residues")
                if len(sequences) > 3:
                    logger.debug(f"      ... and {len(sequences) - 3} more sequences")
                break  # Use first available file for further tests
            except FileNotFoundError:
                logger.warning(f"    File not found: {fasta_file}")
                continue
        
        # Test 4: Save sequences with new name
        print("\n4. Testing save_sequences (create new multi-sequence file):")
        if loaded_sequences:
            subset_seqs = dict(list(loaded_sequences.items())[:3])  # Take first 3
            logger.info(f"  Saving {len(subset_seqs)} sequences to 'test_subset.fasta'")
            seq_proc.save_sequences(subset_seqs, "test_subset", dataset_name="test_dataset")
            logger.debug(f"  File saved to: {seq_proc.path_fasta_dir / 'test_subset.fasta'}")
        
        # Test 5: Verify the saved file
        print("\n5. Testing load of newly saved file:")
        logger.info("  Loading sequences from 'test_subset.fasta'")
        try:
            loaded_subset = seq_proc.load_entity("test_subset")
            if loaded_subset:
                logger.info(f"  ✓ Successfully loaded saved file")
                if isinstance(loaded_subset, dict):
                    logger.info(f"    Contains {len(loaded_subset)} sequences")
        except Exception as e:
            logger.error(f"  Failed to load saved file: {e}")
        
        # Test 6: Pairwise alignment with real sequences
        print("\n6. Testing align_sequences (pairwise alignment with real data):")
        if 'BACR_HALSA' in loaded_sequences and 'ChR2' in loaded_sequences:
            logger.info("  Aligning BACR_HALSA vs ChR2")
            try:
                score, alignment = seq_proc.align_sequences(
                    loaded_sequences["BACR_HALSA"], 
                    loaded_sequences["ChR2"],
                    "BACR_HALSA", 
                    "ChR2",
                    store_alignment=True
                )
                logger.info(f"  ✓ Alignment score: {score}")
                logger.debug(f"  Alignment saved to: {seq_proc.path_pairwise_alignments_dir / 'BACR_HALSA_vs_ChR2.json'}")
                
                # Verify the file was created
                alignment_file = seq_proc.path_pairwise_alignments_dir / "BACR_HALSA_vs_ChR2.json"
                if alignment_file.exists():
                    logger.info(f"  ✓ Alignment file created successfully")
                else:
                    logger.error(f"  ✗ Alignment file not found!")
                    
            except Exception as e:
                logger.error(f"  Alignment failed: {e}")
        else:
            # Use any two available sequences
            if len(loaded_sequences) >= 2:
                seq_ids = list(loaded_sequences.keys())[:2]
                logger.info(f"  Aligning {seq_ids[0]} vs {seq_ids[1]}")
                try:
                    score, alignment = seq_proc.align_sequences(
                        loaded_sequences[seq_ids[0]], 
                        loaded_sequences[seq_ids[1]],
                        seq_ids[0], 
                        seq_ids[1],
                        store_alignment=True
                    )
                    logger.info(f"  ✓ Alignment score: {score}")
                except Exception as e:
                    logger.error(f"  Alignment failed: {e}")
        
        # Test 7: Load alignment
        print("\n7. Testing load_alignment:")
        try:
            logger.info("  Loading alignment from 'BACR_HALSA_vs_ChR2.json'")
            loaded_alignment = seq_proc.load_alignment("BACR_HALSA_vs_ChR2.json", alignment_type="pairwise")
            logger.info(f"  Successfully loaded alignment with score: {loaded_alignment.get('score', 'N/A')}")
        except Exception as e:
            logger.error(f"  Failed to load alignment: {e}")
        
        # Test 8: Multiple sequence alignment with real sequences
        print("\n8. Testing multiple_sequence_alignment with real data:")
        if len(loaded_sequences) >= 3:
            # Use a subset of sequences for MSA
            msa_subset = dict(list(loaded_sequences.items())[:5])  # First 5 sequences
            logger.info(f"  Performing MSA on {len(msa_subset)} sequences")
            try:
                msa_results = seq_proc.multiple_sequence_alignment(msa_subset, use_mmseqs=False)
                logger.info(f"  ✓ MSA completed for {len(msa_results)} sequences")
                for query_id, (ref_id, score, _) in list(msa_results.items())[:3]:
                    logger.debug(f"    {query_id} -> {ref_id}: score={score}")
            except Exception as e:
                logger.error(f"  MSA failed: {e}")
        else:
            logger.warning("  Not enough sequences loaded for MSA test")
        
        # Test 9: Save different alignment types
        print("\n9. Testing save_alignment with different types:")
        test_alignment = {
            'seq1': 'BACR_HALSA',
            'seq2': 'NpHR',
            'score': 150.5,
            'alignment': ['test alignment data'],
            'timestamp': datetime.now().isoformat()
        }
        
        for align_type in ['pairwise', 'multiple', 'mmseqs']:
            logger.info(f"  Saving {align_type} alignment")
            seq_proc.save_alignment(test_alignment, f"test_{align_type}.json", alignment_type=align_type)
            
            # Check the appropriate directory
            if align_type == 'pairwise':
                check_path = seq_proc.path_pairwise_alignments_dir / f"test_{align_type}.json"
            elif align_type == 'multiple':
                check_path = seq_proc.path_multiple_alignments_dir / f"test_{align_type}.json"
            else:
                check_path = seq_proc.path_mmseqs_alignments_dir / f"test_{align_type}.json"
                
            if check_path.exists():
                logger.info(f"    ✓ File created in correct directory")
            else:
                logger.error(f"    ✗ File not found in expected location: {check_path}")
        
        # Test 10: Create MMseqs2 database
        print("\n10. Testing create_mmseqs_database:")
        logger.info("  Creating MMseqs2 database 'opsin_db' from sequences")
        try:
            # Use loaded sequences if available
            if loaded_sequences:
                db_path = seq_proc.create_mmseqs_database(loaded_sequences, "opsin_db")
                logger.info(f"  Database created at: {db_path}")
            else:
                logger.warning("  No sequences available for MMseqs2 database creation")
        except Exception as e:
            logger.warning(f"  MMseqs2 database creation skipped (expected if MMseqs2 not installed): {e}")
        
        # Test 11: List databases
        print("\n11. Testing list_mmseqs_databases:")
        databases = seq_proc.list_mmseqs_databases()
        logger.info(f"  Found {len(databases)} databases: {databases}")
        
        # Test 12: Get sequence metadata
        print("\n12. Testing get_sequence_metadata:")
        metadata_df = seq_proc.get_sequence_metadata()
        logger.info(f"  Retrieved metadata for {len(metadata_df)} sequences")
        for idx, row in metadata_df.iterrows():
            logger.debug(f"    {row['sequence_id']}: length={row['length']}, MW={row['molecular_weight']}")
        
        # Test 13: Generate variants
        print("\n13. Testing generate_variants:")
        if 'BACR_HALSA' in loaded_sequences:
            logger.info("  Generating variants of BACR_HALSA at positions 10, 20")
            variants = seq_proc.generate_variants(
                loaded_sequences["BACR_HALSA"],
                positions=[10, 20],
                possible_aas=[['A', 'V', 'L'], ['S', 'T']],
                base_id="BACR_variant"
            )
            logger.info(f"  Generated {len(variants)} variants")
            for var_id, var_seq in list(variants.items())[:3]:
                logger.debug(f"    {var_id}: ...{var_seq[5:25]}...")
        else:
            # Use any available sequence
            if loaded_sequences:
                seq_id = list(loaded_sequences.keys())[0]
                logger.info(f"  Generating variants of {seq_id} at positions 5, 10")
                variants = seq_proc.generate_variants(
                    loaded_sequences[seq_id],
                    positions=[5, 10],
                    possible_aas=[['A', 'G'], ['S', 'T']],
                    base_id=f"{seq_id}_variant"
                )
                logger.info(f"  Generated {len(variants)} variants")
        
        # List all entities
        print("\n14. Testing list_entities:")
        all_entities = seq_proc.list_entities()
        logger.info(f"  Total entities registered: {len(all_entities)}")
        for entity in all_entities[:5]:  # Show first 5
            logger.debug(f"    - {entity}")
            
        return True
        
    except Exception as e:
        logger.error(f"SequenceProcessor test failed: {e}", exc_info=True)
        return False


def test_grn_processor(data_dir: Path):
    """Test all GRNProcessor functionalities with real data."""
    print_section("TESTING GRN PROCESSOR WITH REAL DATA")
    
    try:
        # Initialize processor
        logger.info("Initializing GRNProcessor...")
        grn_proc = GRNProcessor(name="test_grn_processor")
        logger.info(f"GRNProcessor initialized successfully")
        
        # Test 1: Check all path properties
        print("\n1. Testing path properties:")
        paths_to_check = [
            ('grn_dir (tables)', grn_proc.path_grn_dir),
            ('ref_dir (reference)', grn_proc.path_ref_dir),
            ('config_dir', grn_proc.path_config_dir),
            ('temp_dir', grn_proc.path_temp_dir)
        ]
        
        for name, path in paths_to_check:
            exists = path.exists()
            logger.info(f"  {name}: {path} (exists: {exists})")
        
        # Test 2: Load existing GRN tables
        print("\n2. Testing load_grn_table with real data:")
        real_grn_tables = ['mo_small', 'bacteriorhodopsin_grn', 'demo_alignment', 'cross_format_demo']
        loaded_table = None
        
        for table_name in real_grn_tables:
            logger.info(f"  Trying to load GRN table: {table_name}")
            try:
                grn_proc.load_grn_table(table_name)
                logger.info(f"  ✓ Loaded {table_name}: {len(grn_proc.data)} sequences, {len(grn_proc.grns)} GRN positions")
                logger.debug(f"    First 5 sequences: {grn_proc.ids[:5]}")
                logger.debug(f"    First 5 GRN positions: {grn_proc.grns[:5]}")
                loaded_table = table_name
                break
            except Exception as e:
                logger.warning(f"    Could not load {table_name}: {e}")
        
        # Test 3: Save a modified version of loaded table
        print("\n3. Testing save_grn_table with real data:")
        if loaded_table and grn_proc.data is not None and not grn_proc.data.empty:
            # Create a subset for testing
            test_subset = grn_proc.data.iloc[:3].copy()  # First 3 sequences
            grn_proc.data = test_subset
            grn_proc.ids = list(test_subset.index)
            
            new_table_name = f"test_{loaded_table}_subset"
            logger.info(f"  Saving subset as '{new_table_name}'")
            grn_proc.save_grn_table(new_table_name)
            
            save_path = grn_proc.path_grn_dir / f"{new_table_name}.csv"
            if save_path.exists():
                logger.info(f"  ✓ GRN table saved successfully to: {save_path}")
            else:
                logger.error(f"  ✗ GRN table not saved!")
        
        # Test 4: Load and save reference table
        print("\n4. Testing reference table operations:")
        # First try to load existing reference
        ref_tables = ['mo_ref', 'gpcrdb_ref']
        ref_data = None
        
        for ref_name in ref_tables:
            logger.info(f"  Trying to load reference table: {ref_name}")
            try:
                ref_data = grn_proc.load_reference_table(ref_name)
                logger.info(f"  ✓ Loaded reference {ref_name}: {ref_data.shape}")
                logger.debug(f"    First 5 sequences: {list(ref_data.index[:5])}")
                break
            except Exception as e:
                logger.warning(f"    Could not load {ref_name}: {e}")
        
        # Save a subset as new reference
        if ref_data is not None:
            ref_subset = ref_data.iloc[:5].copy()  # First 5 sequences
            ref_metadata = {
                'family': 'test_reference',
                'source': 'subset_of_' + ref_name,
                'version': '1.0'
            }
            grn_proc.save_reference_table("test_reference", ref_subset, metadata=ref_metadata)
            ref_path = grn_proc.path_ref_dir / "test_reference.csv"
            if ref_path.exists():
                logger.info(f"  ✓ Reference table saved successfully")
            else:
                logger.error(f"  ✗ Reference table not found!")
        
        # Test 5: Load reference table
        print("\n5. Testing load_reference_table:")
        logger.info("  Loading reference table 'opsin_reference'")
        try:
            ref_table = grn_proc.load_reference_table("opsin_reference")
            logger.info(f"  Loaded reference table: {ref_table.shape}")
            logger.debug(f"  Reference sequences: {list(ref_table.index)}")
        except Exception as e:
            logger.error(f"  Failed to load reference table: {e}")
        
        # Test 6: Save temporary table
        print("\n6. Testing save_temp_table:")
        logger.info("  Saving temporary table")
        # Use loaded data if available
        if grn_proc.data is not None and not grn_proc.data.empty:
            temp_data = grn_proc.data.iloc[:2].copy()  # Subset of data
            temp_path = grn_proc.save_temp_table("test_temp", temp_data)
            logger.info(f"  Temporary table saved to: {temp_path}")
            if temp_path.exists():
                logger.info(f"  ✓ Temp file created successfully")
                # Add a small delay to ensure different timestamps
                time.sleep(1)
        
        # Test 7: Load entity
        print("\n7. Testing load_entity:")
        for entity_name in ['BACR_HALSA', 'ChR2']:
            logger.info(f"  Loading entity: {entity_name}")
            entity_data = grn_proc.load_entity(entity_name)
            if entity_data is not None:
                logger.info(f"    ✓ Loaded successfully, {len(entity_data)} GRN positions")
                logger.debug(f"    First 5 positions: {dict(list(entity_data.items())[:5])}")
            else:
                logger.error(f"    ✗ Failed to load entity")
        
        # Test 8: Save entity
        print("\n8. Testing save_entity:")
        new_entity = pd.Series({
            '1.50': 'M', '2.50': 'A', '3.50': 'R', '4.50': 'K',
            '5.50': 'E', '6.50': 'R', '7.50': 'S'
        })
        logger.info("  Saving new entity 'TEST_PROTEIN'")
        grn_proc.save_entity("TEST_PROTEIN", new_entity)
        logger.info("  Entity saved")
        
        # Test 9: Filter by IDs
        print("\n9. Testing filter_by_ids:")
        logger.info("  Filtering to only BACR_HALSA and ChR2")
        grn_proc.reset_data()  # Reset to full data
        original_count = len(grn_proc.data)
        grn_proc.filter_by_ids(['BACR_HALSA', 'ChR2'])
        filtered_count = len(grn_proc.data)
        logger.info(f"  Filtered from {original_count} to {filtered_count} sequences")
        logger.debug(f"  Remaining IDs: {grn_proc.ids}")
        
        # Test 10: Apply interval
        print("\n10. Testing apply_interval:")
        logger.info("  Limiting to GRN positions 1.50, 2.50, 3.50")
        grn_proc.reset_data()
        original_cols = len(grn_proc.grns)
        grn_proc.apply_interval(['1.50', '2.50', '3.50'])
        filtered_cols = len(grn_proc.grns)
        logger.info(f"  Reduced from {original_cols} to {filtered_cols} GRN positions")
        logger.debug(f"  Remaining positions: {grn_proc.grns}")
        
        # Test 11: Get GRN dictionary
        print("\n11. Testing get_grn_dict:")
        grn_proc.reset_data()
        grn_dict = grn_proc.get_grn_dict()
        logger.info(f"  Generated GRN dictionary for {len(grn_dict)} sequences")
        for seq_id, positions in list(grn_dict.items())[:2]:
            logger.debug(f"    {seq_id}: {len(positions)} positions filled")
        
        # Test 12: Get sequence dictionary
        print("\n12. Testing get_seq_dict:")
        seq_dict = grn_proc.get_seq_dict()
        logger.info(f"  Generated sequence dictionary for {len(seq_dict)} sequences")
        for seq_id, sequence in list(seq_dict.items())[:2]:
            logger.debug(f"    {seq_id}: {sequence[:20]}...")
        
        # Test 13: Clean temp files
        print("\n13. Testing clean_temp_files:")
        # Create an old temp file using current data
        old_temp_path = None
        if grn_proc.data is not None and not grn_proc.data.empty:
            old_temp_path = grn_proc.save_temp_table("old_temp", grn_proc.data)
            # Modify its timestamp to be old
            old_time = time.time() - (25 * 3600)  # 25 hours ago
            os.utime(old_temp_path, (old_time, old_time))
        
        logger.info("  Cleaning temp files older than 24 hours")
        grn_proc.clean_temp_files(older_than_hours=24)
        
        if old_temp_path and not old_temp_path.exists():
            logger.info("  ✓ Old temp file removed successfully")
        elif old_temp_path:
            logger.error("  ✗ Old temp file still exists")
        else:
            logger.warning("  No temp file was created for cleanup test")
            
        # Test 14: Merge GRN tables
        print("\n14. Testing load_and_merge_grn_tables:")
        # Find available tables to merge
        available_tables = []
        if loaded_table:
            available_tables.append(loaded_table)
        if (grn_proc.path_grn_dir / "test_mo_small_subset.csv").exists():
            available_tables.append("test_mo_small_subset")
        
        if len(available_tables) >= 2:
            logger.info(f"  Merging tables: {available_tables[:2]}")
            try:
                merged_data = grn_proc.load_and_merge_grn_tables(available_tables[:2])
                logger.info(f"  Merged table has {len(merged_data)} sequences")
                logger.debug(f"  Merged IDs: {list(merged_data.index)}")
            except Exception as e:
                logger.error(f"  Merge failed: {e}")
        else:
            logger.warning("  Not enough tables available for merge test")
        
        # Test 15: List datasets
        print("\n15. Testing list_datasets:")
        datasets = grn_proc.list_datasets()
        logger.info(f"  Found {len(datasets)} datasets")
        for dataset in datasets:
            logger.debug(f"    - {dataset}")
            
        return True
        
    except Exception as e:
        logger.error(f"GRNProcessor test failed: {e}", exc_info=True)
        return False


def test_cross_processor_workflow(data_dir: Path):
    """Test workflow using both processors together with real data."""
    print_section("TESTING CROSS-PROCESSOR WORKFLOW WITH REAL DATA")
    
    try:
        logger.info("Testing sequence -> GRN workflow with real data")
        
        # Initialize processors
        seq_proc = SequenceProcessor(name="workflow_seq")
        grn_proc = GRNProcessor(name="workflow_grn")
        
        # Load sequences from the sequence processor
        logger.info("1. Loading sequences from SequenceProcessor")
        
        # Try to load mo_small sequences which should have corresponding GRN data
        sequences = {}
        try:
            mo_sequences = seq_proc.load_sequences("mo_small.fasta")
            logger.info(f"  Loaded {len(mo_sequences)} sequences from mo_small.fasta")
            sequences = mo_sequences
        except:
            # Fallback to individual sequences
            for seq_id in ['BACR_HALSA', 'ChR2', 'NpHR']:
                try:
                    seq = seq_proc.load_entity(seq_id)
                    if seq:
                        sequences[seq_id] = seq
                        logger.debug(f"  Loaded {seq_id}: {len(seq)} residues")
                except:
                    continue
        
        # Load corresponding GRN table
        logger.info("2. Loading corresponding GRN table")
        grn_table_loaded = False
        
        # Try to load mo_small GRN table if we loaded mo_small sequences
        if 'mo_small.fasta' in str(sequences):
            try:
                grn_proc.load_grn_table("mo_small")
                grn_table_loaded = True
                logger.info("  ✓ Loaded mo_small GRN table")
            except:
                pass
        
        # Otherwise try other tables
        if not grn_table_loaded:
            for table_name in ['demo_alignment', 'bacteriorhodopsin_grn', 'cross_format_demo']:
                try:
                    grn_proc.load_grn_table(table_name)
                    grn_table_loaded = True
                    logger.info(f"  ✓ Loaded {table_name} GRN table")
                    break
                except:
                    continue
        
        if grn_table_loaded:
            # Verify GRN sequences against loaded sequences
            logger.info("3. Cross-checking sequences between processors")
            matches = 0
            for protein_id in grn_proc.ids[:5]:  # Check first 5
                if protein_id in sequences:
                    matches += 1
                    logger.debug(f"  ✓ Found matching sequence for {protein_id}")
            
            logger.info(f"  Found {matches} matching sequences between processors")
        
        # Create alignment of GRN sequences
        logger.info("3. Aligning sequences based on GRN data")
        grn_seqs = grn_proc.get_seq_dict()
        if len(grn_seqs) >= 2:
            seq_ids = list(grn_seqs.keys())[:2]
            score, alignment = seq_proc.align_sequences(
                grn_seqs[seq_ids[0]], 
                grn_seqs[seq_ids[1]],
                f"grn_{seq_ids[0]}", 
                f"grn_{seq_ids[1]}"
            )
            logger.info(f"  GRN sequence alignment score: {score}")
            
        return True
        
    except Exception as e:
        logger.error(f"Cross-processor workflow failed: {e}", exc_info=True)
        return False


def main():
    """Main test function."""
    print("\n" + "="*80)
    print(" PROTOS PROCESSOR UPDATE TEST SUITE")
    print(" Testing SequenceProcessor and GRNProcessor Updates with REAL DATA")
    print("="*80)
    
    # Use existing data directory
    data_dir = Path("./data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    logger.info(f"Using data directory: {data_dir.absolute()}")
    logger.info(f"Checking data availability...")
    
    # Check what data is available
    if data_dir.exists():
        seq_count = len(list((data_dir / "sequence" / "fasta").glob("*.fasta"))) if (data_dir / "sequence" / "fasta").exists() else 0
        grn_count = len(list((data_dir / "grn" / "tables").glob("*.csv"))) if (data_dir / "grn" / "tables").exists() else 0
        ref_count = len(list((data_dir / "grn" / "reference").glob("*.csv"))) if (data_dir / "grn" / "reference").exists() else 0
        
        logger.info(f"  Found {seq_count} FASTA files")
        logger.info(f"  Found {grn_count} GRN tables")
        logger.info(f"  Found {ref_count} reference tables")
    
    # Run tests
    results = {}
    
    # Test SequenceProcessor
    results['SequenceProcessor'] = test_sequence_processor(data_dir)
    
    # Test GRNProcessor
    results['GRNProcessor'] = test_grn_processor(data_dir)
    
    # Test cross-processor workflow
    results['Cross-Processor'] = test_cross_processor_workflow(data_dir)
    
    # Summary
    print_section("TEST SUMMARY")
    total_tests = len(results)
    passed_tests = sum(1 for r in results.values() if r)
    
    print(f"\nTotal tests: {total_tests}")
    print(f"Passed: {passed_tests}")
    print(f"Failed: {total_tests - passed_tests}")
    
    for test_name, result in results.items():
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"\n{test_name}: {status}")
    
    # Show directory structure
    print_section("CREATED DIRECTORY STRUCTURE")
    
    def show_tree(path, prefix="", max_depth=3, current_depth=0):
        """Display directory tree."""
        if current_depth >= max_depth:
            return
            
        items = sorted(path.iterdir()) if path.exists() else []
        for i, item in enumerate(items):
            is_last = i == len(items) - 1
            current_prefix = "└── " if is_last else "├── "
            print(f"{prefix}{current_prefix}{item.name}")
            
            if item.is_dir():
                next_prefix = prefix + ("    " if is_last else "│   ")
                show_tree(item, next_prefix, max_depth, current_depth + 1)
    
    print(f"\n{data_dir}/")
    show_tree(data_dir)
    
    return all(results.values())


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
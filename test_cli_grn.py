#!/usr/bin/env python3
"""
Test script for GRN CLI functionality.
"""

import os
import sys
from pathlib import Path

# Set up environment
data_root = "data/"
os.environ["PROTOS_DATA_ROOT"] = data_root
print(f"Data root: {data_root}")

# Test 1: Direct GRNBaseProcessor loading
print("\n=== Test 1: Loading GRN reference table ===")
try:
    from protos.processing.grn.grn_base_processor import GRNBaseProcessor
    
    # Try loading without preload first
    processor = GRNBaseProcessor(
        name="test_processor",
        data_root=data_root,
        processor_data_dir="grn"
    )
    print("✓ Processor created")
    
    # Now try loading the table
    processor.load_grn_table("ref/mo_ref")
    print(f"✓ Loaded GRN table with {len(processor.data)} proteins and {len(processor.data.columns)} GRN positions")
    print(f"  First 5 proteins: {processor.data.index[:5].tolist()}")
    print(f"  First 5 GRNs: {processor.data.columns[:5].tolist()}")
    
except Exception as e:
    print(f"✗ Error loading GRN table: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Sequence extraction
print("\n=== Test 2: Sequence extraction ===")
try:
    seq_dict = processor.get_seq_dict()
    print(f"✓ Extracted {len(seq_dict)} sequences")
    for protein_id, seq in list(seq_dict.items())[:3]:
        print(f"  {protein_id}: {seq[:50]}...")
        
except Exception as e:
    print(f"✗ Error extracting sequences: {type(e).__name__}: {e}")

# Test 3: Check FASTA files
print("\n=== Test 3: Available FASTA files ===")
fasta_dir = Path(data_root) / "sequence" / "fasta" / "processed"
if fasta_dir.exists():
    fasta_files = list(fasta_dir.glob("*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files:")
    for f in fasta_files:
        print(f"  - {f.name}")
else:
    print(f"✗ FASTA directory not found: {fasta_dir}")

# Test 4: Try loading a FASTA file
print("\n=== Test 4: Loading FASTA sequences ===")
try:
    from protos.io.fasta_utils import read_fasta
    
    test_fasta = fasta_dir / "mo_test.fasta"
    if test_fasta.exists():
        sequences = read_fasta(str(test_fasta))
        print(f"✓ Loaded {len(sequences)} sequences from mo_test.fasta")
        for seq_id, seq in sequences.items():
            print(f"  {seq_id}: {len(seq)} residues")
    else:
        print(f"✗ Test FASTA not found: {test_fasta}")
        
except Exception as e:
    print(f"✗ Error loading FASTA: {type(e).__name__}: {e}")

# Test 5: Check if MMseqs2 is available
print("\n=== Test 5: MMseqs2 availability ===")
try:
    import subprocess
    result = subprocess.run(["mmseqs", "--help"], capture_output=True, text=True)
    if result.returncode == 0:
        print("✓ MMseqs2 is available")
    else:
        print("✗ MMseqs2 not found or error running it")
except FileNotFoundError:
    print("✗ MMseqs2 not installed or not in PATH")
except Exception as e:
    print(f"✗ Error checking MMseqs2: {type(e).__name__}: {e}")

# Test 6: Test diagnose CLI
print("\n=== Test 6: Testing diagnose_grn CLI ===")
try:
    from protos.cli.grn.diagnose_grn import diagnose_grn_table
    
    results = diagnose_grn_table(
        grn_table_path="ref/mo_ref",
        protein_family="microbial_opsins",
        use_processor=True,
        data_root=data_root
    )
    
    print(f"✓ Diagnostic completed:")
    print(f"  - Proteins: {results['protein_count']}")
    print(f"  - GRN positions: {results['grn_count']}")
    print(f"  - Invalid GRNs: {results['summary']['invalid_grn_count']}")
    print(f"  - Schiff base issues: {results['summary']['schiff_base_issue_count']}")
    
except Exception as e:
    print(f"✗ Error running diagnostics: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

# Test 7: Test visualize CLI
print("\n=== Test 7: Testing visualize_grn CLI ===")
try:
    from protos.cli.grn.visualize_grn import visualize_grn_table
    
    output_path = Path(data_root) / "output" / "test_heatmap.png"
    output_path.parent.mkdir(exist_ok=True)
    
    # Try without matplotlib requirement
    print("  Note: Visualization requires matplotlib")
    print(f"  Would save to: {output_path}")
    
except Exception as e:
    print(f"✗ Error with visualization: {type(e).__name__}: {e}")

print("\n=== Summary ===")
print(f"Data directory is set up at: {data_root}")
print(f"To run CLI commands, use:")
print(f"  export PROTOS_DATA_ROOT={data_root}")
print(f"  python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t ref/mo_ref")
print(f"  python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_test")
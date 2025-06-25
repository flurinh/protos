"""
Diagnostic tools for GRN assignment debugging.

This module provides tools for debugging and visualizing GRN assignment results.
"""

import os
import sys
import logging
import argparse
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Union, Optional, Any

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string
)

# Configure logging
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

def diagnose_grn_table(
    grn_table_path: str,
    protein_family: str,
    output_path: Optional[str] = None,
    diagnose_loops: bool = True,
    diagnose_tails: bool = True,
    check_schiff_base: bool = True,
    use_processor: bool = True,
    data_root: Optional[str] = None
) -> Dict[str, Any]:
    """
    Diagnose issues in a GRN table.
    
    Args:
        grn_table_path: Dataset ID (if use_processor) or file path
        protein_family: Protein family (e.g., 'gpcr_a', 'microbial_opsins')
        output_path: Path to save diagnostic report
        diagnose_loops: Whether to diagnose loop regions
        diagnose_tails: Whether to diagnose N and C-terminal regions
        check_schiff_base: Whether to check for Schiff base lysine
        use_processor: Whether to use GRNBaseProcessor for loading
        data_root: Root data directory (if None, uses PROTOS_DATA_ROOT)
        
    Returns:
        Dictionary with diagnostic results
    """
    # Load the GRN table
    logger.info(f"Loading GRN table from {grn_table_path}")
    
    if use_processor:
        # Use GRNBaseProcessor
        if data_root is None:
            data_root = os.environ.get('PROTOS_DATA_ROOT', 'data')
        
        processor = GRNBaseProcessor(
            name='diagnose_grn',
            data_root=data_root,
            processor_data_dir='grn',
            dataset=grn_table_path,
            preload=True
        )
        df = processor.data
        
        # Get actual column names (may be dot or x notation)
        columns = df.columns.tolist()
    else:
        # Load directly from file
        df = pd.read_csv(grn_table_path, index_col=0)
        columns = df.columns.tolist()
    
    # Check if table is empty
    if df.empty:
        logger.error("GRN table is empty")
        return {"status": "error", "message": "GRN table is empty"}
    
    results = {
        "status": "success",
        "table_path": grn_table_path,
        "protein_count": len(df),
        "grn_count": len(df.columns),
        "invalid_grns": [],
        "normalized_grns": [],
        "loop_issues": [],
        "tail_issues": [],
        "schiff_base_issues": []
    }
    
    # Check GRN format validity
    logger.info("Checking GRN format validity")
    for column in df.columns:
        is_valid, message = validate_grn_string(column)
        if not is_valid:
            results["invalid_grns"].append({"grn": column, "message": message})
        elif "normalization" in message.lower():
            normalized = normalize_grn_format(column)
            results["normalized_grns"].append({"grn": column, "normalized": normalized})
    
    # Check loop region consistency
    if diagnose_loops:
        logger.info("Diagnosing loop regions")
        loop_grns = [col for col in df.columns if '.' in col and col[0].isdigit() and col[1].isdigit()]
        
        for loop_grn in loop_grns:
            # Split into components
            parts = loop_grn.split('.')
            if len(parts) != 2:
                results["loop_issues"].append({
                    "grn": loop_grn,
                    "issue": "Invalid loop format, expected AB.CCC"
                })
                continue
            
            helix_pair, distance_str = parts
            
            # Check if helices are in correct range
            if len(helix_pair) != 2 or not helix_pair.isdigit():
                results["loop_issues"].append({
                    "grn": loop_grn,
                    "issue": "Invalid helix pair format, expected two digits"
                })
                continue
            
            closer_helix = int(helix_pair[0])
            further_helix = int(helix_pair[1])
            
            if not (1 <= closer_helix <= 8) or not (1 <= further_helix <= 8):
                results["loop_issues"].append({
                    "grn": loop_grn,
                    "issue": f"Helix values out of range: {closer_helix}, {further_helix}"
                })
            
            # Check distance format
            if len(distance_str) != 3 or not distance_str.isdigit():
                results["loop_issues"].append({
                    "grn": loop_grn,
                    "issue": f"Invalid distance format: {distance_str}, expected 3 digits"
                })
    
    # Check N/C terminal regions
    if diagnose_tails:
        logger.info("Diagnosing terminal regions")
        n_term_grns = [col for col in df.columns if col.startswith('n.')]
        c_term_grns = [col for col in df.columns if col.startswith('c.')]
        
        # Validate format
        for n_grn in n_term_grns:
            if not n_grn[2:].isdigit():
                results["tail_issues"].append({
                    "grn": n_grn,
                    "issue": "Invalid N-terminal format, expected n.XX"
                })
        
        for c_grn in c_term_grns:
            if not c_grn[2:].isdigit():
                results["tail_issues"].append({
                    "grn": c_grn,
                    "issue": "Invalid C-terminal format, expected c.XX"
                })
    
    # Check for Schiff base lysine
    if check_schiff_base:
        logger.info("Checking for Schiff base lysine")
        # Handle both x and dot notation
        if protein_family == 'microbial_opsins':
            # Look for both notations
            schiff_col = None
            if '7x50' in df.columns:
                schiff_col = '7x50'
            elif '7.50' in df.columns:
                schiff_col = '7.50'
            
            if schiff_col:
                missing_k = df[~df[schiff_col].str.contains('K', na=False)].index.tolist()
                if missing_k:
                    results["schiff_base_issues"] = [{
                        "protein_id": protein_id,
                        "issue": f"Missing lysine at Schiff base position {schiff_col}",
                        "actual": df.loc[protein_id, schiff_col] if protein_id in df.index else "N/A"
                    } for protein_id in missing_k]
        elif protein_family == 'gpcr_a':
            # Look for both notations
            schiff_col = None
            if '7x43' in df.columns:
                schiff_col = '7x43'
            elif '7.43' in df.columns:
                schiff_col = '7.43'
                
            if schiff_col:
                missing_k = df[~df[schiff_col].str.contains('K', na=False)].index.tolist()
                if missing_k:
                    results["schiff_base_issues"] = [{
                        "protein_id": protein_id,
                        "issue": f"Missing lysine at conserved position {schiff_col}",
                        "actual": df.loc[protein_id, schiff_col] if protein_id in df.index else "N/A"
                    } for protein_id in missing_k]
    
    # Create a summary with issue counts
    results["summary"] = {
        "invalid_grn_count": len(results["invalid_grns"]),
        "normalized_grn_count": len(results["normalized_grns"]),
        "loop_issue_count": len(results["loop_issues"]),
        "tail_issue_count": len(results["tail_issues"]),
        "schiff_base_issue_count": len(results["schiff_base_issues"])
    }
    
    # Save diagnostic report if output path is provided
    if output_path:
        import json
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Diagnostic report saved to {output_path}")
    
    return results

def print_diagnostic_summary(results: Dict[str, Any]):
    """
    Print a summary of diagnostic results.
    
    Args:
        results: Diagnostic results from diagnose_grn_table
    """
    print("\n======= GRN Table Diagnostic Summary =======")
    print(f"Table: {results['table_path']}")
    print(f"Proteins: {results['protein_count']}")
    print(f"GRNs: {results['grn_count']}")
    print("\n--- Issues found ---")
    
    if results['invalid_grns']:
        print(f"Invalid GRNs: {len(results['invalid_grns'])}")
        for issue in results['invalid_grns'][:5]:  # Show first 5 for brevity
            print(f"  - {issue['grn']}: {issue['message']}")
        if len(results['invalid_grns']) > 5:
            print(f"  - ... and {len(results['invalid_grns']) - 5} more")
    else:
        print("Invalid GRNs: None")
        
    if results['normalized_grns']:
        print(f"Normalized GRNs: {len(results['normalized_grns'])}")
        for issue in results['normalized_grns'][:5]:  # Show first 5 for brevity
            print(f"  - {issue['grn']} -> {issue['normalized']}")
        if len(results['normalized_grns']) > 5:
            print(f"  - ... and {len(results['normalized_grns']) - 5} more")
    else:
        print("Normalized GRNs: None")
        
    if results['loop_issues']:
        print(f"Loop region issues: {len(results['loop_issues'])}")
        for issue in results['loop_issues'][:5]:  # Show first 5 for brevity
            print(f"  - {issue['grn']}: {issue['issue']}")
        if len(results['loop_issues']) > 5:
            print(f"  - ... and {len(results['loop_issues']) - 5} more")
    else:
        print("Loop region issues: None")
        
    if results['tail_issues']:
        print(f"Terminal region issues: {len(results['tail_issues'])}")
        for issue in results['tail_issues'][:5]:  # Show first 5 for brevity
            print(f"  - {issue['grn']}: {issue['issue']}")
        if len(results['tail_issues']) > 5:
            print(f"  - ... and {len(results['tail_issues']) - 5} more")
    else:
        print("Terminal region issues: None")
        
    if results['schiff_base_issues']:
        print(f"Schiff base issues: {len(results['schiff_base_issues'])}")
        for issue in results['schiff_base_issues'][:5]:  # Show first 5 for brevity
            print(f"  - {issue['protein_id']}: {issue['issue']}")
        if len(results['schiff_base_issues']) > 5:
            print(f"  - ... and {len(results['schiff_base_issues']) - 5} more")
    else:
        print("Schiff base issues: None")
        
    print("\n===========================================")

def main():
    """Main function when run as script."""
    parser = argparse.ArgumentParser(
        description='Diagnose issues in a GRN table',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Diagnose a GRN table using dataset ID
  python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t ref/mo_ref
  
  # Diagnose with file path
  python -m protos.cli.grn.diagnose_grn -p gpcr_a -t grn_table.csv --use-file
  
  # Save diagnostic report
  python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t mo_table -o report.json
  
  # Skip certain checks
  python -m protos.cli.grn.diagnose_grn -p gpcr_a -t ref_table --no-loops --no-schiff
        """
    )
    parser.add_argument('-p', '--protein_family', required=True, 
                        help='Protein family (e.g., gpcr_a, microbial_opsins)')
    parser.add_argument('-t', '--table', required=True, 
                        help='Dataset ID (or file path if --use-file)')
    parser.add_argument('-o', '--output', 
                        help='Path to save diagnostic report')
    parser.add_argument('--no-loops', action='store_true',
                        help='Skip loop region diagnostics')
    parser.add_argument('--no-tails', action='store_true',
                        help='Skip terminal region diagnostics')
    parser.add_argument('--no-schiff', action='store_true',
                        help='Skip Schiff base checking')
    parser.add_argument('--use-file', action='store_true',
                        help='Treat table argument as file path instead of dataset ID')
    parser.add_argument('--data-root', type=str, default=None,
                        help='Root data directory (default: uses PROTOS_DATA_ROOT env var)')
    
    args = parser.parse_args()
    
    # Set environment variable if data root provided
    if args.data_root:
        os.environ['PROTOS_DATA_ROOT'] = os.path.abspath(args.data_root)
    
    results = diagnose_grn_table(
        grn_table_path=args.table,
        protein_family=args.protein_family,
        output_path=args.output,
        diagnose_loops=not args.no_loops,
        diagnose_tails=not args.no_tails,
        check_schiff_base=not args.no_schiff,
        use_processor=not args.use_file,
        data_root=args.data_root
    )
    
    print_diagnostic_summary(results)
    
    # Return success if no critical issues
    return 0 if len(results['invalid_grns']) == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
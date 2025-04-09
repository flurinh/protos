"""
Implementation of the expand_annotation function for GRN assignment.

This module provides an implementation of the expand_annotation function,
which is used to expand GRN annotations to cover an entire protein sequence.
"""

import logging
from typing import Dict, List, Tuple, Union, Optional, Any
import pandas as pd

# Configure logging
logger = logging.getLogger(__name__)

def expand_annotation(
    new_row: pd.Series, 
    query_seq: str, 
    alignment: Tuple[str, str, str], 
    protein_family: str,
    max_alignment_gap: int = 1, 
    verbose: int = 0
) -> Tuple[List[str], List[str], List[int]]:
    """
    Expand GRN annotations to the full query sequence.
    
    This function takes an initial set of aligned GRNs and expands them to
    annotate the entire query sequence, including loops, tails, and gaps.
    
    Args:
        new_row: Series with initial GRN annotations
        query_seq: Query protein sequence
        alignment: Alignment tuple from format_alignment
        protein_family: Protein family (e.g., 'gpcr_a', 'microbial_opsins')
        max_alignment_gap: Maximum allowed gap in alignment
        verbose: Verbosity level
        
    Returns:
        Tuple of (grns, rns, missing) containing:
            - grns: List of GRN annotations
            - rns: List of residue position labels
            - missing: List of positions without GRN assignment
    """
    try:
        # Import required modules and functions
        from protos.processing.schema.grn_utils_updated import (
            parse_grn_str2float, 
            parse_grn_float2str,
            normalize_grn_format
        )
        from protos.processing.grn.grn_assignment import (
            calculate_missing_gene_numbers,
            assign_missing_std_grns,
            annotate_gaps_and_loops
        )
        from protos.processing.grn.grn_utils import GRNConfigManager
        
        # Initialize config for the protein family
        config_manager = GRNConfigManager(protein_family=protein_family)
        grn_config = config_manager.get_config(strict=False)  # Use non-strict for full coverage
        
        # Step 1: Create query gene numbers (e.g., "A1", "G2", etc.)
        # Convert query sequence to a list of residues with positions
        all_gene_numbers = [f"{res}{i+1}" for i, res in enumerate(query_seq)]
        all_query_gene_numbers = all_gene_numbers.copy()
        
        # Step 2: Extract aligned GRNs from alignment
        # Get dictionary with aligned GRNs from reference by analyzing the alignment
        if verbose > 0:
            logger.info("Getting correctly aligned GRNs from alignment")
        
        # Prepare reference GRN dict from new_row
        reference_grn_dict = new_row.to_dict()
        
        # Extract aligned GRNs from alignment
        try:
            from protos.processing.schema.interface_definitions import get_correctly_aligned_grns
            aligned_grns = get_correctly_aligned_grns(
                all_query_gene_numbers=all_query_gene_numbers,
                reference_grn_dict=reference_grn_dict,
                alignment=alignment,
                max_alignment_gap=max_alignment_gap
            )
        except ImportError:
            # If function not found in interface_definitions, try from grn_table_utils
            from protos.processing.grn.grn_table_utils import get_correctly_aligned_grns
            aligned_grns = get_correctly_aligned_grns(
                all_query_gene_numbers=all_query_gene_numbers,
                reference_grn_dict=reference_grn_dict,
                alignment=alignment,
                max_alignment_gap=max_alignment_gap
            )
        
        if not aligned_grns:
            logger.warning("No GRNs could be aligned. Check if sequences are similar enough.")
            return [], [], list(range(1, len(query_seq) + 1))
            
        # Step 3: Identify missing positions
        missing_gene_numbers = calculate_missing_gene_numbers(
            all_gene_numbers=all_gene_numbers,
            aligned_grns=aligned_grns
        )
        
        if verbose > 0:
            logger.info(f"Found {len(missing_gene_numbers)} missing positions")
        
        # Step 4: Prepare data for further processing
        query_gene_len = len(query_seq)
        
        # Create a list of (gene_number, grn) pairs for aligned positions
        present_seq_nr_grn_list = [(gene_num, grn) for gene_num, grn in 
                                 zip(list(aligned_grns.keys()), list(aligned_grns.values()))]
        
        # Convert GRNs to float for sorting and processing
        grns_float = [parse_grn_str2float(grn) for grn in aligned_grns.values()]
        
        # Convert to string format for consistency
        grns_str = [parse_grn_float2str(grn_float) for grn_float in grns_float]
        
        # Step 5: Annotate N-terminal region
        try:
            from protos.processing.schema.interface_definitions import calculate_missing_ntail_grns
            n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(
                aligned_grns=aligned_grns,
                missing_gene_numbers=missing_gene_numbers,
                grns_float=grns_float
            )
        except ImportError:
            # If function not found in interface_definitions, implement here
            n_tail_list = []
            first_gene_number_int = min([int(gn[1:]) for gn in aligned_grns.keys() if isinstance(gn, str) and gn[0].isalpha()], default=1)
        
        if verbose > 0 and n_tail_list:
            logger.info(f"Annotated {len(n_tail_list)} N-terminal residues")
        
        # Step 6: Annotate C-terminal region
        try:
            from protos.processing.schema.interface_definitions import calculate_missing_ctail_grns
            c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(
                aligned_grns=aligned_grns,
                missing_gene_numbers=missing_gene_numbers,
                query_gene_len=query_gene_len,
                grns_float=grns_float
            )
        except ImportError:
            # If function not found in interface_definitions, implement here
            c_tail_list = []
            last_gene_number_int = max([int(gn[1:]) for gn in aligned_grns.keys() if isinstance(gn, str) and gn[0].isalpha()], default=query_gene_len)
        
        if verbose > 0 and c_tail_list:
            logger.info(f"Annotated {len(c_tail_list)} C-terminal residues")
        
        # Step 7: Identify missing standard GRNs
        # Extract only standard GRNs (e.g., 1x50, 2x50)
        std_grns = [g for g in grns_str if ('x' in g and g[0] != '0' and len(g.split('x')[0]) == 1)]
        missing_std_grns = [g for g in std_grns if g not in aligned_grns.values()]
        
        # Assign missing standard GRNs where possible
        std_grns_assignment, missing = assign_missing_std_grns(
            missing_std_grns=missing_std_grns,
            present_seq_nr_grn_list=present_seq_nr_grn_list,
            query_seq=query_seq,
            missing=missing_gene_numbers,
            grns_str=grns_str
        )
        
        if verbose > 0:
            logger.info(f"Assigned {len(std_grns_assignment)} missing standard GRNs")
            logger.info(f"Still have {len(missing)} missing positions")
        
        # Step 8: Annotate loops and gaps
        nloop, gaps, cloop = annotate_gaps_and_loops(
            present_seq_nr_grn_list=present_seq_nr_grn_list,
            missing=missing,
            query_seq=query_seq,
            grn_config=grn_config,
            grns_str=grns_str
        )
        
        if verbose > 0:
            logger.info(f"Annotated {len(nloop)} N-terminal loop residues")
            logger.info(f"Annotated {len(gaps)} gap residues")
            logger.info(f"Annotated {len(cloop)} C-terminal loop residues")
        
        # Step 9: Combine all regions and sort
        all_grns = n_tail_list + list(aligned_grns.items()) + std_grns_assignment + \
                   gaps + nloop + cloop + c_tail_list
        
        # Step 10: Normalize all GRN formats for consistency
        normalized_grns = []
        for gene_num, grn in all_grns:
            try:
                from protos.processing.schema.grn_utils_updated import normalize_grn_format
                normalized_grn = normalize_grn_format(grn)
            except ImportError:
                # Simple fallback normalization if module not available
                # Convert 1.50 to 1x50 for standard format
                if '.' in grn and len(grn.split('.')[0]) == 1:
                    helix, pos = grn.split('.')
                    normalized_grn = f"{helix}x{int(pos):02d}"
                else:
                    normalized_grn = grn
            normalized_grns.append((gene_num, normalized_grn))
        
        # Extract results
        grns = [g[1] for g in normalized_grns]
        rns = [g[0] for g in normalized_grns]
        
        if verbose > 0:
            logger.info(f"Final annotation has {len(grns)} positions with GRNs and {len(missing)} missing positions")
        
        return grns, rns, missing
    
    except Exception as e:
        logger.error(f"Error in expand_annotation: {e}")
        logger.info("Falling back to legacy implementation")
        
        # Fall back to legacy implementation if there's an error
        try:
            from protos.processing.grn.grn_table_utils import expand_annotation as legacy_expand_annotation
            
            # Convert new_row to the format expected by legacy_expand_annotation
            new_row_dict = new_row.to_dict()
            new_row_list = [(k, v) for k, v in new_row_dict.items()]
            
            # Call the legacy function
            return legacy_expand_annotation(new_row_list, query_seq, alignment, protein_family, max_alignment_gap, verbose)
        except Exception as fallback_error:
            logger.error(f"Legacy fallback also failed: {fallback_error}")
            return [], [], list(range(1, len(query_seq) + 1))
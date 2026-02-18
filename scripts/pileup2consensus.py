#!/usr/bin/env python
"""
Generates a consensus sequence by applying substitutions and multi-base indels 
from a samtools mpileup summary file to a reference FASTA. This version uses a 
more accurate method for calling indels by comparing their support to the 
coverage of the affected downstream region.

Written by Kirsty Kwok, last updated on 2025/09/26 in Glasgow
Usage: python pileup2consensus.py <reference.fasta> <samtools_mpileup_summary.txt> -o <output_directory>
"""

import argparse
import os
import sys
import pyfastx
from collections import defaultdict

# A dictionary to map a combination of nucleotide bases to their corresponding IUPAC code.
# For example, 'AG' represents a position that could be Adenine or Guanine, coded as 'R'.
iupac_code = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N'}

def _parse_indel_string(indel_str: str):
    """
    ## Helper Function: Parse Indel Information
    
    This function takes the raw indel string from the mpileup summary 
    (e.g., '24:TTG|4:AG') and converts it into a more usable Python dictionary.
    
    For '24:TTG|4:AG', it returns: {'TTG': 24, 'AG': 4}
    """
    if not indel_str or indel_str.upper() == 'NA':
        return {}
    
    variants = {}
    parts = indel_str.split('|')
    for part in parts:
        try:
            count, seq = part.split(':', 1)
            variants[seq] = int(count)
        except ValueError:
            continue
    return variants

def get_downstream_coverage(pos: int, length: int, pileup_data: dict):
    """
    ## Helper Function: Look Ahead for Coverage
    
    For a potential deletion, calculates the *average* coverage of all the 
    reference bases that would be removed.
    """
    total_coverage = 0
    actual_positions_checked = 0
    for i in range(length):
        future_pos = pos + i
        if future_pos in pileup_data:
            future_info = pileup_data[future_pos]
            counts = list(map(int, future_info[3:11]))
            total_coverage += sum(counts)
            actual_positions_checked += 1
    
    return total_coverage / actual_positions_checked if actual_positions_checked > 0 else 0

def consensus_call(position_info: list, pileup_data: dict, min_total_depth: int):
    """
    ## Core Logic: The Decision-Making Engine
    
    This function decides the most likely event (Substitution, Insertion, or Deletion)
    at a single genomic position.
    """
    # --- 1. GATHER EVIDENCE ---
    pos = int(position_info[1])
    ref_base = position_info[2]
    counts = list(map(int, position_info[3:11]))
    
    base_counts = {
        'A': counts[0] + counts[4], 'T': counts[1] + counts[5],
        'C': counts[2] + counts[6], 'G': counts[3] + counts[7]
    }
    substitution_support = sum(base_counts.values())

    ins_variants = _parse_indel_string(position_info[11])
    del_variants = _parse_indel_string(position_info[12])
    insertion_support = sum(ins_variants.values())
    deletion_support = sum(del_variants.values())
    
    total_coverage_at_pos = substitution_support + insertion_support + deletion_support
    
    if total_coverage_at_pos < min_total_depth:
        return ('SUB', 'N')

    # --- 2. COMPARE EVENTS FAIRLY ---
    dominant_ins, _ = sorted(ins_variants.items(), key=lambda x:x[1], reverse=True)[0] if ins_variants else (None, 0)
    dominant_del_seq, _ = sorted(del_variants.items(), key=lambda x:x[1], reverse=True)[0] if del_variants else (None, 0)

    del_downstream_coverage = get_downstream_coverage(pos + 1, len(dominant_del_seq) if dominant_del_seq else 0, pileup_data)

    event_scores = []
    if insertion_support > 0 and total_coverage_at_pos > 0:
        event_scores.append((insertion_support / total_coverage_at_pos, 'INS', dominant_ins))
    if deletion_support > 0 and del_downstream_coverage > 0:
        event_scores.append((deletion_support / del_downstream_coverage, 'DEL', dominant_del_seq))
    
    if total_coverage_at_pos > 0:
        event_scores.append((substitution_support / total_coverage_at_pos, 'SUB', None))

    # --- 3. DECLARE THE WINNER ---
    if not event_scores:
        return ('SUB', ref_base)
        
    _, winning_type, winning_value = sorted(event_scores, key=lambda x: x[0], reverse=True)[0]
    
    # --- 4. FORMAT THE OUTPUT ---
    if winning_type == 'INS':
        return ('INS', winning_value)
    elif winning_type == 'DEL':
        return ('DEL', winning_value)
    else: # If the winner is substitution, perform the new "perfect tie" logic.
        
        # --- NEW LOGIC START ---
        
        sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
        
        # Handle cases with no or only one base type
        if len(sorted_bases) < 2 or substitution_support == 0:
            top_base = sorted_bases[0][0] if sorted_bases else ref_base
            return ('SUB', top_base)

        # Get the top two bases
        top_base, top_count = sorted_bases[0]
        second_base, second_count = sorted_bases[1]

        # Calculate their percentage support, rounded to the nearest integer
        top_percent = round((top_count / substitution_support) * 100)
        second_percent = round((second_count / substitution_support) * 100)
        
        # Check for a "perfect tie" based on the rounded percentages
        if top_percent == second_percent:
            # It's a tie, so use an IUPAC ambiguity code for the top two bases.
            consensus_key = ''.join(sorted([top_base, second_base]))
            iupac_result = iupac_code.get(consensus_key, top_base)
            return ('SUB', iupac_result)
        else:
            # There is a clear majority winner, so take the top base.
            return ('SUB', top_base)
        
        # --- NEW LOGIC END ---

def build_coordinate_mapping(chrom, ref_length, substitutions, insertions, deletions):
    """
    ## NEW: Build Coordinate Mapping

    Creates a mapping from consensus positions to reference (scaffold) positions.

    Key insight: We start with reference sequence and apply the SAME transformations
    as the assembly process, but instead of transforming bases, we transform position labels.

    Returns:
        dict: consensus_pos -> reference_pos mapping
    """
    # Create array where each index represents a position, and value is the reference position
    # Initially: position i maps to reference position i
    position_map = list(range(ref_length + 1))  # [0, 1, 2, 3, ..., ref_length]

    # Collect indel events
    indel_events = []
    for pos, seq in insertions.get(chrom, {}).items():
        indel_events.append((pos, 'INS', seq))
    for pos, length in deletions.get(chrom, {}).items():
        indel_events.append((pos, 'DEL', length))

    # Apply indels in REVERSE order (same as assembly, line 268)
    indel_events.sort(key=lambda x: x[0], reverse=True)

    for pos, event_type, value in indel_events:
        if event_type == 'INS':
            # Insertion at position 'pos': insert the sequence
            # In position_map, we insert len(seq) new positions, all mapping to ref pos 'pos'
            insertion_seq = value
            for i in range(len(insertion_seq)):
                position_map.insert(pos, pos)  # All inserted bases map to position 'pos'

        elif event_type == 'DEL':
            # Deletion starting at position 'pos': remove 'length' bases
            # Note: In pileup2consensus.py line 273, deletion is at pos+1
            # But in our position map, we use pos-1 as start_index (line 272)
            start_index = pos - 1
            end_index = start_index + value
            if start_index < len(position_map):
                del position_map[start_index:end_index]

    # Convert to dict: consensus_pos -> reference_pos
    mapping = {}
    for consensus_pos in range(1, len(position_map)):
        ref_pos = position_map[consensus_pos]
        mapping[consensus_pos] = ref_pos if ref_pos > 0 else None

    return mapping

def main(args):
    """
    ## Main Function: Orchestrator

    This function coordinates the entire process: reading files, calling the decision-
    making function for each position, and assembling the final sequence.
    """
    # --- 1. SETUP ---
    print("INFO: Loading reference FASTA...")
    ref_fasta = pyfastx.Fasta(args.reference_fasta)

    print("INFO: Reading and indexing mpileup summary...")
    pileup_data = defaultdict(dict)
    with open(args.input_mpileup_summary) as f:
        mpileup_lines = [line.strip().split('\t') for line in f if line.strip() and not line.strip().startswith('chr')]
        for line in mpileup_lines:
            if len(line) < 13: continue
            chrom = line[0]
            pos = int(line[1])
            pileup_data[chrom][pos] = line

    consensus_seqs = {seq.name: list(seq.seq) for seq in ref_fasta}
    ref_lengths = {seq.name: len(seq.seq) for seq in ref_fasta}
    substitutions = defaultdict(dict)
    insertions = defaultdict(dict)
    deletions = defaultdict(dict)

    # --- 2. ANALYSIS ---
    print("INFO: Collecting consensus calls using look-ahead logic...")
    for chrom, positions in pileup_data.items():
        if chrom not in consensus_seqs:
            print(f"WARNING: Chromosome '{chrom}' from pileup not found in reference. Skipping.")
            continue
        for pos, position_info in positions.items():
            call_type, value = consensus_call(position_info, pileup_data[chrom], args.min_depth)

            if call_type == 'SUB':
                substitutions[chrom][pos] = value
            elif call_type == 'INS':
                insertions[chrom][pos] = value
            elif call_type == 'DEL':
                deletions[chrom][pos + 1] = len(value)

    # --- 3. ASSEMBLY ---
    print("INFO: Applying changes to reference sequence...")
    for chrom in consensus_seqs:
        for pos, base in substitutions[chrom].items():
            if pos <= len(consensus_seqs[chrom]):
                consensus_seqs[chrom][pos - 1] = base

        indel_events = []
        for pos, seq in insertions[chrom].items():
            indel_events.append((pos, 'INS', seq))
        for pos, length in deletions[chrom].items():
            indel_events.append((pos, 'DEL', length))
            
        for pos, event_type, value in sorted(indel_events, key=lambda x: x[0], reverse=True):
            if event_type == 'INS':
                consensus_seqs[chrom].insert(pos, value)
            elif event_type == 'DEL':
                start_index = pos - 1
                end_index = start_index + value
                if start_index < len(consensus_seqs[chrom]):
                    del consensus_seqs[chrom][start_index:end_index]

    # --- 4. OUTPUT ---
    print("INFO: Writing final consensus FASTA...")
    # Determine output path - if it ends with .fasta/.fa, treat as file; otherwise as directory
    if args.output.endswith(('.fasta', '.fa')):
        output_path = args.output
        # Create parent directory if needed
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        # Treat as directory
        output_path = os.path.join(args.output, "consensus_sequence.fasta")
        if not os.path.exists(args.output):
            os.makedirs(args.output)

    with open(output_path, "w") as f:
        for chrom, seq_list in consensus_seqs.items():
            f.write(f">{chrom}\n")
            f.write("".join(seq_list) + "\n")

    # --- 5. COORDINATE MAPPING OUTPUT ---
    if args.output_mapping:
        print("INFO: Writing coordinate mapping file...")
        mapping_path = args.output_mapping

        # Create parent directory if needed
        mapping_dir = os.path.dirname(mapping_path)
        if mapping_dir and not os.path.exists(mapping_dir):
            os.makedirs(mapping_dir)

        with open(mapping_path, 'w') as f:
            f.write("# Consensus-to-Reference Coordinate Mapping\n")
            f.write("# Generated during consensus calling\n")
            f.write("# Consensus_Pos: Position in consensus sequence (output FASTA)\n")
            f.write("# Reference_Pos: Position in reference/scaffold (BAM reference, pileup coordinates)\n")
            f.write("# Event: Type of variant at this position\n")
            f.write("#\n")
            f.write("# IMPORTANT: Use Reference_Pos when querying BAM/pileup files!\n")
            f.write("# The consensus has indels relative to the reference, causing coordinate shifts.\n")
            f.write("#\n")
            f.write("Chromosome\tConsensus_Pos\tReference_Pos\tConsensus_Base\tEvent\n")

            for chrom in consensus_seqs:
                consensus_seq_str = "".join(consensus_seqs[chrom])
                ref_length = ref_lengths.get(chrom, len(consensus_seq_str))

                # Build mapping for this chromosome
                mapping = build_coordinate_mapping(chrom, ref_length, substitutions, insertions, deletions)

                # Write mapping for each consensus position
                for cons_pos in range(1, len(consensus_seq_str) + 1):
                    ref_pos = mapping.get(cons_pos, cons_pos)  # Default to identity if no mapping
                    cons_base = consensus_seq_str[cons_pos - 1]

                    # Determine event type
                    event = "MATCH"
                    if ref_pos is None:
                        event = "DELETION"
                        ref_pos = "-"
                    elif cons_pos in insertions.get(chrom, {}):
                        event = "INSERTION"
                    elif ref_pos in deletions.get(chrom, {}):
                        event = "DELETION"
                    elif cons_pos in substitutions.get(chrom, {}):
                        event = "SUBSTITUTION"

                    f.write(f"{chrom}\t{cons_pos}\t{ref_pos}\t{cons_base}\t{event}\n")

        print(f"INFO: Coordinate mapping written to: {mapping_path}")

    print("INFO: Done.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Generates a consensus sequence from a samtools mpileup summary. 
                       It calls variants by comparing the evidence score for substitutions, insertions, and deletions. 
                       For deletions, it uses a look-ahead method, comparing deletion-supporting reads to the 
                       average coverage of the region that would be removed.''',
        formatter_class=argparse.RawTextHelpFormatter 
    )
    
    parser.add_argument('reference_fasta', help='Path to the reference FASTA file')
    parser.add_argument('input_mpileup_summary', help='Path to the samtools mpileup summary file')
    parser.add_argument('-o', '--output', default='.',
                        help='Output file path (if ends with .fasta/.fa) or directory (default: current directory)')

    parser.add_argument('-m', '--output-mapping', default=None,
                        help='''Output file for consensus-to-reference coordinate mapping (optional).
This mapping is essential for validation, as it tracks how indels affect coordinate systems.
The BAM/pileup uses reference (scaffold) coordinates, while consensus has indels applied.''')

    parser.add_argument('-d', '--min_depth', type=int, default=10,
                        help='''The minimum total read depth required at a position to attempt a call.
This acts as an initial filter before variant scores are compared. (default: 10)''')
    
    args = parser.parse_args()
    
    main(args)
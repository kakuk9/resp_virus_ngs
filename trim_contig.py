#!/usr/bin/env python
"""
Trim and mask viral contigs based on sequencing depth and generate a summary report.
Last updated: 2025-09-29 in Glasgow
"""

import argparse
import os
import pandas as pd
import pyfastx
import sys

def parse_low_coverage_intervals(depth_file, min_depth):
    """
    Parses a bedGraph-style file to find and merge regions of low coverage.
    """
    print(f"INFO: Reading bedGraph file '{depth_file}' and finding regions with depth < {min_depth}...")
    try:
        df = pd.read_csv(
            depth_file, sep='\t', header=None, usecols=[0, 1, 2, 3],
            names=['contig', 'start', 'end', 'depth'],
            dtype={'contig': str, 'start': int, 'end': int, 'depth': int}
        )
    except (pd.errors.EmptyDataError, ValueError):
        # This case is now largely handled by the os.path.getsize() check in main(),
        # but is kept as a safeguard for malformed (but non-empty) files.
        print("INFO: Depth file is empty or malformed. No contigs to process.")
        return {}

    low_depth_df = df[df['depth'] < min_depth].copy()
    if low_depth_df.empty:
        print("INFO: No regions found below the specified depth threshold.")
        return {}

    low_depth_df = low_depth_df.sort_values(['contig', 'start'])
    low_depth_df['group'] = ((low_depth_df['contig'] != low_depth_df['contig'].shift(1)) |
                             (low_depth_df['start'] != low_depth_df['end'].shift(1))).cumsum()
    
    merged_intervals = low_depth_df.groupby(['contig', 'group']).agg(
        start=('start', 'min'), end=('end', 'max')
    ).reset_index()

    low_coverage_regions = {}
    for _, row in merged_intervals.iterrows():
        contig_name = row['contig']
        if contig_name not in low_coverage_regions:
            low_coverage_regions[contig_name] = []
        low_coverage_regions[contig_name].append((row['start'], row['end']))
    
    print(f"INFO: Found and merged low-coverage regions across {len(low_coverage_regions)} contigs.")
    return low_coverage_regions

def main(args):
    """Main function to orchestrate the trimming, masking, and summarizing process."""
    # --- Pre-flight checks for file existence and emptiness ---
    for f in [args.input_depth, args.input_fasta]:
        if not os.path.exists(f):
            print(f"ERROR: Input file not found at '{f}'", file=sys.stderr)
            sys.exit(1)

    if os.path.getsize(args.input_depth) == 0 or os.path.getsize(args.input_fasta) == 0:
        print("INFO: An input file (FASTA or depth) is empty. Generating empty output files.")
        open(args.output_fasta, 'w').close()  # Create empty fasta
        if args.summary:
            print(f"INFO: Writing empty summary file to '{args.summary}'...")
            summary_columns = [
                'contig_name', 'status', 'original_length', 'final_length',
                'bases_trimmed_5prime', 'bases_trimmed_3prime', 'internal_bases_masked',
                'final_n_count', 'final_n_percent'
            ]
            pd.DataFrame(columns=summary_columns).to_csv(args.summary, sep='\t', index=False)
        print("---")
        print("INFO: Processing complete.")
        sys.exit(0)

    # --- Proceed with processing since files are valid ---
    low_coverage_intervals = parse_low_coverage_intervals(args.input_depth, args.min_depth)

    fasta_in = pyfastx.Fasta(args.input_fasta, build_index=False)

    print(f"INFO: Processing FASTA file '{args.input_fasta}'...")
    print(f"INFO: Writing trimmed/masked contigs to '{args.output_fasta}'...")

    summary_data = [] # List to store summary info for each contig

    with open(args.output_fasta, 'w') as f_out:
        for seq in fasta_in:
            # Handle both indexed (seq object) and non-indexed (tuple) modes
            if isinstance(seq, tuple):
                seq_name, seq_sequence = seq[0], seq[1]
            else:
                seq_name, seq_sequence = seq.name, seq.seq

            original_length = len(seq_sequence)
            # Make a copy of intervals so we can pop from it without affecting other contigs
            intervals = list(low_coverage_intervals.get(seq_name, []))

            sequence_chars = list(seq_sequence)
            
            # --- Initialize stats for the summary ---
            bases_trimmed_5prime = 0
            bases_trimmed_3prime = 0
            internal_bases_masked = 0
            
            trim_start, trim_end = 0, original_length

            # --- Handle 5' Trimming ---
            if intervals and intervals[0][0] == 0:
                trim_start = intervals[0][1]
                bases_trimmed_5prime = trim_start
                intervals.pop(0)

            # --- Handle 3' Trimming ---
            if intervals and intervals[-1][1] == original_length:
                trim_end = intervals[-1][0]
                bases_trimmed_3prime = original_length - trim_end
                intervals.pop(-1)

            # --- Handle Internal Masking ---
            for start, end in intervals:
                mask_len = end - start
                sequence_chars[start:end] = 'N' * mask_len
                internal_bases_masked += mask_len

            # --- Finalize sequence and stats ---
            final_sequence = "".join(sequence_chars[trim_start:trim_end])
            final_length = len(final_sequence)
            
            if final_length > 0:
                final_n_count = final_sequence.count('N')
                final_n_percent = (final_n_count / final_length) * 100
            else:
                final_n_count = 0
                final_n_percent = 0 # An empty sequence has 0% Ns

            # --- Apply filters and determine final status ---
            if final_length < args.min_len:
                status = "discarded_too_short"
            elif final_n_percent >= args.max_n_percent:
                status = "discarded_high_n_content"
            else:
                is_trimmed = bases_trimmed_5prime > 0 or bases_trimmed_3prime > 0
                is_masked = internal_bases_masked > 0
                
                if not is_trimmed and not is_masked:
                    change = "u"  # unchanged
                else:
                    change = ""
                    if is_trimmed:
                        change += "t"
                    if is_masked:
                        change += "m"
                
                status = f"kept_{change}"
                name = f"{seq_name}_{change}_{final_length}"
                f_out.write(f">{name}\n{final_sequence}\n")

            # --- Append all stats to our summary list ---
            summary_data.append({
                'contig_name': seq_name,
                'status': status,
                'original_length': original_length,
                'final_length': final_length,
                'bases_trimmed_5prime': bases_trimmed_5prime,
                'bases_trimmed_3prime': bases_trimmed_3prime,
                'internal_bases_masked': internal_bases_masked,
                'final_n_count': final_n_count,
                'final_n_percent': round(final_n_percent, 2)
            })

    # --- Write the summary file if requested ---
    if args.summary:
        print(f"INFO: Writing summary file to '{args.summary}'...")
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(args.summary, sep='\t', index=False)
    
    print("---")
    print("INFO: Processing complete.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Trim and mask contigs based on sequencing depth and generate a summary report.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_depth', help='Input depth file in bedGraph format (from `bedtools genomecov -bga`)')
    parser.add_argument('input_fasta', help='Input fasta file to be processed')
    parser.add_argument('-o', '--output_fasta', required=True, help='Path to the output fasta file')
    parser.add_argument('-s', '--summary', help='(Optional) Path to the output summary TSV file.')
    parser.add_argument('-d', '--min_depth', help='Minimum depth to keep (default: 10)', default=10, type=int)
    parser.add_argument('-l', '--min_len', help='Minimum contig length to keep after trimming (default: 200)', default=200, type=int)
    parser.add_argument('-n', '--max_n_percent', help='Maximum percentage of Ns allowed (default: 20.0)', default=20.0, type=float)
    args = parser.parse_args()
    main(args)

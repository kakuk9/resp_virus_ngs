#!/usr/bin/env python
#############################################
# Date: 2025-10-03 in Glasgow
# Description:
# A robust Python script to parse a multi-sample samtools mpileup file.
# It correctly handles complex indel and quality string formats.
# The script counts high-quality bases, insertions, and deletions for
# each sample at each genomic position and writes the output to
# separate, tab-delimited files.
#############################################

import sys
import argparse
from collections import Counter

def parse_pileup_entry(
    bases_string: str,
    qualities_string: str,
    ref_base: str,
    bq_cutoff: int,
    offset: int
) -> str:
    """
    Robustly parses the 'bases' and 'qualities' strings from a single pileup column.
    This function correctly handles the complex mpileup format where the length
    of the bases string does not match the qualities string.

    Returns a formatted string of counts. Returns "0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA" if no coverage.
    """
    if bases_string == "*" or not bases_string:
        # Return zero counts for positions with no coverage (needed for -aa flag)
        return "0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA"

    # --- Data structures to hold counts ---
    base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, 'a': 0, 't': 0, 'c': 0, 'g': 0, 'n': 0}
    insertions = Counter()
    deletions = Counter()

    # --- Walk through the bases string, keeping track of position in qualities ---
    base_idx = 0
    qual_idx = 0
    while base_idx < len(bases_string):
        char = bases_string[base_idx]

        if char == '^':  # Signifies start of a read; the next char is mapping quality
            base_idx += 2
            continue
        elif char == '$':  # Signifies end of a read
            base_idx += 1
            continue
        elif char in '+-': # An indel
            # Find the number of bases in the indel (can be multi-digit)
            num_str = ""
            i = base_idx + 1
            while i < len(bases_string) and bases_string[i].isdigit():
                num_str += bases_string[i]
                i += 1

            if not num_str:
                base_idx += 1
                continue # Malformed indel, skip

            num = int(num_str)
            indel_seq = bases_string[i : i + num]

            if char == '+':
                insertions[indel_seq.upper()] += 1
            elif char == '-':
                deletions[indel_seq.upper()] += 1

            # Skip past the entire indel notation in the bases string
            base_idx = i + num
            continue

        # --- If we reach here, it must be a real base call ---
        base_call = char

        # Check if there is a corresponding quality score
        if qual_idx < len(qualities_string):
            score = ord(qualities_string[qual_idx]) - offset

            if score >= bq_cutoff:
                if base_call == '.':
                    base_counts[ref_base.upper()] += 1
                elif base_call == ',':
                    base_counts[ref_base.lower()] += 1
                elif base_call.upper() in 'ATCG':
                    base_counts[base_call] += 1

        # Advance both pointers
        base_idx += 1
        qual_idx += 1

    # --- Format the output string ---
    counts_str = "\t".join(str(c) for c in [
        base_counts['A'], base_counts['T'], base_counts['C'], base_counts['G'],
        base_counts['a'], base_counts['t'], base_counts['c'], base_counts['g']
    ])

    ins_str = "NA"
    if insertions:
        ins_str = "|".join([f"{count}:{seq}" for seq, count in insertions.most_common()])

    del_str = "NA"
    if deletions:
        del_str = "|".join([f"{count}:{seq}" for seq, count in deletions.most_common()])

    return f"{counts_str}\t{ins_str}\t{del_str}"

def main():
    """Main function to orchestrate file parsing and writing."""
    parser = argparse.ArgumentParser(
        description="Parse a multi-sample samtools mpileup file to count bases, insertions, and deletions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input", help="Input pileup file from one or more samples.")
    parser.add_argument("-p", "--prefix", default="sample", help="Output file prefix. Default is 'sample'.")
    parser.add_argument("-out", "--output", default=None,
                        help="Direct output file path (for single-sample pileup). If specified, overrides prefix for sample 1.")
    parser.add_argument("-bq", "--bq_cutoff", type=int, default=0,
                        help="Base quality score cutoff. Default is 0 (no filter).")
    parser.add_argument("-o", "--offset", type=int, default=33,
                        help="ASCII to base quality score offset. Default is 33 (Sanger).")
    args = parser.parse_args()

    try:
        infile = open(args.input, 'r')
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"[{__import__('datetime').datetime.now()}] Begin parsing...")

    # --- Read the first line to determine the number of samples ---
    first_line = infile.readline()
    if not first_line:
        print("Input file is empty.", file=sys.stderr)
        sys.exit(1)

    columns = first_line.strip().split('\t')
    num_samples = (len(columns) - 3) // 3
    if num_samples == 0:
        print("Error: Could not determine number of samples from pileup file.", file=sys.stderr)
        sys.exit(1)

    print(f"Detected {num_samples} sample(s).")

    # --- Open output files and write headers ---
    out_files = {}
    header = "chr\tloc\tref\tA\tT\tC\tG\ta\tt\tc\tg\tInsertion\tDeletion\n"
    for i in range(1, num_samples + 1):
        try:
            # Use --output for first sample if provided, otherwise use prefix
            if i == 1 and args.output:
                output_filename = args.output
            else:
                output_filename = f"{args.prefix}{i}.txt"

            fh = open(output_filename, 'w')
            fh.write(header)
            out_files[i] = fh
        except IOError as e:
            print(f"Error: Could not open output file {output_filename}. Reason: {e}", file=sys.stderr)
            sys.exit(1)

    # --- Process all lines in the file ---
    infile.seek(0) # Rewind to start of file to process the first line again
    for line_num, line in enumerate(infile, 1):
        columns = line.strip().split('\t')
        if len(columns) < 3:
            print(f"Warning: Skipping malformed line {line_num}.", file=sys.stderr)
            continue

        chrom, pos, ref = columns[0:3]

        for i in range(num_samples):
            try:
                sample_cols = columns[3 + i * 3 : 3 + (i + 1) * 3]
                depth, bases, qualities = sample_cols
            except ValueError:
                continue # Not enough columns for this sample on this line

            result_str = parse_pileup_entry(bases, qualities, ref, args.bq_cutoff, args.offset)

            # Write all positions, including those with 0 coverage (for -aa flag support)
            out_files[i + 1].write(f"{chrom}\t{pos}\t{ref}\t{result_str}\n")

    # --- Clean up ---
    infile.close()
    for fh in out_files.values():
        fh.close()

    print(f"[{__import__('datetime').datetime.now()}] Finished.")

if __name__ == "__main__":
    main()

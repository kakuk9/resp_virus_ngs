#!/usr/bin/env python3
"""
filter_nucmer_alignments.py

Parse and filter nucmer alignment coordinates.
Remove self-alignments and within-sample comparisons.

Author: Kirsty
Date: 2026-01-08
"""

import argparse
import pandas as pd


def parse_coords(coords_file):
    """
    Parse nucmer show-coords output (-r -c -l -T format).

    Columns:
    [S1]  [E1]  [S2]  [E2]  [LEN 1]  [LEN 2]  [% IDY]  [LEN R]  [LEN Q]  [COV R]  [COV Q]  [TAGS]
    """
    # Read coords file, skip header lines
    with open(coords_file, 'r') as f:
        lines = f.readlines()

    # Find where data starts (after header lines)
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('='):
            data_start = i + 1
            break

    if data_start == 0:
        # No alignments found
        return pd.DataFrame()

    # Parse data lines
    alignments = []
    for line in lines[data_start:]:
        if not line.strip():
            continue

        parts = line.strip().split('\t')
        if len(parts) < 13:
            continue

        try:
            alignments.append({
                'ref_start': int(parts[0]),
                'ref_end': int(parts[1]),
                'query_start': int(parts[2]),
                'query_end': int(parts[3]),
                'ref_align_len': int(parts[4]),
                'query_align_len': int(parts[5]),
                'percent_identity': float(parts[6]),
                'ref_length': int(parts[7]),
                'query_length': int(parts[8]),
                'ref_coverage': float(parts[9]),
                'query_coverage': float(parts[10]),
                'ref_contig': parts[11],
                'query_contig': parts[12]
            })
        except (ValueError, IndexError):
            continue

    if not alignments:
        return pd.DataFrame()

    return pd.DataFrame(alignments)


def filter_alignments(df, min_identity, min_length):
    """
    Filter alignments based on identity and length thresholds.
    Remove self-alignments and within-sample comparisons.
    """
    if df.empty:
        return df

    # Remove self-alignments (same contig)
    df = df[df['query_contig'] != df['ref_contig']].copy()

    # Extract sample names from contig names (format: sample|contig_id)
    df['query_sample'] = df['query_contig'].str.split('|').str[0]
    df['ref_sample'] = df['ref_contig'].str.split('|').str[0]

    # Remove within-sample comparisons
    df = df[df['query_sample'] != df['ref_sample']].copy()

    # Filter by identity
    df = df[df['percent_identity'] >= min_identity].copy()

    # Filter by alignment length (use the smaller of the two)
    df['align_length'] = df[['ref_align_len', 'query_align_len']].min(axis=1)
    df = df[df['align_length'] >= min_length].copy()

    return df


def main():
    parser = argparse.ArgumentParser(
        description='Filter nucmer alignment coordinates'
    )
    parser.add_argument(
        'coords_file',
        help='Input coords file from show-coords'
    )
    parser.add_argument(
        'output_file',
        help='Output filtered alignments TSV'
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=95.0,
        help='Minimum percent identity (default: 95.0)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum alignment length (default: 100)'
    )

    args = parser.parse_args()

    # Parse coords
    df = parse_coords(args.coords_file)

    if df.empty:
        # Create empty output with header
        pd.DataFrame(columns=[
            'query_contig', 'ref_contig', 'query_sample', 'ref_sample',
            'query_start', 'query_end', 'ref_start', 'ref_end',
            'align_length', 'percent_identity',
            'query_length', 'ref_length', 'query_coverage', 'ref_coverage'
        ]).to_csv(args.output_file, sep='\t', index=False)
        return

    # Filter alignments
    filtered = filter_alignments(df, args.min_identity, args.min_length)

    if filtered.empty:
        # Create empty output with header
        pd.DataFrame(columns=[
            'query_contig', 'ref_contig', 'query_sample', 'ref_sample',
            'query_start', 'query_end', 'ref_start', 'ref_end',
            'align_length', 'percent_identity',
            'query_length', 'ref_length', 'query_coverage', 'ref_coverage'
        ]).to_csv(args.output_file, sep='\t', index=False)
        return

    # Select and reorder columns for output
    output_cols = [
        'query_contig', 'ref_contig', 'query_sample', 'ref_sample',
        'query_start', 'query_end', 'ref_start', 'ref_end',
        'align_length', 'percent_identity',
        'query_length', 'ref_length', 'query_coverage', 'ref_coverage'
    ]

    filtered[output_cols].to_csv(args.output_file, sep='\t', index=False)

    print(f"Filtered {len(filtered)} alignments from {len(df)} total alignments")
    print(f"Thresholds: identity >= {args.min_identity}%, length >= {args.min_length} bp")


if __name__ == '__main__':
    main()

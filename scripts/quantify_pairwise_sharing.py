#!/usr/bin/env python3
"""
Quantify pairwise sharing of identical reads between samples.

Input: Duplicate reads TSV from find_duplicate_reads.py
Output: Pairwise matrix showing shared read counts between sample pairs

Author: Kirsty Kwok
Date: 2026-01-07
"""

import sys
import argparse
import pandas as pd
from itertools import combinations
from collections import defaultdict


def parse_duplicates(duplicates_file):
    """
    Parse duplicates file and extract sample-pair sharing information.

    Returns:
        - pairwise_sharing: dict {(sample_A, sample_B): {'reads': count, 'unique_seqs': count}}
        - sample_totals: dict {sample: total_duplicate_read_count}
    """
    df = pd.read_csv(duplicates_file, sep='\t')

    pairwise_sharing = defaultdict(lambda: {'reads': 0, 'unique_seqs': 0})
    sample_totals = defaultdict(int)

    for _, row in df.iterrows():
        samples = row['sample_list'].split('|')
        counts = list(map(int, row['count_list'].split('|')))

        # Update sample totals
        for sample, count in zip(samples, counts):
            sample_totals[sample] += count

        # For each pair of samples sharing this sequence
        for i, sample_a in enumerate(samples):
            for j, sample_b in enumerate(samples):
                if i < j:  # Avoid duplicates (A-B is same as B-A)
                    # Count the minimum (conservative estimate of shared reads)
                    shared_count = min(counts[i], counts[j])
                    pair = tuple(sorted([sample_a, sample_b]))
                    pairwise_sharing[pair]['reads'] += shared_count
                    pairwise_sharing[pair]['unique_seqs'] += 1  # Each row = 1 unique sequence

    return pairwise_sharing, sample_totals


def create_pairwise_matrix(pairwise_sharing, sample_totals, genus):
    """
    Create a pairwise matrix showing shared reads and unique sequences between samples.

    Returns:
        DataFrame with columns:
        - sample_A, sample_B, genus, shared_read_count, identical_unique_seqs,
        - sample_A_total, sample_B_total,
        - percent_of_A, percent_of_B, percent_of_smaller
    """
    rows = []

    for (sample_a, sample_b), sharing_data in pairwise_sharing.items():
        shared_reads = sharing_data['reads']
        shared_unique_seqs = sharing_data['unique_seqs']

        total_a = sample_totals.get(sample_a, 0)
        total_b = sample_totals.get(sample_b, 0)

        percent_a = (shared_reads / total_a * 100) if total_a > 0 else 0.0
        percent_b = (shared_reads / total_b * 100) if total_b > 0 else 0.0

        smaller_total = min(total_a, total_b)
        percent_smaller = (shared_reads / smaller_total * 100) if smaller_total > 0 else 0.0

        rows.append({
            'sample_A': sample_a,
            'sample_B': sample_b,
            'genus': genus,
            'shared_read_count': shared_reads,
            'identical_unique_seqs': shared_unique_seqs,
            'sample_A_total': total_a,
            'sample_B_total': total_b,
            'percent_of_A': round(percent_a, 2),
            'percent_of_B': round(percent_b, 2),
            'percent_of_smaller': round(percent_smaller, 2)
        })

    df = pd.DataFrame(rows)

    # Sort by shared_read_count (descending)
    if not df.empty:
        df = df.sort_values('shared_read_count', ascending=False)

    return df


def main():
    parser = argparse.ArgumentParser(
        description='Quantify pairwise sharing of identical reads'
    )
    parser.add_argument('input_duplicates', help='Input duplicates TSV file')
    parser.add_argument('output_pairwise', help='Output pairwise matrix TSV file')
    parser.add_argument('--genus', required=True, help='Genus name')

    args = parser.parse_args()

    # Parse duplicates
    print(f"Processing genus: {args.genus}")
    pairwise_sharing, sample_totals = parse_duplicates(args.input_duplicates)

    print(f"Total sample pairs with shared reads: {len(pairwise_sharing)}")

    # Create pairwise matrix
    pairwise_df = create_pairwise_matrix(pairwise_sharing, sample_totals, args.genus)

    # Write output
    pairwise_df.to_csv(args.output_pairwise, sep='\t', index=False)

    print(f"Pairwise matrix written to: {args.output_pairwise}")

    if not pairwise_df.empty:
        print(f"\nTop 5 sample pairs by shared reads:")
        for idx, row in pairwise_df.head(5).iterrows():
            print(f"  {row['sample_A']} <-> {row['sample_B']}: "
                  f"{row['shared_read_count']} reads "
                  f"({row['percent_of_smaller']:.1f}% of smaller sample)")


if __name__ == '__main__':
    main()

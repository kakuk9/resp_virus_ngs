#!/usr/bin/env python3
"""
Find duplicate (identical) reads within a genus and track sample sources.

Input: Combined FASTQ file with sample-prefixed read names (sample_id|read_name)
Output:
  - TSV file listing each unique sequence that appears in 2+ samples
  - Stats file with per-sample read counts

Author: Kirsty Kwok
Date: 2026-01-07
"""

import sys
import gzip
import argparse
from collections import defaultdict


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def canonicalize_sequence(seq):
    """
    Return the lexicographically smaller of the sequence and its reverse complement.
    This ensures that a sequence and its RC are treated as the same.
    """
    rc = reverse_complement(seq)
    return seq if seq < rc else rc


def parse_fastq_gz(fastq_file):
    """
    Parse gzipped FASTQ file and yield (read_name, sequence) tuples.
    """
    with gzip.open(fastq_file, 'rt') as f:
        while True:
            # Read 4 lines per FASTQ entry
            header = f.readline().strip()
            if not header:
                break

            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            # Extract sample ID from read name (format: sample_id|read_name)
            if '|' in header:
                sample_id = header.split('|')[0].lstrip('@')
                read_name = header.split('|', 1)[1]
            else:
                # Fallback if no sample prefix
                sample_id = 'unknown'
                read_name = header.lstrip('@')

            yield sample_id, read_name, sequence


def find_duplicates(fastq_file, genus, use_revcomp=True):
    """
    Find sequences that appear in multiple samples.
    Optionally treats a sequence and its reverse complement as identical.

    Returns:
        - cross_sample_duplicates: dict {sequence: {sample1: count, sample2: count, ...}}
        - sample_totals: dict {sample: total_read_count}
        - within_sample_duplicates: dict {sample: duplicate_read_count}
        - sample_unique_seqs: dict {sample: set of unique sequences}
    """
    # Track sequence -> sample -> count
    seq_to_samples = defaultdict(lambda: defaultdict(int))
    sample_totals = defaultdict(int)

    print(f"Parsing FASTQ file for genus: {genus}")
    if use_revcomp:
        print(f"Treating sequences and their reverse complements as identical")

    for sample_id, read_name, sequence in parse_fastq_gz(fastq_file):
        # Canonicalize sequence (use lexicographically smaller of seq/RC)
        if use_revcomp:
            canonical_seq = canonicalize_sequence(sequence)
        else:
            canonical_seq = sequence

        seq_to_samples[canonical_seq][sample_id] += 1
        sample_totals[sample_id] += 1

    print(f"Total unique sequences: {len(seq_to_samples)}")
    print(f"Total samples: {len(sample_totals)}")

    # Build per-sample unique sequence sets
    sample_unique_seqs = defaultdict(set)
    for seq, sample_counts in seq_to_samples.items():
        for sample in sample_counts.keys():
            sample_unique_seqs[sample].add(seq)

    # Calculate within-sample duplicates (sequences that appear >1 time in same sample)
    within_sample_duplicates = defaultdict(int)
    for seq, sample_counts in seq_to_samples.items():
        for sample, count in sample_counts.items():
            if count > 1:
                # This sequence appears multiple times within this sample
                within_sample_duplicates[sample] += count

    # Filter to only sequences present in 2+ samples (cross-sample duplicates)
    cross_sample_duplicates = {}
    for seq, sample_counts in seq_to_samples.items():
        if len(sample_counts) >= 2:
            cross_sample_duplicates[seq] = sample_counts

    print(f"Cross-sample duplicate sequences (present in 2+ samples): {len(cross_sample_duplicates)}")

    return cross_sample_duplicates, sample_totals, within_sample_duplicates, sample_unique_seqs


def write_duplicates(duplicates, output_file, genus):
    """
    Write duplicate sequences to TSV file.

    Format:
    sequence    genus    sample_list    count_list    total_occurrences    num_samples
    """
    with open(output_file, 'w') as f:
        # Header
        f.write('\t'.join([
            'sequence',
            'genus',
            'sample_list',
            'count_list',
            'total_occurrences',
            'num_samples'
        ]) + '\n')

        # Sort by number of samples (descending), then by total occurrences
        sorted_duplicates = sorted(
            duplicates.items(),
            key=lambda x: (len(x[1]), sum(x[1].values())),
            reverse=True
        )

        for sequence, sample_counts in sorted_duplicates:
            samples = sorted(sample_counts.keys())
            counts = [str(sample_counts[s]) for s in samples]

            f.write('\t'.join([
                sequence,
                genus,
                '|'.join(samples),
                '|'.join(counts),
                str(sum(sample_counts.values())),
                str(len(sample_counts))
            ]) + '\n')


def write_stats(sample_totals, cross_sample_duplicates, within_sample_duplicates,
                sample_unique_seqs, output_file, genus):
    """
    Write per-sample statistics including duplicates and unique sequence metrics.

    Format:
    sample  genus  total_reads  unique_sequences  sample_specific_unique_seqs  shared_unique_seqs
    cross_sample_dup_reads  within_sample_dup_reads  percent_cross_sample  percent_within_sample
    percent_sample_specific  library_complexity
    """
    # Count cross-sample duplicate reads per sample
    cross_sample_dup_counts = defaultdict(int)
    for seq, sample_counts in cross_sample_duplicates.items():
        for sample, count in sample_counts.items():
            cross_sample_dup_counts[sample] += count

    # Determine which sequences are shared vs sample-specific
    sample_specific_seqs = {}
    shared_seqs = {}

    for sample, seqs in sample_unique_seqs.items():
        sample_specific = set()
        shared = set()

        for seq in seqs:
            # Check if this sequence appears in any other sample
            appears_in_others = False
            for other_sample, other_seqs in sample_unique_seqs.items():
                if other_sample != sample and seq in other_seqs:
                    appears_in_others = True
                    break

            if appears_in_others:
                shared.add(seq)
            else:
                sample_specific.add(seq)

        sample_specific_seqs[sample] = sample_specific
        shared_seqs[sample] = shared

    with open(output_file, 'w') as f:
        # Header
        f.write('\t'.join([
            'sample',
            'genus',
            'total_reads',
            'unique_sequences',
            'sample_specific_unique_seqs',
            'shared_unique_seqs',
            'cross_sample_dup_reads',
            'within_sample_dup_reads',
            'percent_cross_sample',
            'percent_within_sample',
            'percent_sample_specific',
            'library_complexity'
        ]) + '\n')

        for sample in sorted(sample_totals.keys()):
            total = sample_totals[sample]
            unique_seqs = len(sample_unique_seqs.get(sample, set()))
            sample_specific = len(sample_specific_seqs.get(sample, set()))
            shared = len(shared_seqs.get(sample, set()))
            cross_dup = cross_sample_dup_counts.get(sample, 0)
            within_dup = within_sample_duplicates.get(sample, 0)

            percent_cross = (cross_dup / total * 100) if total > 0 else 0.0
            percent_within = (within_dup / total * 100) if total > 0 else 0.0
            percent_specific = (sample_specific / unique_seqs * 100) if unique_seqs > 0 else 0.0
            complexity = (unique_seqs / total) if total > 0 else 0.0

            f.write('\t'.join([
                sample,
                genus,
                str(total),
                str(unique_seqs),
                str(sample_specific),
                str(shared),
                str(cross_dup),
                str(within_dup),
                f'{percent_cross:.2f}',
                f'{percent_within:.2f}',
                f'{percent_specific:.2f}',
                f'{complexity:.4f}'
            ]) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Find duplicate reads within a genus across samples'
    )
    parser.add_argument('input_fastq', help='Combined FASTQ file (gzipped)')
    parser.add_argument('output_duplicates', help='Output TSV file for duplicates')
    parser.add_argument('output_stats', help='Output TSV file for statistics')
    parser.add_argument('--genus', required=True, help='Genus name')
    parser.add_argument('--no-revcomp', action='store_true',
                        help='Disable reverse complement checking (default: enabled)')

    args = parser.parse_args()

    # Find duplicates (with or without reverse complement)
    use_revcomp = not args.no_revcomp
    cross_sample_duplicates, sample_totals, within_sample_duplicates, sample_unique_seqs = find_duplicates(
        args.input_fastq, args.genus, use_revcomp
    )

    # Write outputs
    write_duplicates(cross_sample_duplicates, args.output_duplicates, args.genus)
    write_stats(sample_totals, cross_sample_duplicates, within_sample_duplicates,
                sample_unique_seqs, args.output_stats, args.genus)

    print(f"Cross-sample duplicate reads written to: {args.output_duplicates}")
    print(f"Statistics (including unique sequences and library complexity) written to: {args.output_stats}")


if __name__ == '__main__':
    main()

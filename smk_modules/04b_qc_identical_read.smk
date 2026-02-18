"""
04b_qc_identical_read.smk

Targeted read-level validation of identical contigs between samples.
For contig pairs identified as 100% identical in 04a, this workflow:
- Extracts reads that mapped to each contig in the pair
- Compares if the reads themselves are identical between samples
- Provides evidence for true contamination vs. convergent assembly

Author: Kirsty K
Date: 2026-01-09
"""

import pandas as pd
import os
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

BATCH = config.get("batch", "batch")
RUN_DIR = config.get("run_directory", ".")

# Input: identical contig report from 04a
IDENTICAL_CONTIG_REPORT = f"04a_qc/nucmer_output/{BATCH}.identical_contig_report.tsv"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_identical_contigs():
    """
    Parse the identical contig report from 04a to get contig pairs to validate.
    Deduplicates reversed pairs (A vs B and B vs A are the same comparison).
    Returns list of dicts with sample/contig info, lengths, and alignment coordinates.
    """
    if not os.path.exists(IDENTICAL_CONTIG_REPORT):
        print(f"Warning: Identical contig report not found: {IDENTICAL_CONTIG_REPORT}")
        return []

    df = pd.read_csv(IDENTICAL_CONTIG_REPORT, sep='\t')

    if df.empty:
        print("No identical contigs found in 04a report")
        return []

    contig_pairs = []
    seen_pairs = set()  # Track unique pairs to avoid duplicates

    for idx, row in df.iterrows():
        sample1 = row['sample1_id']
        sample2 = row['sample2_id']

        # Clean contig IDs (remove sample prefix)
        contig1_full = row['contig1_id']
        contig2_full = row['contig2_id']

        contig1 = contig1_full.replace(f"{sample1}_", "")
        contig2 = contig2_full.replace(f"{sample2}_", "")

        genus = row['sample1_genus']

        # Get contig lengths from 04a report
        contig1_len = row['contig1_len']
        contig2_len = row['contig2_len']

        # Get alignment coordinates from 04a report
        contig1_aln_start = row['contig1_start']
        contig1_aln_end = row['contig1_end']
        contig2_aln_start = row['contig2_start']
        contig2_aln_end = row['contig2_end']

        # Create a canonical pair key (sorted to catch reversals)
        # e.g., both (A,X,B,Y) and (B,Y,A,X) become the same key
        pair_key = tuple(sorted([(sample1, contig1), (sample2, contig2)]))

        # Skip if we've already seen this pair (in either direction)
        if pair_key in seen_pairs:
            continue

        seen_pairs.add(pair_key)

        # Create unique pair ID for tracking
        pair_id = f"{sample1}_{contig1}_vs_{sample2}_{contig2}"

        contig_pairs.append({
            'pair_id': pair_id,
            'sample1': sample1,
            'contig1': contig1,
            'sample2': sample2,
            'contig2': contig2,
            'genus': genus,
            'contig1_len': contig1_len,
            'contig2_len': contig2_len,
            'contig1_aln_start': contig1_aln_start,
            'contig1_aln_end': contig1_aln_end,
            'contig2_aln_start': contig2_aln_start,
            'contig2_aln_end': contig2_aln_end
        })

    print(f"Found {len(contig_pairs)} unique identical contig pairs to validate (after deduplication)")
    return contig_pairs

# Parse contig pairs
CONTIG_PAIRS = parse_identical_contigs()
PAIR_IDS = [pair['pair_id'] for pair in CONTIG_PAIRS]

# Create lookup dictionaries for wildcards
PAIR_TO_INFO = {pair['pair_id']: pair for pair in CONTIG_PAIRS}

# =============================================================================
# RULE ALL
# =============================================================================

rule all:
    input:
        f"04b_qc/read_validation_report/{BATCH}.read_validation_summary.tsv"

# =============================================================================
# RULE 1: EXTRACT READS MAPPED TO CONTIG1
# =============================================================================

rule extract_reads_contig1:
    """
    Extract read names that mapped to contig1 from sample1's BAM file.
    Uses BAM files directly from untrimmed_viral_contig_bowtie2.
    """
    input:
        bam = lambda wildcards: f"untrimmed_viral_contig_bowtie2/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}.bam"
    output:
        read_names = "04b_qc/contig_reads/{pair_id}.sample1.read_names.txt"
    params:
        contig_id = lambda wildcards: PAIR_TO_INFO[wildcards.pair_id]['contig1']
    threads:
        4
    shell:
        """
        samtools view -F 4 {input.bam} {params.contig_id} -@ {threads}| cut -f1 | sort -u > {output.read_names}
        """

# =============================================================================
# RULE 2: EXTRACT READS MAPPED TO CONTIG2
# =============================================================================

rule extract_reads_contig2:
    """
    Extract read names that mapped to contig2 from sample2's BAM file.
    Uses BAM files directly from untrimmed_viral_contig_bowtie2.
    """
    input:
        bam = lambda wildcards: f"untrimmed_viral_contig_bowtie2/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}.bam"
    output:
        read_names = "04b_qc/contig_reads/{pair_id}.sample2.read_names.txt"
    params:
        contig_id = lambda wildcards: PAIR_TO_INFO[wildcards.pair_id]['contig2']
    threads:
        4
    shell:
        """
        samtools view -F 4 {input.bam} {params.contig_id} -@ {threads}| cut -f1 | sort -u > {output.read_names}
        """

# =============================================================================
# RULE 3: EXTRACT AND COMBINE FASTQ SEQUENCES FOR SAMPLE1 READS
# =============================================================================

rule extract_fastq_sample1:
    """
    Extract read sequences from host-depleted FASTQ files for sample1.
    Combines R1 and R2 into a single FASTQ for simpler comparison.
    """
    input:
        read_names = "04b_qc/contig_reads/{pair_id}.sample1.read_names.txt",
        fq1 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}_1.fastq",
        fq2 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}_2.fastq"
    output:
        combined = temp("04b_qc/contig_reads/{pair_id}.sample1.combined.fastq.gz")
    threads:
        4
    shell:
        """
        # Extract reads from both R1 and R2, combine them
        seqkit grep -f {input.read_names} {input.fq1} > temp_{wildcards.pair_id}_s1_r1.fastq
        seqkit grep -f {input.read_names} {input.fq2} > temp_{wildcards.pair_id}_s1_r2.fastq
        cat temp_{wildcards.pair_id}_s1_r1.fastq temp_{wildcards.pair_id}_s1_r2.fastq | gzip > {output.combined}
        rm temp_{wildcards.pair_id}_s1_r1.fastq temp_{wildcards.pair_id}_s1_r2.fastq
        """

# =============================================================================
# RULE 4: EXTRACT AND COMBINE FASTQ SEQUENCES FOR SAMPLE2 READS
# =============================================================================

rule extract_fastq_sample2:
    """
    Extract read sequences from host-depleted FASTQ files for sample2.
    Combines R1 and R2 into a single FASTQ for simpler comparison.
    """
    input:
        read_names = "04b_qc/contig_reads/{pair_id}.sample2.read_names.txt",
        fq1 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}_1.fastq",
        fq2 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}_2.fastq"
    output:
        combined = temp("04b_qc/contig_reads/{pair_id}.sample2.combined.fastq.gz")
    threads:
        4
    shell:
        """
        # Extract reads from both R1 and R2, combine them
        seqkit grep -f {input.read_names} {input.fq1} > temp_{wildcards.pair_id}_s2_r1.fastq
        seqkit grep -f {input.read_names} {input.fq2} > temp_{wildcards.pair_id}_s2_r2.fastq
        cat temp_{wildcards.pair_id}_s2_r1.fastq temp_{wildcards.pair_id}_s2_r2.fastq | gzip > {output.combined}
        rm temp_{wildcards.pair_id}_s2_r1.fastq temp_{wildcards.pair_id}_s2_r2.fastq
        """

# =============================================================================
# RULE 5: COMPARE READ SEQUENCES BETWEEN SAMPLES
# =============================================================================

rule compare_reads:
    """
    Compare read sequences between the two samples to identify shared reads.
    Uses combined R1+R2 FASTQs for simpler, more accurate comparison.
    Also calculates contig lengths and median depth for edge case detection.
    """
    input:
        s1_combined = "04b_qc/contig_reads/{pair_id}.sample1.combined.fastq.gz",
        s2_combined = "04b_qc/contig_reads/{pair_id}.sample2.combined.fastq.gz",
        s1_fq1 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}_1.fastq",
        s1_fq2 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}_2.fastq",
        s2_fq1 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}_1.fastq",
        s2_fq2 = lambda wildcards: f"hostDepleted/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}_2.fastq",
        s1_depth = lambda wildcards: f"untrimmed_viral_contig_bowtie2/{PAIR_TO_INFO[wildcards.pair_id]['sample1']}.depth",
        s2_depth = lambda wildcards: f"untrimmed_viral_contig_bowtie2/{PAIR_TO_INFO[wildcards.pair_id]['sample2']}.depth"
    output:
        comparison = "04b_qc/read_comparisons/{pair_id}.comparison.tsv"
    params:
        pair_info = lambda wildcards: PAIR_TO_INFO[wildcards.pair_id]
    run:
        import gzip
        from itertools import islice
        import statistics

        def parse_fastq_sequences(filename):
            """
            Memory-efficient FASTQ parser.
            Yields only sequences (not read names).
            """
            if filename.endswith('.gz'):
                f = gzip.open(filename, 'rt')
            else:
                f = open(filename, 'r')

            try:
                while True:
                    record = list(islice(f, 4))
                    if not record:
                        break
                    if len(record) == 4:
                        seq = record[1].strip()
                        yield seq
            finally:
                f.close()

        def count_total_reads_in_fastq(fq1, fq2):
            """Count total reads in a pair of FASTQ files (R1 + R2)."""
            count = 0
            for _ in parse_fastq_sequences(fq1):
                count += 1
            for _ in parse_fastq_sequences(fq2):
                count += 1
            return count

        def get_contig_median_depth(depth_file, contig_id):
            """
            Calculate median depth for a specific contig from bedtools genomecov output.
            Depth file format: contig, start, end, depth (BED format from genomecov -bga)
            """
            depths = []
            with open(depth_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4 and parts[0] == contig_id:
                        start = int(parts[1])
                        end = int(parts[2])
                        depth = int(parts[3])
                        # Expand to per-base depth for accurate median
                        depths.extend([depth] * (end - start))

            if depths:
                return statistics.median(depths)
            return 0

        # Count total host-depleted reads for each sample
        sample1_total_host_depleted = count_total_reads_in_fastq(input.s1_fq1, input.s1_fq2)
        sample2_total_host_depleted = count_total_reads_in_fastq(input.s2_fq1, input.s2_fq2)

        # Collect sequences from contig-mapped reads (read once into list)
        sample1_seqs_list = list(parse_fastq_sequences(input.s1_combined))
        sample2_seqs_list = list(parse_fastq_sequences(input.s2_combined))

        # Convert to sets for fast lookup
        sample1_seqs_set = set(sample1_seqs_list)
        sample2_seqs_set = set(sample2_seqs_list)

        # Count contig-mapped reads
        sample1_contig_mapped = len(sample1_seqs_list)
        sample2_contig_mapped = len(sample2_seqs_list)

        # Find shared sequences
        shared_seqs = sample1_seqs_set & sample2_seqs_set

        # Count how many reads from each sample have shared sequences
        sample1_shared_count = sum(1 for seq in sample1_seqs_list if seq in shared_seqs)
        sample2_shared_count = sum(1 for seq in sample2_seqs_list if seq in shared_seqs)

        # Calculate percentages
        sample1_shared_pct = (sample1_shared_count / sample1_contig_mapped * 100) if sample1_contig_mapped > 0 else 0
        sample2_shared_pct = (sample2_shared_count / sample2_contig_mapped * 100) if sample2_contig_mapped > 0 else 0

        # Determine likely contributor and recipient
        # Higher shared % = more likely to be the recipient (contaminated sample)
        # Lower shared % = more likely to be the contributor (source of contamination)
        if sample2_shared_pct > sample1_shared_pct * 2:
            likely_recipient = params.pair_info['sample2']
            likely_contributor = params.pair_info['sample1']
            # Calculate unique reads in recipient
            recipient_unique_seqs = sample2_seqs_set - shared_seqs
            recipient_unique_distinct_seqs = len(recipient_unique_seqs)
            recipient_unique_total_reads = sample2_contig_mapped - sample2_shared_count
        elif sample1_shared_pct > sample2_shared_pct * 2:
            likely_recipient = params.pair_info['sample1']
            likely_contributor = params.pair_info['sample2']
            # Calculate unique reads in recipient
            recipient_unique_seqs = sample1_seqs_set - shared_seqs
            recipient_unique_distinct_seqs = len(recipient_unique_seqs)
            recipient_unique_total_reads = sample1_contig_mapped - sample1_shared_count
        else:
            likely_recipient = 'unclear'
            likely_contributor = 'unclear'
            recipient_unique_distinct_seqs = 'NA'
            recipient_unique_total_reads = 'NA'

        # Get contig lengths from pair_info (from 04a report)
        sample1_contig_len = params.pair_info['contig1_len']
        sample2_contig_len = params.pair_info['contig2_len']

        # Get alignment coordinates from pair_info (from 04a report)
        contig1_aln_start = params.pair_info['contig1_aln_start']
        contig1_aln_end = params.pair_info['contig1_aln_end']
        contig2_aln_start = params.pair_info['contig2_aln_start']
        contig2_aln_end = params.pair_info['contig2_aln_end']

        # Calculate median depth for each contig
        sample1_median_depth = get_contig_median_depth(input.s1_depth, params.pair_info['contig1'])
        sample2_median_depth = get_contig_median_depth(input.s2_depth, params.pair_info['contig2'])

        # Write results with new column structure
        with open(output.comparison, 'w') as f:
            f.write("pair_id\tsample1_id\tsample1_contig\tsample2_id\tsample2_contig\tgenus\t"
                   "sample1_total_host_depleted_reads\tsample2_total_host_depleted_reads\t"
                   "sample1_contig_mapped_reads\tsample2_contig_mapped_reads\t"
                   "sample1_contig_len\tsample2_contig_len\t"
                   "contig1_aln_start\tcontig1_aln_end\tcontig2_aln_start\tcontig2_aln_end\t"
                   "sample1_median_depth\tsample2_median_depth\t"
                   "sample1_shared_reads\tsample1_shared_pct\t"
                   "sample2_shared_reads\tsample2_shared_pct\t"
                   "likely_contributor\tlikely_recipient\t"
                   "recipient_unique_distinct_seqs\trecipient_unique_total_reads\n")

            f.write(f"{params.pair_info['pair_id']}\t"
                   f"{params.pair_info['sample1']}\t{params.pair_info['contig1']}\t"
                   f"{params.pair_info['sample2']}\t{params.pair_info['contig2']}\t"
                   f"{params.pair_info['genus']}\t"
                   f"{sample1_total_host_depleted}\t{sample2_total_host_depleted}\t"
                   f"{sample1_contig_mapped}\t{sample2_contig_mapped}\t"
                   f"{sample1_contig_len}\t{sample2_contig_len}\t"
                   f"{contig1_aln_start}\t{contig1_aln_end}\t{contig2_aln_start}\t{contig2_aln_end}\t"
                   f"{sample1_median_depth}\t{sample2_median_depth}\t"
                   f"{sample1_shared_count}\t{sample1_shared_pct:.2f}\t"
                   f"{sample2_shared_count}\t{sample2_shared_pct:.2f}\t"
                   f"{likely_contributor}\t{likely_recipient}\t"
                   f"{recipient_unique_distinct_seqs}\t{recipient_unique_total_reads}\n")

# =============================================================================
# RULE 6: AGGREGATE READ VALIDATION RESULTS
# =============================================================================

rule aggregate_read_validation:
    """
    Combine all pairwise read comparisons into a final summary report.
    Sorts by contamination percentage (descending) for easier review.
    """
    input:
        comparisons = expand("04b_qc/read_comparisons/{pair_id}.comparison.tsv",
                           pair_id=PAIR_IDS)
    output:
        summary = "04b_qc/read_validation_report/{batch}.read_validation_summary.tsv"
    params:
        batch = BATCH
    run:
        all_results = []

        for comp_file in input.comparisons:
            if os.path.exists(comp_file) and os.path.getsize(comp_file) > 1:
                df = pd.read_csv(comp_file, sep='\t')
                all_results.append(df)

        if all_results:
            combined = pd.concat(all_results, ignore_index=True)

            # Add helper column for sorting: max shared percentage between samples
            combined['max_shared_pct'] = combined[['sample1_shared_pct', 'sample2_shared_pct']].max(axis=1)

            # Sort by maximum shared percentage (descending) - highest contamination first
            combined = combined.sort_values('max_shared_pct', ascending=False)

            # Write to output
            combined.to_csv(output.summary, sep='\t', index=False)
        else:
            # Create empty file with header if no comparisons
            pd.DataFrame(columns=[
                'pair_id', 'sample1_id', 'sample1_contig', 'sample2_id', 'sample2_contig', 'genus',
                'sample1_total_host_depleted_reads', 'sample2_total_host_depleted_reads',
                'sample1_contig_mapped_reads', 'sample2_contig_mapped_reads',
                'sample1_contig_len', 'sample2_contig_len',
                'contig1_aln_start', 'contig1_aln_end', 'contig2_aln_start', 'contig2_aln_end',
                'sample1_median_depth', 'sample2_median_depth',
                'sample1_shared_reads', 'sample1_shared_pct',
                'sample2_shared_reads', 'sample2_shared_pct',
                'likely_contributor', 'likely_recipient',
                'recipient_unique_distinct_seqs', 'recipient_unique_total_reads', 'max_shared_pct'
            ]).to_csv(output.summary, sep='\t', index=False)

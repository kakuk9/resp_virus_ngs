"""
04c_final_detection_report.smk

Final integrated detection report combining:
- Unique read sequence counting per sample-genus
- Contig metadata from viral_contig_summary
- Within-sample nucmer deduplication (overlapping contigs)
- Cross-sample contamination detection from 04a (contig-level) and 04b (read-level)
- Consensus availability from 03_generate_consensus
- Genome coverage estimates

Produces final sample-genus level detection report with contamination filtering.

Author: Kirsty K
Date: 2026-01-15
"""

import pandas as pd
import os

# =============================================================================
# CONFIGURATION
# =============================================================================

BATCH = config.get("batch", "batch")

# Input files
VIRAL_CONTIG_SUMMARY = f"viral_contig_summary/batch/{BATCH}.genus.viral_contig_summary_all_samples.tsv"
COVERAGE_SUMMARY = f"coverage_estimate_summary/{BATCH}.project_summary.tsv"
READ_VALIDATION_REPORT_04B = f"04b_qc/read_validation_report/{BATCH}.read_validation_summary.tsv"
CONSENSUS_DIR = "03_generate_consensus/consensus"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_sample_genus_pairs():
    """
    Get all sample-genus pairs from viral_contig_summary.
    Only includes pairs where the BAM file exists.
    Returns list of tuples: (sample_id, genus)
    """
    if not os.path.exists(VIRAL_CONTIG_SUMMARY):
        print(f"Warning: Viral contig summary not found: {VIRAL_CONTIG_SUMMARY}")
        return []

    df = pd.read_csv(VIRAL_CONTIG_SUMMARY, sep='\t')

    # Filter for samples with mapped reads
    df = df[df['mapped_read_count'] > 0]

    pairs = []
    skipped_samples = set()
    for _, row in df.iterrows():
        sample_id = row['library_id']
        genus = row['uniq_value']

        # Check if BAM file exists for this sample
        bam_file = f"untrimmed_viral_contig_bowtie2/{sample_id}.bam"
        if os.path.exists(bam_file):
            pairs.append((sample_id, genus))
        else:
            skipped_samples.add(sample_id)

    if skipped_samples:
        print(f"Warning: Skipping {len(skipped_samples)} samples with missing BAM files: {', '.join(sorted(skipped_samples))}")

    print(f"Found {len(pairs)} sample-genus pairs with viral detections (and BAM files)")
    return pairs

# Get sample-genus pairs
SAMPLE_GENUS_PAIRS = get_sample_genus_pairs()
SAMPLES = list(set([pair[0] for pair in SAMPLE_GENUS_PAIRS]))

# Create lookup for sample-genus pair info
SAMPLE_GENUS_PAIR_IDS = [f"{sample}_{genus}" for sample, genus in SAMPLE_GENUS_PAIRS]

# =============================================================================
# RULE ALL
# =============================================================================

rule all:
    input:
        f"04c_final_report/{BATCH}.final_detection_report.tsv"

# =============================================================================
# RULE 1: WITHIN-SAMPLE NUCMER DEDUPLICATION
# =============================================================================

rule within_sample_nucmer:
    """
    Run nucmer within each sample to identify overlapping/redundant contigs.
    Uses high identity threshold (99%) to catch near-identical overlaps.
    """
    input:
        fasta = "04a_qc/contig_fasta/{sample}.viral_contigs.fasta"
    output:
        delta = "04c_qc/within_sample_nucmer/{sample}.delta",
        coords = "04c_qc/within_sample_nucmer/{sample}.coords"
    params:
        prefix = "04c_qc/within_sample_nucmer/{sample}"
    threads:
        4
    shell:
        """
        # Run nucmer against itself
        nucmer --maxmatch --nosimplify -p {params.prefix} {input.fasta} {input.fasta}

        # Convert to coords format
        show-coords -r -c -l -T {output.delta} > {output.coords}
        """

# =============================================================================
# RULE 2: PARSE WITHIN-SAMPLE OVERLAPS
# =============================================================================

rule parse_within_sample_overlaps:
    """
    Parse nucmer coords to identify redundant contigs within the same sample.
    Conservative deduplication criteria:
    1. Remove contigs that are 100% contained within another contig
    2. Remove contigs whose unique (non-overlapping) portion is < 300bp
    3. Always keep at least the longest contig per genus (even if all are redundant)
    """
    input:
        coords = "04c_qc/within_sample_nucmer/{sample}.coords",
        fasta = "04a_qc/contig_fasta/{sample}.viral_contigs.fasta",
        viral_summary = VIRAL_CONTIG_SUMMARY
    output:
        overlap_report = "04c_qc/within_sample_overlaps/{sample}.overlap_report.tsv"
    params:
        min_idy = 99.0,
        min_aln_len = 100,
        min_unique_len = 300
    run:
        import pandas as pd
        from collections import defaultdict

        # Helper function to merge overlapping intervals
        def merge_intervals(intervals):
            """
            Merge overlapping intervals and return merged list.
            Input: list of (start, end) tuples
            Output: list of merged (start, end) tuples
            """
            if not intervals:
                return []

            # Sort by start position
            intervals = sorted(intervals)
            merged = [intervals[0]]

            for current in intervals[1:]:
                last = merged[-1]
                if current[0] <= last[1]:  # Overlapping
                    # Merge by extending the end position
                    merged[-1] = (last[0], max(last[1], current[1]))
                else:
                    # Non-overlapping, add as new interval
                    merged.append(current)

            return merged

        # Parse FASTA to get contig lengths
        contig_lengths = {}
        if os.path.exists(input.fasta):
            from Bio import SeqIO
            for record in SeqIO.parse(input.fasta, 'fasta'):
                contig_lengths[record.id] = len(record.seq)

        # Parse nucmer coords to build overlap information
        overlaps = []

        if os.path.exists(input.coords) and os.path.getsize(input.coords) > 0:
            # Read nucmer coords
            df = pd.read_csv(input.coords, sep='\t', skiprows=4)

            # Rename columns (nucmer coords format)
            df.columns = ['S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'LEN_R', 'LEN_Q', 'COV_R', 'COV_Q', 'REF', 'QUERY']

            for _, row in df.iterrows():
                contig_ref = row['REF']
                contig_query = row['QUERY']

                # Skip self-alignments
                if contig_ref == contig_query:
                    continue

                identity = row['IDY']
                aln_length = row['LEN1']

                # Only consider high-identity alignments
                if identity >= params.min_idy and aln_length >= params.min_aln_len:
                    overlaps.append({
                        'ref': contig_ref,
                        'query': contig_query,
                        'ref_start': int(row['S1']),
                        'ref_end': int(row['E1']),
                        'query_start': int(row['S2']),
                        'query_end': int(row['E2']),
                        'cov_ref': row['COV_R'],
                        'cov_query': row['COV_Q'],
                        'ref_len': int(row['LEN_R']),
                        'query_len': int(row['LEN_Q'])
                    })

        # Identify contigs to remove
        contigs_to_remove = set()

        # Criterion 1: Remove contigs that are 100% contained within another contig
        for overlap in overlaps:
            # Check if query is 100% contained in ref
            if overlap['cov_query'] >= 99.9:  # 100% coverage of query
                # Query is completely contained in ref
                contigs_to_remove.add(overlap['query'])
            # Check if ref is 100% contained in query
            elif overlap['cov_ref'] >= 99.9:  # 100% coverage of ref
                # Ref is completely contained in query
                contigs_to_remove.add(overlap['ref'])

        # Criterion 2: Calculate unique portions and remove contigs with < 300bp unique sequence
        # Build list of overlapped regions for each contig
        contig_overlaps = defaultdict(list)

        for overlap in overlaps:
            # Add overlapped region to ref contig (as tuple)
            contig_overlaps[overlap['ref']].append(
                (overlap['ref_start'], overlap['ref_end'] + 1)
            )
            # Add overlapped region to query contig (as tuple)
            contig_overlaps[overlap['query']].append(
                (overlap['query_start'], overlap['query_end'] + 1)
            )

        # Calculate unique length for each contig
        for contig_id, contig_len in contig_lengths.items():
            if contig_id in contigs_to_remove:
                continue  # Already marked for removal

            # Calculate total overlapped length
            overlapped_length = 0
            if contig_id in contig_overlaps:
                # Merge overlapping intervals and sum their lengths
                merged_intervals = merge_intervals(contig_overlaps[contig_id])
                overlapped_length = sum(end - start for start, end in merged_intervals)

            # Calculate unique (non-overlapping) length
            unique_length = contig_len - overlapped_length

            # Mark for removal if unique portion is too short
            if unique_length < params.min_unique_len:
                contigs_to_remove.add(contig_id)

        # Load viral summary to get genus information for each contig
        # Ensure at least one contig is kept per genus (the longest one)
        viral_summary = pd.read_csv(input.viral_summary, sep='\t')
        sample_contigs = viral_summary[viral_summary['library_id'] == wildcards.sample]

        # Group contigs by genus and check if all would be removed
        for genus in sample_contigs['uniq_value'].unique():
            genus_rows = sample_contigs[sample_contigs['uniq_value'] == genus]
            # contig_ids is a pipe-separated list, expand and add sample prefix
            genus_contig_ids = set()
            for _, row in genus_rows.iterrows():
                if pd.notna(row['contig_ids']) and row['contig_ids']:
                    for cid in str(row['contig_ids']).split('|'):
                        genus_contig_ids.add(f"{wildcards.sample}_{cid}")

            # Check if all contigs of this genus would be removed
            if genus_contig_ids.issubset(contigs_to_remove):
                # All contigs would be removed - keep the longest one
                longest_contig = None
                longest_length = 0
                for contig_id in genus_contig_ids:
                    contig_len = contig_lengths.get(contig_id, 0)
                    if contig_len > longest_length:
                        longest_length = contig_len
                        longest_contig = contig_id

                if longest_contig:
                    # Remove the longest contig from the removal set
                    contigs_to_remove.discard(longest_contig)

        # Create output report with redundant contigs
        results = []
        for contig in contigs_to_remove:
            results.append({
                'sample_id': wildcards.sample,
                'redundant_contig': contig,
                'contig_length': contig_lengths.get(contig, 0)
            })

        # Write overlap report
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(output.overlap_report, sep='\t', index=False)
        else:
            # Empty file with header
            pd.DataFrame(columns=['sample_id', 'redundant_contig', 'contig_length']).to_csv(
                output.overlap_report, sep='\t', index=False)

# =============================================================================
# RULE 3: COUNT UNIQUE READ SEQUENCES PER SAMPLE-GENUS
# =============================================================================

rule count_unique_reads_per_genus:
    """
    Count unique read sequences mapped to contigs for each sample-genus pair.
    Splits counts into:
    - Total unique reads (all contigs)
    - Reads mapped to clean contigs
    - Reads mapped to contaminated contigs
    """
    input:
        bam = "untrimmed_viral_contig_bowtie2/{sample}.bam",
        viral_summary = VIRAL_CONTIG_SUMMARY,
        read_validation_04b = READ_VALIDATION_REPORT_04B
    output:
        counts = "04c_qc/unique_read_counts/{sample}_{genus}.unique_read_count.tsv"
    params:
        sample = "{sample}",
        genus = "{genus}",
        contamination_threshold = config.get("contamination_threshold", 10.0)
    threads:
        4
    run:
        import subprocess
        import pandas as pd

        def count_unique_reads(bam_file, contig_list, threads):
            """Count unique read sequences mapped to given contigs."""
            if not contig_list:
                return 0
            contigs_str = ' '.join(contig_list)
            cmd = f"samtools view -F 4 -@ {threads} {bam_file} {contigs_str} | cut -f10 | sort -u | wc -l"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return int(result.stdout.strip()) if result.stdout.strip() else 0

        def count_mapped_reads(bam_file, contig_list, threads):
            """Count total mapped reads (non-unique) to given contigs."""
            if not contig_list:
                return 0
            contigs_str = ' '.join(contig_list)
            cmd = f"samtools view -F 4 -c -@ {threads} {bam_file} {contigs_str}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return int(result.stdout.strip()) if result.stdout.strip() else 0

        # Read viral contig summary to get contig IDs for this sample-genus
        df = pd.read_csv(input.viral_summary, sep='\t')

        # Filter for this specific sample and genus
        sample_genus_df = df[
            (df['library_id'] == params.sample) &
            (df['uniq_value'] == params.genus)
        ]

        # Initialize counts
        num_unique_reads = 0
        num_unique_reads_clean_contigs = 0
        num_unique_reads_contaminated_contigs = 0
        mapped_reads_clean_contigs = 0
        mapped_reads_contaminated_contigs = 0

        if not sample_genus_df.empty:
            # Get all contig IDs (with sample prefix)
            contig_ids_str = sample_genus_df['contig_ids'].iloc[0]
            all_contig_ids = contig_ids_str.split('|') if pd.notna(contig_ids_str) else []

            if all_contig_ids:
                # Count total unique reads from all contigs
                num_unique_reads = count_unique_reads(input.bam, all_contig_ids, threads)

                # Load 04b report to identify contaminated contigs
                reads_04b = pd.read_csv(input.read_validation_04b, sep='\t')

                # Find contaminated contigs for this sample-genus
                contaminated_contigs = set()
                if not reads_04b.empty:
                    recipient_rows = reads_04b[
                        (reads_04b['likely_recipient'] == params.sample) &
                        (reads_04b['genus'] == params.genus) &
                        (reads_04b['max_shared_pct'] >= params.contamination_threshold)
                    ]

                    for _, recip_row in recipient_rows.iterrows():
                        # Get the contig ID for this sample (without sample prefix to match all_contig_ids)
                        if recip_row['sample1_id'] == params.sample:
                            contaminated_contig = recip_row['sample1_contig']
                        else:
                            contaminated_contig = recip_row['sample2_contig']
                        contaminated_contigs.add(contaminated_contig)

                # Split contigs into clean vs contaminated
                clean_contigs = [c for c in all_contig_ids if c not in contaminated_contigs]
                contaminated_contig_list = [c for c in all_contig_ids if c in contaminated_contigs]

                # Count reads for each category
                if clean_contigs:
                    num_unique_reads_clean_contigs = count_unique_reads(input.bam, clean_contigs, threads)
                    mapped_reads_clean_contigs = count_mapped_reads(input.bam, clean_contigs, threads)
                if contaminated_contig_list:
                    num_unique_reads_contaminated_contigs = count_unique_reads(input.bam, contaminated_contig_list, threads)
                    mapped_reads_contaminated_contigs = count_mapped_reads(input.bam, contaminated_contig_list, threads)

        # Write TSV output with 5 columns
        with open(output.counts, 'w') as f:
            f.write("num_unique_reads\tnum_unique_reads_clean_contigs\tnum_unique_reads_contaminated_contigs\tmapped_reads_clean_contigs\tmapped_reads_contaminated_contigs\n")
            f.write(f"{num_unique_reads}\t{num_unique_reads_clean_contigs}\t{num_unique_reads_contaminated_contigs}\t{mapped_reads_clean_contigs}\t{mapped_reads_contaminated_contigs}\n")

# =============================================================================
# RULE 4: AGGREGATE ALL DATA FOR FINAL REPORT
# =============================================================================

rule aggregate_detection_data:
    """
    Aggregate all detection data including:
    - Unique read sequence counts
    - Viral contig metadata
    - Within-sample overlaps
    - Cross-sample contamination (04a + 04b)
    - Coverage estimates
    - Consensus availability
    """
    input:
        viral_summary = VIRAL_CONTIG_SUMMARY,
        coverage_summary = COVERAGE_SUMMARY,
        read_validation_04b = READ_VALIDATION_REPORT_04B,
        within_sample_overlaps = expand("04c_qc/within_sample_overlaps/{sample}.overlap_report.tsv",
                                       sample=SAMPLES),
        unique_read_counts = expand("04c_qc/unique_read_counts/{pair_id}.unique_read_count.tsv",
                                   pair_id=SAMPLE_GENUS_PAIR_IDS)
    output:
        aggregated = f"04c_qc/aggregated_data/{BATCH}.aggregated_detection_data.tsv"
    params:
        consensus_dir = CONSENSUS_DIR,
        sample_genus_pairs = SAMPLE_GENUS_PAIRS,
        contamination_threshold = config.get("contamination_threshold", 10.0)
    run:
        import pandas as pd
        import os
        from collections import defaultdict

        # Load viral contig summary
        viral_df = pd.read_csv(input.viral_summary, sep='\t')
        viral_df = viral_df[viral_df['mapped_read_count'] > 0]  # Only samples with detections

        # Load coverage summary
        coverage_df = pd.read_csv(input.coverage_summary, sep='\t')

        # Load 04b read validation (read-level contamination)
        reads_04b = pd.read_csv(input.read_validation_04b, sep='\t')

        # Load all within-sample overlaps
        overlap_dfs = []
        for overlap_file in input.within_sample_overlaps:
            if os.path.exists(overlap_file) and os.path.getsize(overlap_file) > 1:
                df = pd.read_csv(overlap_file, sep='\t')
                if not df.empty:
                    overlap_dfs.append(df)

        within_sample_overlaps = pd.concat(overlap_dfs, ignore_index=True) if overlap_dfs else pd.DataFrame()

        # Load unique read counts (now TSV with 4 columns)
        unique_read_counts = {}
        for pair_file in input.unique_read_counts:
            # Extract sample_genus from filename
            basename = os.path.basename(pair_file)
            pair_id = basename.replace('.unique_read_count.tsv', '')

            counts_df = pd.read_csv(pair_file, sep='\t')
            if not counts_df.empty:
                # Handle both old (4 column) and new (5 column) formats
                unique_read_counts[pair_id] = {
                    'num_unique_reads': int(counts_df['num_unique_reads'].iloc[0]),
                    'num_unique_reads_clean_contigs': int(counts_df['num_unique_reads_clean_contigs'].iloc[0]),
                    'num_unique_reads_contaminated_contigs': int(counts_df['num_unique_reads_contaminated_contigs'].iloc[0]),
                    'mapped_reads_clean_contigs': int(counts_df['mapped_reads_clean_contigs'].iloc[0]),
                    'mapped_reads_contaminated_contigs': int(counts_df['mapped_reads_contaminated_contigs'].iloc[0]) if 'mapped_reads_contaminated_contigs' in counts_df.columns else 0
                }
            else:
                unique_read_counts[pair_id] = {
                    'num_unique_reads': 0,
                    'num_unique_reads_clean_contigs': 0,
                    'num_unique_reads_contaminated_contigs': 0,
                    'mapped_reads_clean_contigs': 0,
                    'mapped_reads_contaminated_contigs': 0
                }

        # Build sample-genus level aggregated data
        results = []

        for _, row in viral_df.iterrows():
            sample_id = row['library_id']
            genus = row['uniq_value']
            segment = row['segment_name']

            # Get unique read counts for this sample-genus
            pair_id = f"{sample_id}_{genus}"
            read_counts = unique_read_counts.get(pair_id, {
                'num_unique_reads': 0,
                'num_unique_reads_clean_contigs': 0,
                'num_unique_reads_contaminated_contigs': 0,
                'mapped_reads_clean_contigs': 0,
                'mapped_reads_contaminated_contigs': 0
            })
            num_unique_reads = read_counts['num_unique_reads']
            num_unique_reads_clean_contigs = read_counts['num_unique_reads_clean_contigs']
            num_unique_reads_contaminated_contigs = read_counts['num_unique_reads_contaminated_contigs']
            mapped_reads_clean_contigs = read_counts['mapped_reads_clean_contigs']
            mapped_reads_contaminated_contigs = read_counts['mapped_reads_contaminated_contigs']

            # Parse contig IDs (pipe-delimited)
            contig_ids_str = row['contig_ids']
            contig_ids = contig_ids_str.split('|') if pd.notna(contig_ids_str) else []

            # Remove sample prefix from contig IDs
            contig_ids_clean = [cid.replace(f"{sample_id}_", "") for cid in contig_ids]

            # Count initial contigs
            n_contigs_initial = len(contig_ids_clean)

            # Start with all contigs for this genus
            contigs_remaining = set(contig_ids_clean)

            # Filter 1: Remove redundant/overlapping contigs
            redundant_contigs = set()
            if not within_sample_overlaps.empty:
                sample_overlaps = within_sample_overlaps[within_sample_overlaps['sample_id'] == sample_id]
                for _, overlap_row in sample_overlaps.iterrows():
                    # Get redundant contig name
                    redundant_contig = overlap_row['redundant_contig']

                    # Remove sample prefix from contig name
                    redundant_contig_clean = redundant_contig.replace(f"{sample_id}_", "")

                    # Only count if this contig belongs to this genus
                    if redundant_contig_clean in contig_ids_clean:
                        redundant_contigs.add(redundant_contig_clean)

            # Remove redundant contigs from remaining set
            contigs_remaining -= redundant_contigs
            n_contigs_non_redundant = len(contigs_remaining)

            # Filter 2: Remove contaminated contigs
            contaminated_contigs = set()
            if not reads_04b.empty:
                # Check if this sample is a likely recipient with max_shared_pct >= threshold
                recipient_rows = reads_04b[
                    (reads_04b['likely_recipient'] == sample_id) &
                    (reads_04b['genus'] == genus) &
                    (reads_04b['max_shared_pct'] >= params.contamination_threshold)
                ]

                for _, recip_row in recipient_rows.iterrows():
                    # Determine which contig column corresponds to this sample
                    if recip_row['sample1_id'] == sample_id:
                        contaminated_contig = recip_row['sample1_contig']
                    else:
                        contaminated_contig = recip_row['sample2_contig']

                    contaminated_contigs.add(contaminated_contig)

            # Remove contaminated contigs from remaining set
            contigs_remaining -= contaminated_contigs
            n_contigs_contaminated = len(contaminated_contigs)

            # Final count after all filters
            n_contigs_after_filter = len(contigs_remaining)

            # Get coverage information for this sample-genus-segment
            coverage_pct = 0.0
            cov_row = coverage_df[coverage_df['sample_id'] == sample_id]
            if not cov_row.empty:
                high_cov_str = cov_row.iloc[0].get('high_coverage_genera_segments', 'None')
                if pd.notna(high_cov_str) and high_cov_str != 'None':
                    # Parse format: "genus_segment:coverage|genus_segment:coverage"
                    segments_list = high_cov_str.split('|')
                    for seg_info in segments_list:
                        if ':' in seg_info:
                            seg_name, seg_cov = seg_info.split(':')
                            if seg_name == f"{genus}_{segment}":
                                coverage_pct = float(seg_cov) * 100  # Convert to percentage
                                break

            # Check consensus availability
            consensus_file = os.path.join(params.consensus_dir, f"{sample_id}--{genus}_{segment}.consensus.fasta")
            consensus_available = os.path.exists(consensus_file)

            # Determine contamination status
            is_recipient = sample_id in reads_04b[reads_04b['max_shared_pct'] >= params.contamination_threshold]['likely_recipient'].values
            contributor_samples = []
            if is_recipient:
                contrib_rows = reads_04b[
                    (reads_04b['likely_recipient'] == sample_id) &
                    (reads_04b['genus'] == genus) &
                    (reads_04b['max_shared_pct'] >= params.contamination_threshold)
                ]
                contributor_samples = contrib_rows['likely_contributor'].unique().tolist()

            contamination_summary = 'Clean'
            if contributor_samples:
                contamination_summary = f"Likely recipient from {','.join(contributor_samples)}"

            # Aggregate results
            results.append({
                'sample_id': sample_id,
                'genus': genus,
                'segment': segment,
                'segment_type': 'non-segmented' if segment == 'non-segmented' else 'segmented',
                'num_unique_reads': num_unique_reads,
                'num_unique_reads_clean_contigs': num_unique_reads_clean_contigs,
                'num_unique_reads_contaminated_contigs': num_unique_reads_contaminated_contigs,
                'mapped_reads_genus_total': row['mapped_read_count'],
                'mapped_reads_clean_contigs': mapped_reads_clean_contigs,
                'mapped_reads_contaminated_contigs': mapped_reads_contaminated_contigs,
                'n_contigs_initial': n_contigs_initial,
                'n_contigs_non_redundant': n_contigs_non_redundant,
                'n_contigs_contaminated': n_contigs_contaminated,
                'n_contigs_after_filter': n_contigs_after_filter,
                'longest_contig_length': row['longest_contig_len'],
                'genome_coverage_pct': coverage_pct,
                'consensus_available': consensus_available,
                'contamination_summary': contamination_summary
            })

        # Write aggregated data
        results_df = pd.DataFrame(results)
        results_df.to_csv(output.aggregated, sep='\t', index=False)

# =============================================================================
# RULE 5: COLLAPSE TO SAMPLE-GENUS LEVEL
# =============================================================================

rule collapse_to_sample_genus:
    """
    Collapse segment-level data to sample-genus level.
    For segmented viruses: aggregate across all segments.
    For non-segmented viruses: one row per sample-genus.
    """
    input:
        aggregated = f"04c_qc/aggregated_data/{BATCH}.aggregated_detection_data.tsv"
    output:
        final_report = f"04c_final_report/{BATCH}.final_detection_report.tsv"
    params:
        min_unique_reads = 2  # minimum unique reads for positive detection
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.aggregated, sep='\t')

        # Group by sample_id and genus
        grouped = df.groupby(['sample_id', 'genus'])

        final_results = []

        for (sample_id, genus), group in grouped:
            # Determine if segmented or non-segmented
            segment_type = group['segment_type'].iloc[0]

            # Get all segments for this genus
            segments_detected = group['segment'].tolist()
            n_segments_detected = len(segments_detected)

            # Expected segments (need to infer from data or use reference)
            # For now, use detected count as proxy
            n_segments_expected = n_segments_detected if segment_type == 'segmented' else 1

            # Get read counts (should be same across all segments for this sample-genus)
            num_unique_reads = group['num_unique_reads'].iloc[0]
            num_unique_reads_clean_contigs = group['num_unique_reads_clean_contigs'].iloc[0]
            num_unique_reads_contaminated_contigs = group['num_unique_reads_contaminated_contigs'].iloc[0]
            mapped_reads_genus_total = group['mapped_reads_genus_total'].sum()
            mapped_reads_clean_contigs = group['mapped_reads_clean_contigs'].iloc[0]
            mapped_reads_contaminated_contigs = group['mapped_reads_contaminated_contigs'].iloc[0]

            # Aggregate metrics
            n_contigs_initial = group['n_contigs_initial'].sum()
            n_contigs_non_redundant = group['n_contigs_non_redundant'].sum()
            n_contigs_contaminated = group['n_contigs_contaminated'].sum()
            n_contigs_after_filter = group['n_contigs_after_filter'].sum()
            longest_contig_length = group['longest_contig_length'].max()

            # Coverage aggregation
            coverage_list = []
            for _, row in group.iterrows():
                seg = row['segment']
                cov = row['genome_coverage_pct']
                if cov > 0:
                    coverage_list.append(f"{seg}:{cov:.1f}%")

            genome_coverage_pct = '|'.join(coverage_list) if coverage_list else 'NA'

            # Consensus availability
            consensus_count = group['consensus_available'].sum()
            consensus_segments_list = group[group['consensus_available'] == True]['segment'].tolist()

            consensus_segments_available = f"{consensus_count} segment{'s' if consensus_count != 1 else ''}"
            consensus_segments_list_str = '|'.join(consensus_segments_list) if consensus_segments_list else ''

            # Segments detected list
            segments_detected_str = '|'.join(segments_detected)

            # Contamination status (take first non-"Clean" summary, or "Clean" if all clean)
            contamination_summaries = group['contamination_summary'].unique()
            contamination_status = 'Clean'
            for summary in contamination_summaries:
                if summary != 'Clean':
                    contamination_status = summary
                    break

            # Detection verdict: Positive if:
            # 1. >= 1 contig after filter AND
            # 2. >= 2 unique reads AND
            # 3. reads from clean contigs > reads from contaminated contigs
            if (n_contigs_after_filter >= 1 and
                num_unique_reads_clean_contigs >= params.min_unique_reads and
                num_unique_reads_clean_contigs > num_unique_reads_contaminated_contigs):
                detection_verdict = 'Positive'
            else:
                detection_verdict = 'Negative'

            final_results.append({
                'sample_id': sample_id,
                'genus': genus,
                'segment_type': segment_type,
                'n_segments_expected': n_segments_expected,
                'segments_detected': segments_detected_str,
                'n_segments_detected': n_segments_detected,
                'genome_coverage_pct': genome_coverage_pct,
                'consensus_segments_available': consensus_segments_available,
                'consensus_segments_list': consensus_segments_list_str,
                'num_unique_reads': num_unique_reads,
                'num_unique_reads_clean_contigs': num_unique_reads_clean_contigs,
                'num_unique_reads_contaminated_contigs': num_unique_reads_contaminated_contigs,
                'mapped_reads_genus_total': mapped_reads_genus_total,
                'mapped_reads_clean_contigs': mapped_reads_clean_contigs,
                'mapped_reads_contaminated_contigs': mapped_reads_contaminated_contigs,
                'n_contigs_initial': n_contigs_initial,
                'n_contigs_non_redundant': n_contigs_non_redundant,
                'n_contigs_contaminated': n_contigs_contaminated,
                'n_contigs_after_filter': n_contigs_after_filter,
                'longest_contig_length': longest_contig_length,
                'detection_verdict': detection_verdict,
                'contamination_status': contamination_status
            })

        # Create final report
        final_df = pd.DataFrame(final_results)

        # Sort by sample_id and genus
        final_df = final_df.sort_values(['sample_id', 'genus'])

        # Write final report
        final_df.to_csv(output.final_report, sep='\t', index=False)

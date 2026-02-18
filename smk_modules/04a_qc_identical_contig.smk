"""
04a_qc_identical_contig.smk

Contig-based similarity detection using nucmer for viral mNGS data.
Compares assembled contigs within the same virus genus to identify:
- Highly similar contigs between samples (potential contamination)
- Strain variants and evolutionary relationships
- Conserved genomic regions

Author: Kirsty K
Date: 2026-01-08
"""

import pandas as pd
import os
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

BATCH = config.get("batch", "batch")

# Nucmer similarity thresholds
MIN_IDENTITY = config.get("min_contig_identity", 95.0)
MIN_ALIGNMENT_LENGTH = config.get("min_alignment_length", 100)
MIN_QUERY_COVERAGE = config.get("min_query_coverage", 80.0)


def get_viral_samples():
    """Get list of samples with viral contigs."""
    summary_file = f"viral_contig_summary/batch/{BATCH}.genus.viral_contig_summary_all_samples.tsv"
    df = pd.read_csv(summary_file, sep="\t")
    viral_samples = df[df['mapped_read_count'] > 0]['library_id'].unique().tolist()
    return viral_samples

VIRAL_SAMPLES = get_viral_samples()
print(f"Viral samples identified for contig similarity analysis: {len(VIRAL_SAMPLES)}; {VIRAL_SAMPLES}")

rule all:
    input:
        expand("04a_qc/nucmer_output/{batch}.identical_contig_report.tsv", batch=BATCH),
        expand("04a_qc/nucmer_output/{batch}.nucmer_output.coords.parsed", batch=BATCH)

rule gather_all_viral_contigs:
    """
    Extract list of viral contig IDs per sample from family summary.
    """
    input:
        viral_contig_summary = "viral_contig_summary/family/{sample}.viral_contig_summary"
    output:
        contig_list = temp("04a_qc/contig_lists/{sample}.contig_list.txt")
    run:
        df = pd.read_csv(input.viral_contig_summary, sep='\t')

        # Collect all contig IDs
        all_contigs = []
        for _, row in df.iterrows():
            if pd.notna(row['contig_ids']) and row['contig_ids'] != 'No_contig_found':
                # Split pipe-delimited contig IDs
                contig_ids = str(row['contig_ids']).split('|')
                all_contigs.extend([c.strip() for c in contig_ids if c.strip()])

        # Write unique contig IDs
        with open(output.contig_list, 'w') as f:
            for contig in sorted(set(all_contigs)):
                f.write(f"{contig}\n")

rule gather_all_viral_contigs_fasta:
    """
    Gather all viral contig sequences into a single FASTA file per sample.
    """
    input:
        contig_list = "04a_qc/contig_lists/{sample}.contig_list.txt",
        contig_fasta = "viral_contig/{sample}.fasta"
    output: 
        fasta_tmp = temp("04a_qc/contig_fasta/{sample}.viral_contigs.fasta.tmp"),
        fasta = "04a_qc/contig_fasta/{sample}.viral_contigs.fasta"
    params:
        sample = "{sample}"
    shell:
        """
        seqkit grep -f {input.contig_list} {input.contig_fasta} > {output.fasta_tmp}
        seqkit replace -p "^" -r "{params.sample}_" {output.fasta_tmp} > {output.fasta}
        """
rule nucmer_contig_similarity:
    """
    Perform pairwise contig similarity analysis using nucmer.
    """
    input:
        expand("04a_qc/contig_fasta/{sample}.viral_contigs.fasta", sample=VIRAL_SAMPLES)
    output:
        master_fasta = temp("04a_qc/nucmer_output/{batch}_combined_contigs.fasta"),
        delta = "04a_qc/nucmer_output/{batch}.nucmer_output",
        coords = "04a_qc/nucmer_output/{batch}.nucmer_output.coords"
    params:
        batch = BATCH
    shell:
        """
        cat {input} > {output.master_fasta}
        nucmer -p {params.batch} {output.master_fasta} {output.master_fasta} --maxmatch
        mv {params.batch}.delta {output.delta}
        show-coords -rcl {output.delta} > {output.coords}
        """
rule identify_identical_contigs_between_samples:
    """
    """
    input:
        coords = "04a_qc/nucmer_output/{batch}.nucmer_output.coords",
        summary = "viral_contig_summary/batch/{batch}.genus.viral_contig_summary_all_samples.tsv"
    output:
        parsed_coords = "04a_qc/nucmer_output/{batch}.nucmer_output.coords.parsed",
        identical_contig_report = "04a_qc/nucmer_output/{batch}.identical_contig_report.tsv"
    params:
        idy_threshold = 100,
        contig_cov_threshold = 100
    run:
        # Load taxonomy info from summary file
        summary_df = pd.read_csv(input.summary, sep='\t')

        # Create lookup dictionary: (sample_id, contig_id) -> (uniq_value, segment_name)
        contig_taxonomy_map = {}
        for _, row in summary_df.iterrows():
            sample_id = row['library_id']
            if pd.notna(row['contig_ids']) and row['contig_ids'] != 'No_contig_found':
                uniq_value = row.get('uniq_value', 'Unknown')
                segment_name = row.get('segment_name', 'Unknown')
                contig_ids = str(row['contig_ids']).split('|')
                for contig_id in contig_ids:
                    contig_id = contig_id.strip()
                    # Map using (sample_id, contig_id) as key
                    contig_taxonomy_map[(sample_id, contig_id)] = (uniq_value, segment_name)

        keep = []

        with open(input.coords, 'r') as f:
            coords = [x.strip().replace("|", "").split() for x in f.readlines()[5:]]
            # Check for 100% identity
            for coord in coords:
                contig1_start, contig1_end, contig2_start, contig2_end, contig1_len, contig2_len, pid, aln_contig1_len, aln_contig2_len, contig1_cov, contig2_cov, contig1_id, contig2_id = coord
                sample1_id = "_".join(contig1_id.split('_')[0:3])
                sample2_id = "_".join(contig2_id.split('_')[0:3])

                # Skip if not 100% identity
                if float(pid) < params.idy_threshold:
                    continue
                # Skip if same contig or same sample
                elif contig1_id == contig2_id or sample1_id == sample2_id:
                    continue
                # Require at least one contig to have 100% coverage (full contig match)
                elif float(contig1_cov) < params.contig_cov_threshold and float(contig2_cov) < params.contig_cov_threshold:
                    continue
                else:
                    # Clean contig IDs by removing sample prefix for lookup
                    contig1_id_clean = contig1_id.replace(f"{sample1_id}_", "")
                    contig2_id_clean = contig2_id.replace(f"{sample2_id}_", "")

                    # Get taxonomy for each contig using (sample_id, clean_contig_id) tuple
                    taxonomy1 = contig_taxonomy_map.get((sample1_id, contig1_id_clean), ('Unknown', 'Unknown'))
                    taxonomy2 = contig_taxonomy_map.get((sample2_id, contig2_id_clean), ('Unknown', 'Unknown'))

                    uniq_value1, segment_name1 = taxonomy1
                    uniq_value2, segment_name2 = taxonomy2

                    coord.extend([sample1_id, sample2_id, uniq_value1, segment_name1, uniq_value2, segment_name2])
                    keep.append(coord)

        header = ["contig1_start", "contig1_end", "contig2_start", "contig2_end", "contig1_len", "contig2_len",
        "percent_identity", "aln_contig1_len", "aln_contig2_len", "contig1_cov", "contig2_cov",
        "contig1_id", "contig2_id", "sample1_id", "sample2_id", "sample1_genus", "sample1_segment_name",
        "sample2_genus", "sample2_segment_name"]

        # Write parsed coords
        with open(output.parsed_coords, 'w') as out_f:
            out_f.write("\t".join(header) + "\n")
            for item in coords:
                out_f.write("\t".join(item) + "\n")

        # Write identical contig report with taxonomy
        header_with_taxonomy = header
        with open(output.identical_contig_report, 'w') as out_f:
            out_f.write("\t".join(header_with_taxonomy) + "\n")
            for item in keep:
                out_f.write("\t".join(item) + "\n")


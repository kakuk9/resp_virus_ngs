"""
Estimate genome coverage for each genus in each sample (segmented virus compatible)
Author: Kirsty K
"""

import os
import glob
import pandas as pd
from snakemake.io import Wildcards

# --- Configuration and Setup ---
INPUT_DIR = config["run_directory"]
GENUS_MAX_SEGMENTS_FILE = config["genus_max_segments"]

# Load genus max segments lookup table
GENUS_MAX_SEGMENTS = {}
if os.path.exists(GENUS_MAX_SEGMENTS_FILE):
    genus_seg_df = pd.read_csv(GENUS_MAX_SEGMENTS_FILE, sep='\t')
    GENUS_MAX_SEGMENTS = dict(zip(genus_seg_df['genus'], genus_seg_df['max_segments']))

# --- Helper Functions ---

def get_samples_from_blast(run_directory):
    """Dynamically finds sample names based on a configurable run directory."""
    search_path = os.path.join(run_directory, "merged_blastn_refseq", "*.merged_blastn_refseq")
    SAMPLES = [os.path.basename(f).replace(".merged_blastn_refseq", "") for f in glob.glob(search_path)]
    return sorted(list(set(SAMPLES)))

SAMPLES = get_samples_from_blast(INPUT_DIR)
print(f"Assessment workflow identified samples: {SAMPLES}")

def get_coverage_targets_for_sample(wildcards):
    """
    Reads the checkpoint metadata for a SINGLE sample and returns a list
    of all coverage reports that need to be generated for it.
    """
    metadata_file = checkpoints.determine_best_ref.get(sample=wildcards.sample).output.best_ref_tsv
    
    if not os.path.exists(metadata_file) or os.path.getsize(metadata_file) == 0:
        return []
        
    df = pd.read_csv(metadata_file, sep='\t')
    df = df[df.get('status') != 'all_contigs_discarded'].copy()

    coverage_reports = []
    for fasta_path in df['fasta_file'].dropna():
        basename = os.path.basename(fasta_path).replace('.fasta', '')
        if '--' in basename:
            virus_segment = basename.split('--', 1)[1]
            coverage_reports.append(f"coverage_report/{wildcards.sample}/{virus_segment}.txt")
    
    return coverage_reports

def get_all_summary_inputs(wildcards):
    """Gathers all inputs needed for the final summarize_assessment rule."""
    all_metadata_tsvs = expand("best_reference/{sample}_best_references_metadata.tsv", sample=SAMPLES)
    
    all_coverage_reports = []
    for sample_id in SAMPLES:
        # Create a temporary wildcards object to pass to the helper function
        w = Wildcards(fromdict={"sample": sample_id})
        all_coverage_reports.extend(get_coverage_targets_for_sample(w))

    # Return a single, flat list of all required files
    return all_metadata_tsvs + all_coverage_reports


# --- Main Workflow ---

rule all:
    """The final goal is the single project-wide summary file."""
    input:
        f"coverage_estimate_summary/{config['batch']}.project_summary.tsv"

# --- STEP 1: Determine Best References (Checkpoint) ---
checkpoint determine_best_ref:
    input:
        merged_blastn_refseq="merged_blastn_refseq/{sample}.merged_blastn_refseq",
        trimming_summary="trimmed_viral_contig/{sample}.summary"
    output:
        best_ref_tsv="best_reference/{sample}_best_references_metadata.tsv"
    params:
        sample_prefix="{sample}",
        scripts_dir=config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        fasta_output_dir="best_reference"
    shell:
        """
        python {params.scripts_dir}/best_ref_mergedblastn.py {input.merged_blastn_refseq} {input.trimming_summary} \
            -o {output.best_ref_tsv} \
            --sample_prefix {params.sample_prefix} \
            --fasta_output_dir {params.fasta_output_dir}
        """


# --- STEP 2: Map Reads to Each Best Reference ---
rule map_to_best_reference:
    input:
        fasta="best_reference/{sample}--{virus_segment}.fasta",
        fq1="hostDepleted/{sample}_1.fastq",
        fq2="hostDepleted/{sample}_2.fastq"
    output:
        bam="map_to_reference/{sample}/{virus_segment}.bam"
    threads: 8
    shell:
        """
        bowtie2-build {input.fasta} {input.fasta} > /dev/null 2>&1
        bowtie2 -q -x {input.fasta} -1 {input.fq1} -2 {input.fq2} --threads {threads} | \
        samtools view -@ {threads} -bS | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# --- STEP 3: Calculate Coverage for Each Mapping ---
rule calculate_coverage:
    input:
        bam="map_to_reference/{sample}/{virus_segment}.bam"
    output:
        coverage_report="coverage_report/{sample}/{virus_segment}.txt"
    shell:
        "bedtools genomecov -bga -ibam {input.bam} > {output.coverage_report}"

# --- STEP 4: Create the Final Summary Table ---
rule summarize_assessment:
    input:
        get_all_summary_inputs
    output:
        project_summary=f"coverage_estimate_summary/{config['batch']}.project_summary.tsv"
    params:
        coverage_cutoff=config.get("coverage_cutoff", 0.70),
        min_depth=config.get("min_depth", 5)
    run:
        # Separate the flat input list back into their respective groups
        all_metadata_tsvs = [f for f in input if f.endswith(".tsv")]
        all_coverage_reports = [f for f in input if f.endswith(".txt")]

        all_sample_results = []

        for sample_id in SAMPLES:
            metadata_file = f"best_reference/{sample_id}_best_references_metadata.tsv"
            total_genera_count = 0
            virus_segment_to_genus = {}
            all_detected_segments = []
            discarded_segments = set()
            found_segments_by_genus = {}

            if metadata_file in all_metadata_tsvs and os.path.exists(metadata_file) and os.path.getsize(metadata_file) > 0:
                meta_df = pd.read_csv(metadata_file, sep='\t')
                total_genera_count = meta_df['genus_b'].nunique()

                for _, row in meta_df.iterrows():
                    genus = row['genus_b']
                    segment = row['segment_key']
                    genus_segment = f"{genus}_{segment}"
                    all_detected_segments.append(genus_segment)

                    # Track found segments per genus
                    if genus not in found_segments_by_genus:
                        found_segments_by_genus[genus] = set()
                    found_segments_by_genus[genus].add(segment)

                    if row.get('status') == 'all_contigs_discarded':
                        discarded_segments.add(genus_segment)
                    elif 'fasta_file' in row.index:
                        basename = os.path.basename(row['fasta_file']).replace('.fasta', '')
                        if '--' in basename:
                            virus_segment = basename.split('--', 1)[1]
                            virus_segment_to_genus[virus_segment] = genus_segment

            high_cov_genera = set()
            high_cov_segments = []
            low_cov_segments = []
            coverage_data = {}

            sample_reports = [r for r in all_coverage_reports if f"/{sample_id}/" in r]

            for report_file in sample_reports:
                virus_segment = os.path.basename(report_file).replace('.txt', '')
                genus_segment = virus_segment_to_genus.get(virus_segment)
                if not genus_segment: continue

                df = pd.read_csv(report_file, sep='\t', header=None, names=['ref', 'start', 'end', 'depth'])
                if df.empty: continue

                genome_length = df['end'].max()
                bases_covered = df[df['depth'] >= params.min_depth]['end'].sum() - df[df['depth'] >= params.min_depth]['start'].sum()

                if genome_length > 0:
                    coverage_fraction = bases_covered / genome_length
                    coverage_data[genus_segment] = coverage_fraction

                    if coverage_fraction >= params.coverage_cutoff:
                        genus = genus_segment.rsplit('_', 1)[0]
                        high_cov_genera.add(genus)
                        high_cov_segments.append(f"{genus_segment}:{coverage_fraction:.2f}")
                    else:
                        low_cov_segments.append(f"{genus_segment}:{coverage_fraction:.2f}")

            # For discarded segments, add "NA" to low coverage
            for discarded in discarded_segments:
                if discarded not in coverage_data:
                    low_cov_segments.append(f"{discarded}:NA")

            # Identify missing segments for known segmented viruses
            missing_segments_info = []
            for genus, found_segs in found_segments_by_genus.items():
                if genus in GENUS_MAX_SEGMENTS:
                    expected_count = GENUS_MAX_SEGMENTS[genus]
                    found_count = len(found_segs)
                    if found_count < expected_count:
                        # Determine which specific segments are missing
                        expected_segs = set([f"segment-{i}" for i in range(1, expected_count + 1)])
                        # Also check for alternative naming patterns like "segment RNA 1"
                        if 'non-segmented' in found_segs:
                            continue  # Skip non-segmented
                        missing = expected_segs - found_segs
                        if missing:
                            missing_segments_info.append(f"{genus}:found_{found_count}/{expected_count}")

            all_sample_results.append({
                "sample_id": sample_id,
                "total_genera_count": total_genera_count,
                "high_coverage_genera_count": len(high_cov_genera),
                "high_coverage_genera_list": "|".join(sorted(list(high_cov_genera))) or "None",
                "all_detected_genera_segments": "|".join(sorted(all_detected_segments)) or "None",
                "high_coverage_genera_segments": "|".join(sorted(high_cov_segments)) or "None",
                "low_coverage_genera_segments": "|".join(sorted(low_cov_segments)) or "None",
                "missing_expected_segments": "|".join(sorted(missing_segments_info)) or "None"
            })

        final_df = pd.DataFrame(all_sample_results)
        final_df.to_csv(output.project_summary, sep='\t', index=False)
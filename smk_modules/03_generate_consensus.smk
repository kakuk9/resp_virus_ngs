"""
Snakemake workflow for generating consensus sequences from mNGS data using RefStitch.
RefStitch-only version (Gapfiller removed based on validation showing RefStitch superiority).
Last updated: 2025-10-24
Author: Kirsty K
"""

import os
import pandas as pd
import pyfastx

# --- Configuration and Setup ---
# Use the config object provided by Snakemake
RUN_DIR = config["run_directory"]
BATCH = config["batch"]
SCRIPT_DIR = config.get("scripts_dir", os.path.join(os.path.dirname(workflow.snakefile), "scripts"))

# --- Helper Functions ---
def get_samples_from_summary(run_directory, batch):
    """
    Reads the project_summary.tsv to get a list of (sample, virus_segment) tuples
    for samples with high coverage genera segments.
    """
    sample_list = []
    summary_file = os.path.join(run_directory, "coverage_estimate_summary", f"{batch}.project_summary.tsv")

    if not os.path.exists(summary_file):
        print(f"Warning: Summary file not found at {summary_file}")
        return sample_list

    df = pd.read_csv(summary_file, sep='\t')

    for _, row in df.iterrows():
        sample_id = row['sample_id']
        high_cov = row.get('high_coverage_genera_segments', 'None')

        # Skip if no high coverage segments
        if high_cov == 'None' or pd.isna(high_cov):
            continue

        # Parse the high_coverage_genera_segments column
        # Format: genus_segment:coverage|genus_segment:coverage
        segments = high_cov.split('|')
        for segment_info in segments:
            if ':' in segment_info:
                virus_segment = segment_info.split(':')[0]  # Extract genus_segment part
                sample_list.append((sample_id, virus_segment))

    if not sample_list:
        print("Warning: No valid (sample, virus_segment) combinations with high coverage were found.")

    return sample_list

# Generate the list of combinations. If it's empty, zip will raise an error, so we handle it.
VALID_COMBINATIONS = get_samples_from_summary(RUN_DIR, BATCH)
if VALID_COMBINATIONS:
    SAMPLES, VIRUS_SEGMENTS = zip(*VALID_COMBINATIONS)
else:
    SAMPLES, VIRUS_SEGMENTS = [], []
print(f"Consensus workflow identified {len(VALID_COMBINATIONS)} sample-virus combinations:")
print(f"Samples: {SAMPLES}")
print(f"Virus segments: {VIRUS_SEGMENTS}")

# --- Snakemake Rules ---

rule all:
    input:
        "03_generate_consensus/consensus/{batch}_consensus_stats_annotated.tsv".format(batch=BATCH),
        "03_generate_consensus/checkv/quality_summary.tsv".format(batch=BATCH),
        "03_generate_consensus/consensus_unscaffolded/{batch}_unscaffolded_stats.tsv".format(batch=BATCH)

rule group_contig:
    input:
        viral_contig_summary = "viral_contig_summary/genus/{sample}.viral_contig_summary",
        fasta = "trimmed_viral_contig/{sample}.fasta"
    output:
        contig_list = temp("03_generate_consensus/trimmed_viral_contig/{sample}--{virus_segment}.txt"),
        fasta = "03_generate_consensus/trimmed_viral_contig/by_genus_segment/{sample}--{virus_segment}.fasta"
    params:
        genus = lambda wildcards: wildcards.virus_segment.rsplit('_', 1)[0],
        segment = lambda wildcards: wildcards.virus_segment.rsplit('_', 1)[1]
    shell:
        """
        # Extract contig IDs for this virus_segment from metadata
        csvtk -t filter2 -f '$uniq_value=="{params.genus}" && $segment_name=="{params.segment}"' {input.viral_contig_summary} | \
        csvtk -t cut -f contig_ids | \
        csvtk -t del-header | \
        tr '|' '\n' > {output.contig_list}

        # Extract those contigs from the trimmed fasta
        seqkit grep -nrif {output.contig_list} {input.fasta} -o {output.fasta}
        """

rule nucmer_scaffold:
    input:
        contig="03_generate_consensus/trimmed_viral_contig/by_genus_segment/{sample}--{virus_segment}.fasta",
        ref="best_reference/{sample}--{virus_segment}.fasta"
    output:
        scaffold="03_generate_consensus/scaffolds/{sample}--{virus_segment}_final_scaffold.fasta",
        tiling="03_generate_consensus/scaffolds/{sample}--{virus_segment}_tiling.tsv",
        summary="03_generate_consensus/scaffolds/{sample}--{virus_segment}_summary.txt"
    params:
        scripts_dir=config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        coverage_threshold = config.get("scaffold_coverage_threshold", 80)
    shell:
        """
        # Run scaffold_nucmer_coords.py with explicit output paths and coverage threshold
        python3 {params.scripts_dir}/scaffold_nucmer_coords.py {input.ref} {input.contig} \
            --fill r \
            --coverage_threshold {params.coverage_threshold} \
            --output_scaffold {output.scaffold} \
            --output_tiling {output.tiling} \
            --output_summary {output.summary}
        """

rule bowtie2_alignment:
    input:
        reference = "03_generate_consensus/scaffolds/{sample}--{virus_segment}_final_scaffold.fasta",
        fastq = os.path.join(RUN_DIR, "hostDepleted/{sample}_1.fastq"),
        fastq2 = os.path.join(RUN_DIR, "hostDepleted/{sample}_2.fastq"),
        singleton = os.path.join(RUN_DIR, "hostDepleted/{sample}_singletons.fastq")
    output:
        bam = "03_generate_consensus/alignment/{sample}--{virus_segment}.bam",
        bam_index = "03_generate_consensus/alignment/{sample}--{virus_segment}.bam.bai"
    params:
        samtools = config.get("samtools", "/software/samtools-v1.16.1/bin/samtools")
    shell:
        """
        bowtie2-build {input.reference} {input.reference}
        bowtie2 -x {input.reference} -1 {input.fastq} -2 {input.fastq2} -U {input.singleton} | {params.samtools} view -bS | {params.samtools} sort -o {output.bam}
        {params.samtools} index {output.bam}
        """

rule mpileup:
    input:
        bam = "03_generate_consensus/alignment/{sample}--{virus_segment}.bam",
        fasta = "03_generate_consensus/scaffolds/{sample}--{virus_segment}_final_scaffold.fasta"
    output:
        pileup = "03_generate_consensus/alignment/{sample}--{virus_segment}.pileup"
    params:
        samtools = config.get("samtools", "/software/samtools-v1.16.1/bin/samtools")
    shell:
        """
        {params.samtools} mpileup -aa -A -d 0 -Q 0 -f {input.fasta} {input.bam} > {output.pileup}
        """

rule mpileup_summary:
    input:
        pileup = "03_generate_consensus/alignment/{sample}--{virus_segment}.pileup"
    output:
        pileup_summary = "03_generate_consensus/alignment/{sample}--{virus_segment}.pileup_summary.txt"
    params:
        script = os.path.join(SCRIPT_DIR, "pileup_summary.py")
    shell:
        """
        python {params.script} {input.pileup} -out {output.pileup_summary}
        """

rule consensus:
    input:
        fasta = "03_generate_consensus/scaffolds/{sample}--{virus_segment}_final_scaffold.fasta",
        pileup_summary = "03_generate_consensus/alignment/{sample}--{virus_segment}.pileup_summary.txt"
    output:
        consensus = temp("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.raw.fasta"),
        mapping = "03_generate_consensus/consensus/{sample}--{virus_segment}.coordinate_map.txt"
    params:
        script = os.path.join(SCRIPT_DIR, "pileup2consensus.py"),
        min_depth = config.get("min_depth", 10)
    shell:
        """
        python {params.script} {input.fasta} {input.pileup_summary} -d {params.min_depth} -o {output.consensus} -m {output.mapping}
        """

rule trim_consensus_terminal_n:
    input:
        consensus = "03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.raw.fasta"
    output:
        trimmed = temp("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.tmp.fasta")
    params:
        script = os.path.join(SCRIPT_DIR, "trim_terminal_n.py")
    shell:
        """
        python {params.script} -i {input.consensus} -o {output.trimmed}
        """

rule rename_consensus_header:
    input:
        consensus = "03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.tmp.fasta"
    output:
        renamed = temp("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.renamed.fasta")
    params:
        new_header = lambda wildcards: f"{wildcards.sample}--{wildcards.virus_segment.replace('_non-segmented', '')}"
    shell:
        """
        # Replace ALL sequence headers with sample--virus_segment (without _non-segmented)
        # This handles both 'final_scaffold' and any other headers like 'NODE_*'
        seqkit replace -p '.+' -r '{params.new_header}' {input.consensus} -o {output.renamed}
        """

rule route_consensus:
    input:
        consensus = "03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.renamed.fasta",
        summary = "03_generate_consensus/scaffolds/{sample}--{virus_segment}_summary.txt"
    output:
        scaffolded = "03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.fasta",
        unscaffolded = "03_generate_consensus/consensus_unscaffolded/{sample}--{virus_segment}.consensus.fasta"
    run:
        import os

        # Read the scaffold summary to check if scaffolding succeeded
        with open(input.summary, 'r') as f:
            summary_text = f.read()

        scaffolding_failed = "Scaffolding failed" in summary_text

        if scaffolding_failed:
            # Scaffolding failed - keep only first sequence (which is the longest from the original contigs)
            # Create empty file for scaffolded output (required by Snakemake)
            shell("touch {output.scaffolded}")
            # Keep first sequence only for unscaffolded (seqkit head doesn't fail on duplicate IDs)
            shell("seqkit head -n 1 {input.consensus} > {output.unscaffolded}")
        else:
            # Scaffolding succeeded - should have single sequence, put in main consensus folder
            # Create empty file for unscaffolded output
            shell("touch {output.unscaffolded}")
            # Copy to scaffolded folder (should already be single sequence)
            shell("cp {input.consensus} {output.scaffolded}")
        

rule consensus_stats:
    input:
        consensus_files=expand("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS),
        reference_files=expand("best_reference/{sample}--{virus_segment}.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS)
    output:
        "03_generate_consensus/consensus/{batch}_consensus_stats.tsv".format(batch=BATCH)
    run:
        import os
        records = []

        # Create a mapping of sample--virus_segment to reference file
        ref_map = {}
        for ref_file in input.reference_files:
            ref_prefix = os.path.basename(ref_file).replace(".fasta", "")
            ref_map[ref_prefix] = ref_file

        for f in input.consensus_files:
            # Skip empty placeholder files (from failed scaffolding)
            if os.path.getsize(f) == 0:
                continue
            prefix = os.path.basename(f).replace(".consensus.fasta", "")
            parts = prefix.split("--")
            sample_id = parts[0]
            virus_segment = parts[1] if len(parts) > 1 else "unknown"

            # Read consensus FASTA (use simple file reading to avoid memory issues)
            consensus_seq = ""
            try:
                with open(f, 'r') as fh:
                    for line in fh:
                        if not line.startswith('>'):
                            consensus_seq += line.strip()
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")
                continue

            # Calculate consensus stats
            consensus_len = len(consensus_seq)
            n_count = consensus_seq.upper().count('N')
            called_bases = consensus_len - n_count
            completeness = (called_bases / consensus_len * 100) if consensus_len > 0 else 0

            # Read reference FASTA to get reference length
            ref_file = ref_map.get(prefix)
            ref_length = 0
            if ref_file and os.path.exists(ref_file):
                try:
                    ref_seq = ""
                    with open(ref_file, 'r') as fh:
                        for line in fh:
                            if not line.startswith('>'):
                                ref_seq += line.strip()
                    ref_length = len(ref_seq)
                except Exception as e:
                    print(f"Warning: Could not read reference {ref_file}: {e}")
                    ref_length = 0

            # Calculate genome coverage (consensus length / reference length)
            genome_coverage = (consensus_len / ref_length * 100) if ref_length > 0 else 0

            records.append({
                'sample_id': sample_id,
                'virus_segment': virus_segment,
                'reference_length': ref_length,
                'consensus_length': consensus_len,
                'called_bases': called_bases,
                'N_count': n_count,
                'completeness_pct': round(completeness, 2),
                'genome_coverage_pct': round(genome_coverage, 2),
                'consensus_file': f
            })

        import pandas as pd
        df = pd.DataFrame(records)
        df.to_csv(output[0], sep='\t', index=False)

rule checkv_qc:
    """Combine successfully scaffolded consensus sequences and run CheckV end_to_end quality control."""
    input:
        consensus_files=expand("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS)
    output:
        combined = temp("03_generate_consensus/checkv/{batch}_all_consensus.fasta".format(batch=BATCH)),
        quality_summary = "03_generate_consensus/checkv/quality_summary.tsv",
        completeness = "03_generate_consensus/checkv/completeness.tsv",
        contamination = "03_generate_consensus/checkv/contamination.tsv",
        complete_genomes = "03_generate_consensus/checkv/complete_genomes.tsv"
    params:
        outdir = "03_generate_consensus/checkv",
        checkv_db = config.get("checkv_db", "checkv-db-v1.5"),
        batch = BATCH
    threads: 8
    run:
        import os
        # Filter out empty placeholder files (failed scaffolds)
        non_empty_files = [f for f in input.consensus_files if os.path.getsize(f) > 0]

        if non_empty_files:
            # Combine non-empty consensus sequences
            shell("cat {non_empty_files} > {output.combined}")
            # Remove old CheckV output if exists
            shell("find {params.outdir} -type f ! -name '{params.batch}_all_consensus.fasta' -delete 2>/dev/null || true")
            # Run CheckV end_to_end
            shell("checkv end_to_end {output.combined} {params.outdir} -t {threads} -d {params.checkv_db}")
        else:
            # No valid consensus files - create empty outputs
            shell("touch {output.combined} {output.quality_summary} {output.completeness} {output.contamination} {output.complete_genomes}")

rule combine_consensus_fasta:
    """Combine successfully scaffolded consensus sequences into a single FASTA for BLASTN."""
    input:
        consensus_files=expand("03_generate_consensus/consensus/{sample}--{virus_segment}.consensus.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS)
    output:
        combined = temp("03_generate_consensus/consensus/{batch}_all_consensus.fasta".format(batch=BATCH))
    run:
        import os
        # Filter out empty placeholder files (failed scaffolds)
        non_empty_files = [f for f in input.consensus_files if os.path.getsize(f) > 0]

        if non_empty_files:
            shell("cat {non_empty_files} > {output.combined}")
        else:
            # Create empty file if no valid consensus
            shell("touch {output.combined}")

rule blastn_consensus:
    """BLASTN consensus sequences against nt database."""
    input:
        fasta = "03_generate_consensus/consensus/{batch}_all_consensus.fasta"
    output:
        blastn = temp("03_generate_consensus/consensus/{batch}_consensus.blastn.tmp")
    params:
        batch = BATCH
    threads: 8
    shell:
        r"""
        blastn -db nt -query {input.fasta} -out {output.blastn} -num_threads {threads} -max_target_seqs 10 -outfmt "6 qseqid sseqid sscinames staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"
        """

rule blastn_filter:
    """Filter BLASTN results to keep top hit only."""
    input:
        blastn = "03_generate_consensus/consensus/{batch}_consensus.blastn.tmp"
    output:
        filtered = "03_generate_consensus/consensus/{batch}_consensus.blastn.filtered"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        batch = BATCH
    shell:
        """
        python {params.scripts_dir}/blastn_filter.py {input.blastn} {output.filtered} --colname qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore -n 1 --skip-exclusion
        """

rule annotate_consensus_stats:
    """Left join BLASTN results with consensus stats."""
    input:
        stats = "03_generate_consensus/consensus/{batch}_consensus_stats.tsv",
        blastn = "03_generate_consensus/consensus/{batch}_consensus.blastn.filtered"
    output:
        annotated = "03_generate_consensus/consensus/{batch}_consensus_stats_annotated.tsv"
    run:
        import pandas as pd

        # Read the consensus stats
        stats_df = pd.read_csv(input.stats, sep='\t')

        # Read the BLASTN results
        blastn_df = pd.read_csv(input.blastn, sep='\t')

        # Create a join key: sample_id--virus_segment (remove _non-segmented suffix to match BLASTN qseqid)
        stats_df['join_key'] = stats_df['sample_id'] + '--' + stats_df['virus_segment'].str.replace('_non-segmented', '')

        # Rename qseqid to join_key for merging
        blastn_df = blastn_df.rename(columns={'qseqid': 'join_key'})

        # Left join: keep all stats, add BLASTN info where available
        merged_df = pd.merge(stats_df, blastn_df, on='join_key', how='left')

        # Drop the join_key column
        merged_df = merged_df.drop(columns=['join_key'])

        # Reorder columns to have BLASTN info after the basic stats
        # Keep original columns first, then add BLASTN columns
        original_cols = list(stats_df.columns)
        original_cols.remove('join_key')
        blastn_cols = [col for col in merged_df.columns if col not in original_cols]
        final_cols = original_cols + blastn_cols
        merged_df = merged_df[final_cols]

        # Save the annotated stats
        merged_df.to_csv(output.annotated, sep='\t', index=False)

rule unscaffolded_consensus_stats:
    """Generate stats for unscaffolded (failed scaffolding) consensus sequences."""
    input:
        consensus_files=expand("03_generate_consensus/consensus_unscaffolded/{sample}--{virus_segment}.consensus.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS),
        reference_files=expand("best_reference/{sample}--{virus_segment}.fasta",
               zip,
               sample=SAMPLES,
               virus_segment=VIRUS_SEGMENTS)
    output:
        "03_generate_consensus/consensus_unscaffolded/{batch}_unscaffolded_stats.tsv".format(batch=BATCH)
    run:
        import os
        records = []

        # Create a mapping of sample--virus_segment to reference file
        ref_map = {}
        for ref_file in input.reference_files:
            ref_prefix = os.path.basename(ref_file).replace(".fasta", "")
            ref_map[ref_prefix] = ref_file

        for f in input.consensus_files:
            # Skip empty placeholder files (successful scaffolds)
            if os.path.getsize(f) == 0:
                continue

            prefix = os.path.basename(f).replace(".consensus.fasta", "")
            parts = prefix.split("--")
            sample_id = parts[0]
            virus_segment = parts[1] if len(parts) > 1 else "unknown"

            # Read consensus FASTA
            consensus_seq = ""
            try:
                with open(f, 'r') as fh:
                    for line in fh:
                        if not line.startswith('>'):
                            consensus_seq += line.strip()
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")
                continue

            # Calculate consensus stats
            consensus_len = len(consensus_seq)
            n_count = consensus_seq.upper().count('N')
            called_bases = consensus_len - n_count
            completeness = (called_bases / consensus_len * 100) if consensus_len > 0 else 0

            # Read reference FASTA to get reference length
            ref_file = ref_map.get(prefix)
            ref_length = 0
            if ref_file and os.path.exists(ref_file):
                try:
                    ref_seq = ""
                    with open(ref_file, 'r') as fh:
                        for line in fh:
                            if not line.startswith('>'):
                                ref_seq += line.strip()
                    ref_length = len(ref_seq)
                except Exception as e:
                    print(f"Warning: Could not read reference {ref_file}: {e}")
                    ref_length = 0

            # Calculate genome coverage
            genome_coverage = (consensus_len / ref_length * 100) if ref_length > 0 else 0

            records.append({
                'sample_id': sample_id,
                'virus_segment': virus_segment,
                'reference_length': ref_length,
                'consensus_length': consensus_len,
                'called_bases': called_bases,
                'N_count': n_count,
                'completeness_pct': round(completeness, 2),
                'genome_coverage_pct': round(genome_coverage, 2),
                'consensus_file': f,
                'scaffold_status': 'FAILED'
            })

        import pandas as pd
        df = pd.DataFrame(records)
        df.to_csv(output[0], sep='\t', index=False)

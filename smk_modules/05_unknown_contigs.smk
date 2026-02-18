"""
Snakemake workflow for identifying contigs without Diamond BLASTX hits.
This workflow cross-checks contigs in filteredContigs against m8 files to identify
contigs that were assembled but have no taxonomic assignment.

Author: Kirsty K
"""

import os
import pandas as pd

# --- Configuration and Setup ---
RUN_DIR = config["run_directory"]
BATCH = config["batch"]
SCRIPT_DIR = config.get("scripts_dir", os.path.join(os.path.dirname(workflow.snakefile), "scripts"))

# --- Helper Functions ---
def get_samples_from_directory(run_directory):
    """
    Get sample list from filteredContigs directory.
    """
    import glob
    samples = []
    pattern = os.path.join(run_directory, "filteredContigs", "*_filterContigs.fasta")

    for filepath in glob.glob(pattern):
        basename = os.path.basename(filepath)
        # Remove _filterContigs.fasta suffix to get sample name
        sample = basename.replace("_filterContigs.fasta", "")
        samples.append(sample)

    if not samples:
        print(f"Warning: No samples found in {run_directory}/filteredContigs/")

    return samples

SAMPLES = get_samples_from_directory(RUN_DIR)
print(f"Unknown contigs workflow identified {len(SAMPLES)} samples:")
print(f"Samples: {SAMPLES}")

# --- Snakemake Rules ---

rule all:
    input:
        expand("04_unknown_contigs_output/unknown_to_diamond_and_blastn/{sample}.true_unknown_contigs.fasta", sample=SAMPLES)
rule extract_contigs_with_hits:
    input:
        m8 = os.path.join(RUN_DIR, "diamond/{sample}.m8")
    output:
        contig_list = temp("04_unknown_contigs_output/temp/{sample}_contigs_with_hits.txt")
    shell:
        """
        # Extract unique contig IDs from m8 file (column 1)
        # If m8 file is empty, create an empty output file
        if [ -s {input.m8} ]; then
            cut -f1 {input.m8} | sort -u > {output.contig_list}
        else
            touch {output.contig_list}
        fi
        """

rule extract_all_contigs:
    input:
        fasta = os.path.join(RUN_DIR, "filteredContigs/{sample}_filterContigs.fasta")
    output:
        contig_list = temp("04_unknown_contigs_output/temp/{sample}_all_contigs.txt")
    shell:
        """
        # Extract all contig headers from FASTA
        grep "^>" {input.fasta} | sed 's/^>//' | cut -d' ' -f1 > {output.contig_list}
        """

rule identify_unknown_contigs:
    input:
        all_contigs = "04_unknown_contigs_output/temp/{sample}_all_contigs.txt",
        contigs_with_hits = "04_unknown_contigs_output/temp/{sample}_contigs_with_hits.txt"
    output:
        unknown_list = "04_unknown_contigs_output/unknown_to_diamond/{sample}_unknown_contigs.txt"
    shell:
        """
        # Find contigs that are in all_contigs but NOT in contigs_with_hits
        # Use comm to find lines unique to file 1 (all_contigs)
        comm -23 <(sort {input.all_contigs}) <(sort {input.contigs_with_hits}) > {output.unknown_list}
        """

rule extract_unknown_contig_fasta:
    input:
        fasta = os.path.join(RUN_DIR, "filteredContigs/{sample}_filterContigs.fasta"),
        unknown_list = "04_unknown_contigs_output/unknown_to_diamond/{sample}_unknown_contigs.txt"
    output:
        fasta = "04_unknown_contigs_output/unknown_to_diamond/{sample}.unknown_contigs_to_diamond.fasta"
    shell:
        """
        # Extract unknown contigs from the original FASTA file
        # If no unknown contigs, create empty file
        if [ -s {input.unknown_list} ]; then
            seqkit grep -nrf {input.unknown_list} {input.fasta} | seqkit sort -l -r -o {output.fasta}
        else
            touch {output.fasta}
        fi
        """

rule unknown_contigs_blastn:
    input:
        unknown_contigs = "04_unknown_contigs_output/unknown_to_diamond/{sample}.unknown_contigs_to_diamond.fasta"
    output:
        blastn_raw = "04_unknown_contigs_output/unknown_to_diamond/blastn/{sample}.unknown2diamond.blast.raw",
        blastn_raw2 = "04_unknown_contigs_output/unknown_to_diamond/blastn/{sample}.unknown2diamond.blast.raw2",
        blastn_filtered = "04_unknown_contigs_output/unknown_to_diamond/blastn/{sample}.unknown2diamond.blastnannotated"
    params:
        scripts_dir=config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        prefix="{sample}"
    threads: 8
    shell:
        r"""
        blastn -db nt -query {input.unknown_contigs} -out {output.blastn_raw} -num_threads {threads} -max_target_seqs 10 -outfmt "6 qseqid sseqid sscinames staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"
        python {params.scripts_dir}/blastn_filter.py {output.blastn_raw} {output.blastn_raw2} --colname qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore -n 5 --skip-exclusion --skip-filtering
        bash {params.scripts_dir}/add_taxon.sh -i {output.blastn_raw2} -f blastn -k staxid -e -c 4 -n qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,library_id,accNum
        cat {params.prefix}.blastnannotated > {output.blastn_filtered}
        rm {params.prefix}.blastnannotated
        """

rule unknown_contigs_to_blastn:
    input:
        blastn_filtered = "04_unknown_contigs_output/unknown_to_diamond/blastn/{sample}.unknown2diamond.blastnannotated",
        fasta = "04_unknown_contigs_output/unknown_to_diamond/{sample}.unknown_contigs_to_diamond.fasta"
    output:
        contigs_with_blastn_hits = temp("04_unknown_contigs_output/unknown_to_diamond/{sample}_contigs_with_blastn_hits.tmp"),
        fasta = "04_unknown_contigs_output/unknown_to_diamond_and_blastn/{sample}.true_unknown_contigs.fasta"
    shell:
        """
        # Extract contigs that have blastn hits
        csvtk cut -t -f qseqid {input.blastn_filtered} | csvtk uniq | csvtk del-header > {output.contigs_with_blastn_hits}

        # Filter out contigs with blastn hits, keeping only true unknowns
        if [ -s {output.contigs_with_blastn_hits} ]; then
            seqkit grep -v -nrf {output.contigs_with_blastn_hits} {input.fasta} | seqkit sort -l -r -o {output.fasta}
        else
            # If no blastn hits found, all contigs are true unknowns
            seqkit sort -l -r {input.fasta} -o {output.fasta}
        fi
        """
        



# rule summarize_unknown_contigs:
#     input:
#         unknown_lists = expand("04_unknown_contigs_output/{sample}_unknown_contigs.txt", sample=SAMPLES),
#         unknown_fastas = expand("04_unknown_contigs_output/{sample}.unknown_contigs_to_diamond.fasta", sample=SAMPLES),
#         all_contig_lists = expand("04_unknown_contigs_output/temp/{sample}_all_contigs.txt", sample=SAMPLES)
#     output:
#         summary = "04_unknown_contigs_output/{batch}_unknown_contigs_summary.tsv"
#     run:
#         records = []

#         for sample in SAMPLES:
#             # Count total contigs
#             all_contig_file = f"04_unknown_contigs_output/temp/{sample}_all_contigs.txt"
#             with open(all_contig_file, 'r') as f:
#                 total_contigs = sum(1 for line in f if line.strip())

#             # Count unknown contigs
#             unknown_file = f"04_unknown_contigs_output/{sample}_unknown_contigs.txt"
#             with open(unknown_file, 'r') as f:
#                 unknown_count = sum(1 for line in f if line.strip())

#             # Calculate stats
#             contigs_with_hits = total_contigs - unknown_count
#             unknown_pct = (unknown_count / total_contigs * 100) if total_contigs > 0 else 0

#             # Get unknown contig fasta path
#             unknown_fasta = f"04_unknown_contigs_output/{sample}.unknown_contigs_to_diamond.fasta"

#             records.append({
#                 'sample_id': sample,
#                 'total_contigs': total_contigs,
#                 'contigs_with_hits': contigs_with_hits,
#                 'unknown_contigs': unknown_count,
#                 'unknown_percentage': round(unknown_pct, 2),
#                 'unknown_contigs_fasta': unknown_fasta
#             })

#         df = pd.DataFrame(records)
#         df.to_csv(output.summary, sep='\t', index=False)

#         print("\n=== Unknown Contigs Summary ===")
#         print(df.to_string(index=False))

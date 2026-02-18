"""
Snakemake workflow for respiratory metagenomic NGS analysis.
Complete pipeline from raw reads to viral contig summaries.
Last updated on 18/12/2025 by Kirsty K.

Pipeline steps:
1. Adapter trimming (TrimGalore)
2. Host depletion (Bowtie2 + Samtools)
3. De novo assembly (metaSPAdes)
4. Contig filtering (SeqKit)
5. Protein BLAST (DIAMOND)
6. Taxonomic annotation (TaxonKit)
7. Viral triage and phage removal
8. BLASTN against nt and RefSeq
9. Contig depth calculation and trimming
10. Summary statistics generation

Usage:
    snakemake -s 01_reads_to_viral_contigs.smk --configfile <config.yaml> --cores <N>
"""


# --- Configuration and Setup ---
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc;")

import os
import re
import csv
import glob
import gzip
import pandas as pd
import subprocess

# Pre-DIAMOND configuration (for assembly and alignment)
RAW_FASTQ_DIR = config['rawFastqDir']
HOST_REF = config['hostRef']
DIAMOND_db = config['DiamondDB']
num_underscore = config['numUnderscore']
Illumina_sequencer = config['IlluminaSequencer']
samtools = config['samtools']

# Post-DIAMOND configuration (for annotation and summary)
VIRUS_TAXON_SEGMENT = pd.read_csv(config["virus_segment_info"], sep="\t")
VIRUS_TAXON_SEGMENT = dict(zip(VIRUS_TAXON_SEGMENT['taxon_name'], VIRUS_TAXON_SEGMENT['segmented']))

# Sequencer-specific FASTQ suffixes
fastq_suffix = ['1', '2', 'singletons']

if Illumina_sequencer.lower() == "nextseq":
    RawFastq_suffix = "/*_R1_001.fastq.gz"
    RFastq1_suffix = "_R1_001.fastq.gz"
    RFastq2_suffix = "_R2_001.fastq.gz"
elif Illumina_sequencer.lower() == "miseq":
    RawFastq_suffix = "/*_L001_R1_001.fastq.gz"
    RFastq1_suffix = "_L001_R1_001.fastq.gz"
    RFastq2_suffix = "_L001_R2_001.fastq.gz"
else:
    print("Please specify the sequencer used for the run in the config file. Options are MiSeq or NextSeq.")
    exit(1)


# --- Helper Functions ---
def list_samples(RawFastqDirectory, RawFastq_suffix):
    """Extract sample names from raw FASTQ files."""
    SAMPLES = []
    PATH = RawFastqDirectory + RawFastq_suffix
    print(f"PATH {PATH}")
    for FILE in glob.glob(PATH):
        base = os.path.basename(FILE)
        # Take up to N underscores as the sample name e.g. NT_4_SS_S1 -> NT_4_SS
        sample = "_".join(base.split("_")[:num_underscore])
        # Don't analyse the undetermined samples
        if "Undetermined" not in sample:
            SAMPLES.append(sample)
    return sorted(list(set(SAMPLES)))

def virus_taxon_segment_info(SPECIES, GENUS, FAMILY):
    """Look up whether a virus is segmented based on species, genus, or family."""
    if SPECIES in VIRUS_TAXON_SEGMENT:
        return VIRUS_TAXON_SEGMENT[SPECIES]
    elif GENUS in VIRUS_TAXON_SEGMENT:
        return VIRUS_TAXON_SEGMENT[GENUS]
    elif FAMILY in VIRUS_TAXON_SEGMENT:
        return VIRUS_TAXON_SEGMENT[FAMILY]
    else:
        return "Unknown"

# Get sample list from raw FASTQ files
SAMPLES = list_samples(RAW_FASTQ_DIR, RawFastq_suffix)
print(SAMPLES)


# --- Master rule ---
rule all:
    input:
        expand("viral_contig_summary/batch/{batch}.family.viral_contig_summary_all_samples.tsv", batch = config.get('batch', 'batch')),
        expand("viral_contig_summary/batch/{batch}.genus.viral_contig_summary_all_samples.tsv", batch = config.get('batch', 'batch')),
        expand("m8_2_triage/{sample}.m8bacteria", sample=SAMPLES),
        expand("m8_2_triage/{sample}.m8eukaryota", sample=SAMPLES)

# =============================================================================
# SECTION 1: ADAPTER TRIMMING AND HOST DEPLETION
# =============================================================================

rule trimGalore_adapterTrimming:
    """Remove adapters and low-quality bases using TrimGalore."""
    input:
        fq = f"{RAW_FASTQ_DIR}/{{sample}}{RFastq1_suffix}",
        fq2 = f"{RAW_FASTQ_DIR}/{{sample}}{RFastq2_suffix}"
    output:
        fq = "trimGalore/{sample}_trimmed_1.fq",
        fq2 = "trimGalore/{sample}_trimmed_2.fq",
        unpaired = "trimGalore/{sample}_trimmed_singletons.fq"
    threads: 8
    shell:
        """
        trim_galore --length 50 -q 30 --stringency 1 --dont_gzip --paired {input.fq} {input.fq2} --retain_unpaired -o trimGalore -j {threads}
        mv trimGalore/{wildcards.sample}*_val_1.fq {output.fq}
        mv trimGalore/{wildcards.sample}*_val_2.fq {output.fq2}
        cat trimGalore/{wildcards.sample}*_unpaired_1.fq trimGalore/{wildcards.sample}*_unpaired_2.fq > {output.unpaired}
        """

rule bowtie2_4hostDepletion:
    """Align reads to host reference genome."""
    input:
        fq = rules.trimGalore_adapterTrimming.output.fq,
        fq2 = rules.trimGalore_adapterTrimming.output.fq2,
        unpaired = rules.trimGalore_adapterTrimming.output.unpaired
    output:
        bam = "bowtie2_4hostDepletion/{sample}_host.bam"
    params:
        ref_fasta = HOST_REF,
        samtools = samtools
    threads: 8
    shell:
        """
        bowtie2 -q -x {params.ref_fasta} -1 {input.fq} -2 {input.fq2} -U {input.unpaired} --threads {threads} | {params.samtools} view -bS | {params.samtools} sort -o {output.bam}
        {params.samtools} index {output.bam}
        """

rule samtools_hostDepletion:
    """Extract unmapped (non-host) reads."""
    input:
        bam = rules.bowtie2_4hostDepletion.output.bam
    output:
        fq = "hostDepleted/{sample}_1.fastq",
        fq2 = "hostDepleted/{sample}_2.fastq",
        singletons = "hostDepleted/{sample}_singletons.fastq"
    params:
        samtools = samtools
    threads: 8
    shell:
        """
        {params.samtools} view -b -f4 {input.bam} | {params.samtools} fastq -@ {threads} -1 {output.fq} -2 {output.fq2} -s {output.singletons} -0 /dev/null -n
        """

# =============================================================================
# SECTION 2: DE NOVO ASSEMBLY AND FILTERING
# =============================================================================

rule metaspades_deNovoAssembly:
    """Assemble host-depleted reads using metaSPAdes."""
    input:
        fq = rules.samtools_hostDepletion.output.fq,
        fq2 = rules.samtools_hostDepletion.output.fq2,
        singletons = rules.samtools_hostDepletion.output.singletons
    output:
        "metaspades/{sample}/contigs.fasta"
    params:
        DIR = "metaspades/{sample}",
        spades = "spades.py"
    threads: 8
    shell:
        "{params.spades} --meta -1 {input.fq} -2 {input.fq2} -s {input.singletons} -t {threads} -k 21,33,55,77,99,127 --only-assembler -o {params.DIR}"

rule seqkit_filterContigs:
    """Filter contigs by minimum length (300 bp)."""
    input:
        rules.metaspades_deNovoAssembly.output
    output:
        "filteredContigs/{sample}_filterContigs.fasta"
    shell:
        "seqkit seq -m 300 {input} -o {output}"

# =============================================================================
# SECTION 3: DIAMOND PROTEIN BLAST
# =============================================================================

rule dimaond_proteinBLAST:
    """Perform DIAMOND BLASTX against protein database."""
    input:
        rules.seqkit_filterContigs.output
    output:
        "diamond/{sample}.m8"
    params:
        db = DIAMOND_db,
        top = 5
    threads: 8
    shell:
        "diamond blastx -d {params.db} -q {input} -o {output} --top {params.top} --threads {threads} -b12 -c1"

# =============================================================================
# SECTION 4: TAXONOMIC ANNOTATION AND VIRAL TRIAGE
# =============================================================================

rule m8_1_filteringAndAnnotating:
    """Filter DIAMOND output and add taxonomic annotations using TaxonKit."""
    input:
        m8 = "diamond/{sample}.m8"
    output:
        m8_1a = temp("m8_1a_filterAndKeepTopHit/{sample}.m8_1aFiltered"),
        m8_1b = "m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        prefix = "{sample}"
    shell:
        """
        python {params.scripts_dir}/diamond_m8_filter.py {input.m8} {output.m8_1a}
        bash {params.scripts_dir}/add_taxon.sh -i {output.m8_1a} -f m8
        mv {params.prefix}.m8annotated {output.m8_1b}
        """

rule m8_2_triage:
    """Extract viral, bacterial, and eukaryotic hits from annotated DIAMOND results."""
    input:
        m8 = "m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated"
    output:
        m8_virus = "m8_2_triage/{sample}.m8virus",
        m8_bacteria = "m8_2_triage/{sample}.m8bacteria",
        m8_eukaryota = "m8_2_triage/{sample}.m8eukaryota"
    run:
        m8 = pd.read_csv(input.m8, sep='\t')
        m8_virus = m8.query('superkingdom == "Viruses"')
        m8_virus.to_csv(output.m8_virus, mode='w', index=False, header=True, sep='\t')
        m8_bacteria = m8.query('superkingdom == "Bacteria"')
        m8_bacteria.to_csv(output.m8_bacteria, mode='w', index=False, header=True, sep='\t')
        m8_eukaryota = m8.query('superkingdom == "Eukaryota"')
        m8_eukaryota.to_csv(output.m8_eukaryota, mode='w', index=False, header=True, sep='\t')

rule m8virus_removePhages:
    """Remove phages and other excluded viral families/species."""
    input:
        m8_virus = "m8_2_triage/{sample}.m8virus",
        exclude = config['exclude_txt']
    output:
        m8_nonphage = "m8_3_virus/{sample}.m8nonPhageVirus"
    shell:
        """
        cat {input.m8_virus} | csvtk grep -f family,species -t -T -v -i -r -P {input.exclude} -o {output.m8_nonphage}
        """

rule collect_viral_contigs4blastn:
    """Extract viral contig sequences for further BLASTN analysis."""
    input:
        m8 = "m8_3_virus/{sample}.m8nonPhageVirus"
    output:
        contig_list = temp("viral_contig/{sample}.txt"),
        fasta = "viral_contig/{sample}.fasta"
    params:
        contigs = "filteredContigs/{sample}_filterContigs.fasta"
    shell:
        """
        cat {input.m8} | csvtk cut -t -f qseqid | csvtk uniq | csvtk del-header -o {output.contig_list}
        seqkit grep -f {output.contig_list} {params.contigs} -o {output.fasta}
        """

# =============================================================================
# SECTION 5: BLASTN ANNOTATION AGAINST NT AND REFSEQ
# =============================================================================

rule blastn_and_filter_annotate:
    """BLASTN viral contigs against nt database and add taxonomy."""
    input:
        viral_contigs = "viral_contig/{sample}.fasta"
    output:
        blastn_filtered = "blastn/{sample}.blastnannotated"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        prefix = "{sample}"
    threads: 8
    shell:
        r"""
        blastn -db nt -query {input.viral_contigs} -out {wildcards.sample}.blastn.tmp -num_threads {threads} -max_target_seqs 10 -outfmt "6 qseqid sseqid sscinames staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore"
        python {params.scripts_dir}/blastn_filter.py {wildcards.sample}.blastn.tmp {wildcards.sample}.blastn.filtered.tmp --colname qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore -n 5
        bash {params.scripts_dir}/add_taxon.sh -i {wildcards.sample}.blastn.filtered.tmp -f blastn -k staxid -e -c 4 -n qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,library_id,accNum
        cat {params.prefix}.blastnannotated | csvtk filter2 -t -f '$superkingdom=="Viruses"' > {output.blastn_filtered}
        rm {wildcards.sample}.blastn.tmp {wildcards.sample}.blastn.filtered.tmp {params.prefix}.blastnannotated
        """

rule refseq_annotate:
    """BLASTN viral contigs against RefSeq viral database."""
    input:
        viral_contigs = "viral_contig/{sample}.fasta"
    output:
        blastn = "blastn/refseq/{sample}.refseqannotated"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        refseq_db = config['refseq_db'],
        refseq_segment_info = config['refseq_segment_info']
    threads: 8
    run:
        blastn_cmd = f"blastn -query {input.viral_contigs} -db {params.refseq_db} -out {wildcards.sample}.refseq.tmp -num_threads {threads} -word_size 11 -max_target_seqs 10 -outfmt '6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore'"
        subprocess.run(blastn_cmd, shell=True, check=True)

        filter_cmd = f"python {params.scripts_dir}/blastn_filter.py {wildcards.sample}.refseq.tmp {wildcards.sample}.refseq.filtered.tmp --colname qseqid,sseqid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"
        subprocess.run(filter_cmd, shell=True, check=True)

        blastn_df = pd.read_csv(f"{wildcards.sample}.refseq.filtered.tmp", sep='\t')
        if not blastn_df.empty:
            refseq_info = pd.read_csv(params.refseq_segment_info, sep='\t')
            final_df = pd.merge(blastn_df, refseq_info, left_on='accNum', right_on='acc', how='left').drop(columns=['acc'])
            final_df.to_csv(output.blastn, sep='\t', index=False)
        else:
            pd.DataFrame(columns=['library_id', 'qseqid', 'accNum', 'species', 'genus', 'family']).to_csv(output.blastn, sep='\t', index=False)

        subprocess.run(f"rm {wildcards.sample}.refseq.tmp {wildcards.sample}.refseq.filtered.tmp", shell=True)

checkpoint blastn_refseq_merge:
    """Merge BLASTN results from nt and RefSeq databases."""
    input:
        blastn = "blastn/{sample}.blastnannotated",
        refseq = "blastn/refseq/{sample}.refseqannotated"
    output:
        merged_blastn_refseq = "merged_blastn_refseq/{sample}.merged_blastn_refseq"
    run:
        blastn = pd.read_csv(input.blastn, sep='\t')
        refseq = pd.read_csv(input.refseq, sep='\t')

        if 'file_id' in blastn.columns:
            blastn.drop(columns=['file_id'], inplace=True)

        if len(blastn) > 0:
            blastn['segmented'] = blastn.apply(lambda x: virus_taxon_segment_info(x.get('species'), x.get('genus'), x.get('family')), axis=1)
        if len(refseq) > 0:
            refseq['segmented'] = refseq.apply(lambda x: virus_taxon_segment_info(x.get('species'), x.get('genus'), x.get('family')), axis=1)

        merged = pd.DataFrame()
        if len(blastn) > 0 and len(refseq) > 0:
            merged = pd.merge(blastn, refseq, on=['library_id', 'qseqid'], how='left', suffixes=('_b', '_r'))
            if 'family_b' in merged.columns and 'family_r' in merged.columns:
                initial_rows = len(merged)
                merged = merged[merged['family_r'].isnull() | (merged['family_b'] == merged['family_r'])]
                print(f"Filtered {initial_rows - len(merged)} rows with mismatched families for sample {wildcards.sample}.")
        # This part ensures that even if a contig is only found in one search,
        # its data is still kept and formatted consistently with the _b or _r suffixes for the downstream rules.
        elif len(blastn) > 0:
            merged = blastn.add_suffix('_b').rename(columns={'library_id_b': 'library_id', 'qseqid_b': 'qseqid'})
        elif len(refseq) > 0:
            merged = refseq.add_suffix('_r').rename(columns={'library_id_r': 'library_id', 'qseqid_r': 'qseqid'})
            rename_dict = {f'{col}_r': f'{col}_b' for col in refseq.columns if col not in ['library_id', 'qseqid']}
            merged.rename(columns=rename_dict, inplace=True)

        if merged.empty:
            merged = pd.DataFrame(columns=['library_id', 'qseqid'])

        merged.to_csv(output.merged_blastn_refseq, mode='w', index=False, header=True, sep='\t')

# =============================================================================
# SECTION 6: VIRAL CONTIG DEPTH CALCULATION AND TRIMMING
# =============================================================================

rule untrimmed_viral_contig_bowtie2:
    """Map reads to viral contigs to calculate depth."""
    input:
        fasta = "viral_contig/{sample}.fasta",
        fq = "hostDepleted/{sample}_1.fastq",
        fq2 = "hostDepleted/{sample}_2.fastq",
        singletons = "hostDepleted/{sample}_singletons.fastq"
    output:
        bam = "untrimmed_viral_contig_bowtie2/{sample}.bam",
        bam_index = "untrimmed_viral_contig_bowtie2/{sample}.bam.bai",
        idxstats = "untrimmed_viral_contig_bowtie2/{sample}.idxstats"
    params:
        samtools = config['samtools']
    threads: 2
    run:
        import os
        import subprocess
        if os.path.getsize(input.fasta) == 0:
            with open(output.bam, 'w') as f: f.write("")
            with open(output.bam_index, 'w') as f: f.write("")
            with open(output.idxstats, 'w') as f: f.write("")
        else:
            # Create output dir if it doesn't exist
            os.makedirs(os.path.dirname(output.bam), exist_ok=True)
            # Using shell=True for simplicity with piping
            build_cmd = f"bowtie2-build {input.fasta} {input.fasta}"
            map_cmd = f"bowtie2 --threads 1 -q -x {input.fasta} -1 {input.fq} -2 {input.fq2} -U {input.singletons} | {params.samtools} view -bS | {params.samtools} sort -o {output.bam}"
            subprocess.run(build_cmd, shell=True, check=True, capture_output=True)
            subprocess.run(map_cmd, shell=True, check=True)
            subprocess.run([params.samtools, "index", output.bam], check=True)
            subprocess.run([params.samtools, "idxstats", output.bam], stdout=open(output.idxstats, 'w'), check=True)

rule viral_contig_depth_initial:
    """Calculate depth coverage for viral contigs."""
    input:
        bam = "untrimmed_viral_contig_bowtie2/{sample}.bam"
    output:
        depth = "untrimmed_viral_contig_bowtie2/{sample}.depth"
    shell:
        "bedtools genomecov -bga -ibam {input.bam} > {output.depth}"

rule viral_contig_trim:
    """Trim viral contigs based on depth threshold."""
    input:
        depth = "untrimmed_viral_contig_bowtie2/{sample}.depth",
        fasta = "viral_contig/{sample}.fasta"
    output:
        trimmed_fasta = "trimmed_viral_contig/{sample}.fasta",
        trimming_summary = "trimmed_viral_contig/{sample}.summary"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        depth = config.get('min_depth', 10)
    shell:
        "python {params.scripts_dir}/trim_contig.py {input.depth} {input.fasta} -o {output.trimmed_fasta} -d {params.depth} -s {output.trimming_summary}"

# =============================================================================
# SECTION 7: SUMMARY STATISTICS GENERATION
# =============================================================================

rule viral_read_statistics:
    """Generate per-sample viral contig statistics at genus and family levels."""
    input:
        merged_blastn_refseq = "merged_blastn_refseq/{sample}.merged_blastn_refseq",
        depth = "untrimmed_viral_contig_bowtie2/{sample}.depth",
        idxstats = "untrimmed_viral_contig_bowtie2/{sample}.idxstats",
        trimming_summary = "trimmed_viral_contig/{sample}.summary",
        fq = "hostDepleted/{sample}_1.fastq",
        fq2 = "hostDepleted/{sample}_2.fastq",
        singletons = "hostDepleted/{sample}_singletons.fastq"
    output:
        genus_summary = "viral_contig_summary/genus/{sample}.viral_contig_summary",
        family_summary = "viral_contig_summary/family/{sample}.viral_contig_summary"
    params:
        scripts_dir = config.get('scripts_dir', os.path.join(os.path.dirname(workflow.snakefile), 'scripts')),
        sample_id = "{sample}"
    shell:
        """
        python {params.scripts_dir}/summarise_contigs_per_taxon.py {input.merged_blastn_refseq} {input.depth} {input.idxstats} -t {input.trimming_summary} -o {output.genus_summary} -p {params.sample_id} -c genus_b --fq {input.fq},{input.fq2},{input.singletons}
        python {params.scripts_dir}/summarise_contigs_per_taxon.py {input.merged_blastn_refseq} {input.depth} {input.idxstats} -t {input.trimming_summary} -o {output.family_summary} -p {params.sample_id} -c family_b --fq {input.fq},{input.fq2},{input.singletons}
        """

rule batch_viral_read_statistics:
    """Combine all sample summaries into batch-level reports."""
    input:
        family = expand("viral_contig_summary/family/{sample}.viral_contig_summary", sample=SAMPLES),
        genus = expand("viral_contig_summary/genus/{sample}.viral_contig_summary", sample=SAMPLES)
    output:
        family_batch = expand("viral_contig_summary/batch/{batch}.family.viral_contig_summary_all_samples.tsv", batch = config.get('batch', 'batch')),
        genus_batch = expand("viral_contig_summary/batch/{batch}.genus.viral_contig_summary_all_samples.tsv", batch = config.get('batch', 'batch'))
    run:
        family_df = pd.DataFrame()
        for f in input.family:
            df = pd.read_csv(f, sep="\t", header=0)
            family_df = pd.concat([family_df, df], ignore_index=True)
        family_df.to_csv(output.family_batch[0], sep="\t", index=False)

        genus_df = pd.DataFrame()
        for f in input.genus:
            df = pd.read_csv(f, sep="\t", header=0)
            genus_df = pd.concat([genus_df, df], ignore_index=True)
        genus_df.to_csv(output.genus_batch[0], sep="\t", index=False)

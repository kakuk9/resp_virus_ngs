"""
Snakemake workflow name: mngs_pre_diamond.smk
Last updated by Kirsty Kwok on 2024-03-07
Description: This snakemake script covers steps from adapter trimming to protein BLAST using DIAMOND.
NB: Run this script in the directory where you want the resulting files to be saved in!
Usage: snakemake --snakefile mngs_pre_diamond.smk --configfile YourConfigYamlFile --cores 8
"""


### Set shell and import packages ###
shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc;")

import os
import csv
import glob
import gzip
import pandas as pd

### config file ###
RAW_FASTQ_DIR = config['rawFastqDir']
HOST_REF = config['hostRef']
DIAMOND_db=config['DiamondDB']
num_underscore=config['numUnderscore']
Illumina_sequencer=config['IlluminaSequencer']
samtools=config['samtools']

### presets ###
fastq_suffix=['1', '2', 'singletons']

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

### Get sample name ###
def list_samples(RawFastqDirectory, RawFastq_suffix):
    SAMPLES=[]
    PATH=RawFastqDirectory+RawFastq_suffix
    print(f"PATH {PATH}")
    for file in glob.glob(PATH):
        base=os.path.basename(file)
        #Take up to three underscores as the sample name e.g. NT_4_SS_S1 -> NT_4_S1. Revise line below if needed.
        sample=""
        #sample = (base.split("_")[0]+"_"+base.split("_")[1]+"_"+base.split("_")[2])
        for i in range(num_underscore):
            if i!=num_underscore-1:
                sample+=base.split("_")[i]+"_"
            else:
                sample+=base.split("_")[i]
        # don't analyse the undetermined samples
        if not "Undetermined" in sample:
            SAMPLES.append(sample)
    return(SAMPLES)

SAMPLES=list_samples(RAW_FASTQ_DIR, RawFastq_suffix)
print(list_samples(RAW_FASTQ_DIR, RawFastq_suffix))

rule all:
    input:
        expand("trimGalore/{sample}_trimmed_{suffix}.fq", sample=SAMPLES, suffix=fastq_suffix),
        expand("bowtie2_4hostDepletion/{sample}_host.bam", sample=SAMPLES),
        expand("hostDepleted/{sample}_{suffix}.fastq", sample=SAMPLES, suffix=fastq_suffix),
        expand("metaspades/{sample}/contigs.fasta", sample=SAMPLES),
        expand("filteredContigs/{sample}_filterContigs.fasta", sample=SAMPLES),
        expand("diamond/{sample}.m8", sample=SAMPLES)

ruleorder: trimGalore_adapterTrimming  > bowtie2_4hostDepletion > samtools_hostDepletion > metaspades_deNovoAssembly > seqkit_filterContigs > dimaond_proteinBLAST

rule trimGalore_adapterTrimming:
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
    input:
        bam = rules.bowtie2_4hostDepletion.output.bam
    output:
        fq="hostDepleted/{sample}_1.fastq",
        fq2="hostDepleted/{sample}_2.fastq",
        singletons="hostDepleted/{sample}_singletons.fastq"
    params:
        samtools = samtools
    threads: 8
    shell:
        """
        {params.samtools} view -b -f4 {input.bam} | {params.samtools} fastq -@ {threads} -1 {output.fq} -2 {output.fq2} -s {output.singletons} -0 /dev/null -n   
        """

rule metaspades_deNovoAssembly:
    input:
        fq = rules.samtools_hostDepletion.output.fq,
        fq2 = rules.samtools_hostDepletion.output.fq2,
        singletons = rules.samtools_hostDepletion.output.singletons
    output:
        "metaspades/{sample}/contigs.fasta"
    params:
        dir="metaspades/{sample}",
        spades = "/home3/2893911k/miniforge3/envs/kirsty_env/bin/spades.py"
    threads: 8
    shell:
        "{params.spades} --meta -1 {input.fq} -2 {input.fq2} -s {input.singletons} -t {threads} -k 21,33,55,77,99,127 --only-assembler -o {params.dir}"

rule seqkit_filterContigs:
    input:
        rules.metaspades_deNovoAssembly.output
    output:
        "filteredContigs/{sample}_filterContigs.fasta"
    shell:
        "seqkit seq -m 300 {input} -o {output}"

rule dimaond_proteinBLAST:
    input:
        rules.seqkit_filterContigs.output
    output:
        "diamond/{sample}.m8"
    params:
        db=DIAMOND_db,
        top=5
    threads: 8
    shell:
        "diamond blastx -d {params.db} -q {input} -o {output} --top {params.top} --threads {threads} -b12 -c1"

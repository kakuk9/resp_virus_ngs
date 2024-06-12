"""
"""

import os
import glob
import pandas as pd
import re

blastn = config["segment_virus_blastn"]
blastn = pd.read_csv(blastn, sep="\t")
### Get sample name ###
def list_samples(sample_tsv):
    SAMPLES = set(sample_tsv["file_id"])
    return(SAMPLES)

SAMPLES=list_samples(blastn)
print(SAMPLES)
rule all:
    input:
        expand("segmentedViruses/extracted_contigs/{sample}.{suffix}", sample = SAMPLES, suffix=['id', 'fasta']),
        expand("segmentedViruses/refseq_blastn/{sample}.{suffix2}", sample = SAMPLES, suffix2=['refseq_viral', 'refseq_viral_raw', 'refseq_viral_annotated'])

rule extract_contigs_id:
    input:
        blastn = config["segment_virus_blastn"]
    output:
        "segmentedViruses/extracted_contigs/{sample}.id"
    params:
        sample = "{sample}"
    run:
        blastn = pd.read_csv(input.blastn, sep="\t")
        df = blastn[blastn['file_id'] == params.sample]
        df = df[['qseqid']].reset_index(drop=True)
        df.to_csv(output[0], sep='\t', index=False, header=False)

rule extract_contigs_segmented:
    input:
        id_list = "segmentedViruses/extracted_contigs/{sample}.id",
        contigs = "metaspades/{sample}/contigs.fasta"
    output:
        fasta = "segmentedViruses/extracted_contigs/{sample}.fasta"
    shell:
        "seqkit grep -f {input.id_list} {input.contigs} > {output.fasta}"

rule refseq_viral:
    input:
        segmented_contigs = "segmentedViruses/extracted_contigs/{sample}.fasta"
    output:
        blastn = "segmentedViruses/refseq_blastn/{sample}.refseq_viral_raw",
        blastn_filtered = "segmentedViruses/refseq_blastn/{sample}.refseq_viral"
    params:
        refseq_db = config['refseq_db'],
        blastn_filter_script = config['blastn_filter_script']
    shell:
        """
        blastn -query {input.segmented_contigs} -db {params.refseq_db} -out {output.blastn} -num_threads 8 -max_target_seqs 10 -outfmt '6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore'
        python {params.blastn_filter_script} {output.blastn} {output.blastn_filtered} --colname qseqid,sseqid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
        """
rule refseq_viral_annotate:
	input: 
		blastn="segmentedViruses/refseq_blastn/{sample}.refseq_viral"
	output:
		blastn="segmentedViruses/refseq_blastn/{sample}.refseq_viral_annotated"
	params:
		refseq_segment_info = config['refseq_segment_info'],
		prefix="{sample}"
	run:
		blastn = pd.read_csv(input.blastn,sep='\t')
		if len(blastn) > 0:
			refseq_segment_info = pd.read_csv(params.refseq_segment_info,sep='\t')
			blastn_final = pd.merge(blastn, refseq_segment_info, left_on='accNum', right_on='acc', how='left')
			blastn_final.to_csv(output.blastn, mode='w',index=False,header=True,sep='\t')
		else:
			empty_df = pd.DataFrame()
			empty_df['library_id'] = params.prefix
			empty_df.to_csv(output.blastn, mode='w',index=False,header=True,sep='\t')

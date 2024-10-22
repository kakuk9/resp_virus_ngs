"""
Snakemake workflow name: post_diamond.smk
Description: This snakemake workflow annotate diamond m8 and filter
NB: Run this script in the directory where you want the resulting files to be saved in.
Usage: snakemake --snakefile post_diamond.smk --configfile YourConfigYamlFile --cores 2
Written by Kirsty Kwok, last updated on 15/12/2023
"""

#shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc;")

import os
import csv
import glob
import pandas as pd

INPUT_DIR = config["m8_dir"]
CONTIGS_DIR = config["contigs_directory"]
VIRUS_TAXON_SEGMENT = pd.read_csv(config["virus_segment_info"], sep="\t")
VIRUS_TAXON_SEGMENT=dict(zip(VIRUS_TAXON_SEGMENT['taxon_name'],VIRUS_TAXON_SEGMENT['segmented']))

### Get sample name ###
def list_samples(DIR_M8):
	SAMPLES=[]
	for file in glob.glob(DIR_M8+"/*.m8"):
		base=os.path.basename(file)
        #Take up to the whole name before the "." extension
		sample = base.split(".")[0]
		if not "Undetermined" in sample:
			SAMPLES.append(sample)
	return(SAMPLES)

SAMPLES=list_samples(INPUT_DIR)
print(list_samples(INPUT_DIR))

### Useful functions ###
def virus_taxon_segment_info(SPECIES, GENUS, FAMILY):
	if SPECIES in VIRUS_TAXON_SEGMENT.keys():
		return VIRUS_TAXON_SEGMENT[SPECIES]
	elif GENUS in VIRUS_TAXON_SEGMENT.keys():
		return VIRUS_TAXON_SEGMENT[GENUS]
	elif FAMILY in VIRUS_TAXON_SEGMENT.keys():
		return VIRUS_TAXON_SEGMENT[FAMILY]
	else:
		return "Unknown"

rule all:
	input:
		expand("m8_1a_filterAndKeepTopHit/{sample}.m8_1aFiltered", sample=SAMPLES),
		expand("m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated", sample=SAMPLES),
		expand("m8_2_triage/{sample}.m8{cat}", sample=SAMPLES, cat=['virus','bacteria','eukaryota']),
		expand("m8_1a_filterAndKeepTopHit/{sample}.m8_1aFiltered", sample=SAMPLES),
		expand("m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated", sample=SAMPLES),
		expand("m8_2_triage/{sample}.m8{cat}", sample=SAMPLES, cat=['virus','bacteria','eukaryota']),
		expand("m8_2_triage/batch/{batch}_batchm8{cat}.tsv", batch=config['batch'], cat=['virus','bacteria','eukaryota']),
		expand("m8_3_virus/{sample}.m8nonPhageVirus", sample=SAMPLES),
		expand("m8_3_virus/{sample}.m8phage", sample=SAMPLES),
		expand("viral_contig/{sample}.txt",sample=SAMPLES),
		expand("viral_contig/{sample}.fasta", sample=SAMPLES),
		expand("raw_blastn/{sample}.blastn", sample=SAMPLES),
		expand("blastn/{sample}.blastnannotated", sample=SAMPLES),
		expand("raw_blastn/refseq/{sample}.refseq", sample=SAMPLES),
		expand("blastn/refseq/{sample}.refseqannotated", sample=SAMPLES),
		expand("merged_blastn_refseq/{sample}.merged_blastn_refseq", sample=SAMPLES),
		expand("viral_contig_bowtie2/{sample}.bam", sample=SAMPLES),
		expand("viral_contig_bowtie2/{sample}.depth", sample=SAMPLES),
		expand("viral_contig_summary/genus/{sample}.viral_contig_summary", sample=SAMPLES),
		expand("viral_contig_summary/family/{sample}.viral_contig_summary", sample=SAMPLES)
		#expand("trimmed_viral_contig/{sample}{suffix_t}", sample=SAMPLES, suffix_t=['_trimming.txt','.fasta']),


rule m8_1_filteringAndAnnotating:
	input:
		m8=f"{INPUT_DIR}/{{sample}}.m8"
	output:
		m8_1a="m8_1a_filterAndKeepTopHit/{sample}.m8_1aFiltered",
		m8_1b="m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated"
	params:
		python_script=config['m8_keepTop_script'],
		bash=config['add_taxon_script'],
		prefix="{sample}"
	shell:
		"""
		python {params.python_script} {input.m8} {output.m8_1a}
		bash {params.bash} -i {output.m8_1a} -f m8
		mv {params.prefix}.m8annotated {output.m8_1b}
		"""

rule m8_2_triage:
	input:
		m8 = "m8_1b_annotateTaxonkit/{sample}.m8_1bAnnotated"
	output:
		m8_virus = "m8_2_triage/{sample}.m8virus",
		m8_bacteria = "m8_2_triage/{sample}.m8bacteria",
		m8_eukaryota = "m8_2_triage/{sample}.m8eukaryota"
	run:
		m8 = pd.read_csv(input.m8,sep='\t')
		m8_virus=m8.query('superkingdom == "Viruses"')
		m8_virus.to_csv(output.m8_virus, mode='w',index=False,header=True,sep='\t')
		m8_bacteria=m8.query('superkingdom == "Bacteria"')
		m8_bacteria.to_csv(output.m8_bacteria, mode='w',index=False,header=True,sep='\t')
		m8_eukaryota=m8.query('superkingdom == "Eukaryota"')
		m8_eukaryota.to_csv(output.m8_eukaryota, mode='w',index=False,header=True,sep='\t')

rule m8_2_triage_batch:
	input:
		virus = expand("m8_2_triage/{sample}.m8virus", sample=SAMPLES),
		bacteria = expand("m8_2_triage/{sample}.m8bacteria", sample=SAMPLES),
		eukaryota = expand("m8_2_triage/{sample}.m8eukaryota", sample=SAMPLES)
	output:
		m8_batch_virus = "m8_2_triage/batch/{batch}_batchm8virus.tsv",
		m8_batch_bacteria = "m8_2_triage/batch/{batch}_batchm8bacteria.tsv",
		m8_batch_eukaryota = "m8_2_triage/batch/{batch}_batchm8eukaryota.tsv"
	run:
		virus = pd.concat([pd.read_csv(f,sep='\t') for f in input.virus])
		virus.to_csv(output.m8_batch_virus, mode='w',index=False,header=True,sep='\t')
		bacteria = pd.concat([pd.read_csv(f,sep='\t') for f in input.bacteria])
		bacteria.to_csv(output.m8_batch_bacteria, mode='w',index=False,header=True,sep='\t')
		eukaryota = pd.concat([pd.read_csv(f,sep='\t') for f in input.eukaryota])
		eukaryota.to_csv(output.m8_batch_eukaryota, mode='w',index=False,header=True,sep='\t')

rule m8virus_removePhages:
	input:
		m8_virus = "m8_2_triage/{sample}.m8virus",
		exclude = config['exclude_txt']
	output:
		m8_nonphage = "m8_3_virus/{sample}.m8nonPhageVirus",
		m8_phage = "m8_3_virus/{sample}.m8phage",
	shell:
		"""
		cat {input.m8_virus} | csvtk grep -f family,species -t -T -v -i -r -P  {input.exclude} -o {output.m8_nonphage}
		cat {input.m8_virus} | csvtk grep -f family,species -t -T -i -r -P  {input.exclude} -o {output.m8_phage}
		"""

rule collect_viral_contigs4blastn:
	input:
		m8 = "m8_3_virus/{sample}.m8nonPhageVirus"
	output:
		contig_list = temp("viral_contig/{sample}.txt1"),
		fasta =  temp("viral_contig/{sample}.fasta1")
	params:
		contigs = f"{CONTIGS_DIR}/{{sample}}_filterContigs.fasta"
	shell:
		"""
		cat {input.m8} | csvtk cut -t -f qseqid | csvtk uniq | csvtk del-header -o {output.contig_list}
		seqkit grep -f {output.contig_list} {params.contigs} -o {output.fasta}
		"""

rule blastn_and_filter_annotate:
	input:
		viral_contigs="viral_contig/{sample}.fasta1"
	output:
		blastn="raw_blastn/{sample}.blastn",
		blastn_tmp=temp("raw_blastn/{sample}.blastn.tmp"),
		blastn_filtered_tmp=temp("blastn/{sample}.blastnannotated.tmp"),
		blastn_filtered="blastn/{sample}.blastnannotated",
		blastn_non_viral="blastn/non_viral/{sample}.blastnannotated"
	params:
		blastn_filtered = config['blastn_filter_script'],
		bash=config['add_taxon_script'],
		prefix="{sample}"
	shell:
		"""
		blastn -db nt -query {input.viral_contigs} -out {output.blastn} -num_threads {threads} -max_target_seqs 10 -outfmt "6 qseqid sseqid sscinames staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" 
		python {params.blastn_filtered} {output.blastn} {output.blastn_tmp} --colname qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore -n 5
		bash {params.bash} -i {output.blastn_tmp} -f blastn -k staxid -e -c 4 -n qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,library_id,accNum
		mv {params.prefix}.blastnannotated {output.blastn_filtered_tmp}
		cat {output.blastn_filtered_tmp} | csvtk filter2 -t -f '$superkingdom=="Viruses"' > {output.blastn_filtered}
		cat {output.blastn_filtered_tmp} | csvtk filter2 -t -f '$superkingdom!="Viruses"' > {output.blastn_non_viral}
		"""

rule extract_viralContigs4refseq:
	input:
		blastn="blastn/{sample}.blastnannotated",
		fasta1="viral_contig/{sample}.fasta1"
	output:
		contig_list = "viral_contig/{sample}.txt",
		fasta = "viral_contig/{sample}.fasta"
	shell:
		"""
		cat {input.blastn} | csvtk cut -tf qseqid | csvtk del-header | csvtk uniq > {output.contig_list}
		seqkit grep -f {output.contig_list} {input.fasta1} -o {output.fasta}
		"""

rule refseqviral_and_filter_annotate:
	input:
		viral_contigs="viral_contig/{sample}.fasta"
	output:
		refseq_blastn = "raw_blastn/refseq/{sample}.refseq",
		refseq_blastn_tmp = temp("raw_blastn/refseq/{sample}.refseq.tmp")
	params:
		blastn_filtered = config['blastn_filter_script'],
		bash=config['add_taxon_script'],
		refseq_db = config['refseq_db'],
		prefix="{sample}"
	threads: 8
	shell:
		"""
		blastn -query {input.viral_contigs} -db {params.refseq_db} -out {output.refseq_blastn} -num_threads 8 -word_size 11 -max_target_seqs 10 -outfmt '6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore'
		python {params.blastn_filtered} {output.refseq_blastn} {output.refseq_blastn_tmp} --colname qseqid,sseqid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
		"""

rule refseq_annotate:
	input: 
		blastn="raw_blastn/refseq/{sample}.refseq.tmp"
	output:
		blastn="blastn/refseq/{sample}.refseqannotated"
	params:
		refseq_segment_info = config['refseq_segment_info']
	run:
		blastn = pd.read_csv(input.blastn,sep='\t')
		if len(blastn) > 0:
			refseq_segment_info = pd.read_csv(params.refseq_segment_info,sep='\t')
			blastn_final = pd.merge(blastn, refseq_segment_info, left_on='accNum', right_on='acc', how='left')
			blastn_final = blastn_final.drop(columns=['acc'])
			blastn_final.to_csv(output.blastn, mode='w',index=False,header=True,sep='\t')
		else:
			empty_df = pd.DataFrame(columns=['library_id', 'qseqid'])
			empty_df.to_csv(output.blastn, mode='w',index=False,header=True,sep='\t')

rule blastn_refseq_merge:
	input:
		blastn = "blastn/{sample}.blastnannotated",
		refseq = "blastn/refseq/{sample}.refseqannotated"
	output:
		merged_blastn_refseq = "merged_blastn_refseq/{sample}.merged_blastn_refseq"
	run:
		blastn = pd.read_csv(input.blastn,sep='\t')
		blastn = blastn.drop(columns=['file_id'])
		refseq = pd.read_csv(input.refseq,sep='\t')
		if len(blastn) > 0:
			blastn['segmented'] = blastn.apply(lambda x: virus_taxon_segment_info(x['species'], x['genus'], x['family']), axis=1)
		if len(refseq) > 0:
			refseq['segmented'] = refseq.apply(lambda x: virus_taxon_segment_info(x['species'], x['genus'], x['family']), axis=1)
		if len(blastn) > 0 and len(refseq) > 0:
			merged = blastn.merge(refseq, on=['library_id','qseqid'], how='left', suffixes=('_b', '_r'))
		elif len(blastn) > 0:
			merged = blastn.add_suffix('_b').rename(columns={'library_id_b':'library_id', 'qseqid_b':'qseqid', 'sscinames_b':'sscinames', 'staxid_b':'staxid'})
		elif len(refseq) > 0:
			merged = refseq.add_suffix('_r').rename(columns={'library_id_r':'library_id', 'qseqid_r':'qseqid'})
		else:
			merged = pd.DataFrame(columns=['library_id', 'qseqid'])
		merged.to_csv(output.merged_blastn_refseq, mode='w',index=False,header=True,sep='\t')

rule viral_contig_bowtie2:
	input:
		contig_list = "viral_contig/{sample}.txt",
		fasta = "viral_contig/{sample}.fasta",
		fq = "hostDepleted/{sample}_1.fastq",
        fq2 = "hostDepleted/{sample}_2.fastq",
		singletons = "hostDepleted/{sample}_singletons.fastq"
	output:
		bam = "viral_contig_bowtie2/{sample}.bam"
	params:
		samtools = config['samtools']
	threads:
		8
	run:
		import os
		import subprocess
		if os.path.getsize(input.fasta) == 0:
			with open(output.bam, 'w') as f:
				f.write("")
		else:
			subprocess.run(["bowtie2-build", input.fasta, input.fasta])
			bowtie2_aln = subprocess.Popen(["bowtie2", "-q", "-x", input.fasta, "-1", input.fq, "-2", input.fq2, "-U", input.singletons, "--threads", "8"], stdout=subprocess.PIPE)
			samtools_1 = subprocess.Popen([params.samtools, "view", "-@", "8", "-bS"], stdin=bowtie2_aln.stdout, stdout=subprocess.PIPE)
			subprocess.run([params.samtools, "sort", "-@", "8", "-o", output.bam], stdin=samtools_1.stdout)
			subprocess.run([params.samtools, "index", "-@", "8", output.bam])

rule viral_contig_depth:
	input:
		bam = "viral_contig_bowtie2/{sample}.bam"
	output:
		depth = "viral_contig_bowtie2/{sample}.depth"
	shell:
		"""
		bedtools genomecov -bga -ibam {input.bam} > {output.depth}
		"""
rule viral_read_statistics:
	input:
		blastn = "blastn/{sample}.blastnannotated",
		bam = "viral_contig_bowtie2/{sample}.bam",
		depth = "viral_contig_bowtie2/{sample}.depth"
	output:
		genus_summary = "viral_contig_summary/genus/{sample}.viral_contig_summary",
		family_summary = "viral_contig_summary/family/{sample}.viral_contig_summary"
	params:
		viral_read_statistics = config['viral_read_statistics'],
		sample_id = "{sample}"
	shell:
		"""
		python {params.viral_read_statistics} {input.blastn} {input.bam} {input.depth} -o {output.genus_summary} -p {params.sample_id} -c genus
		python {params.viral_read_statistics} {input.blastn} {input.bam} {input.depth} -o {output.family_summary} -p {params.sample_id} -c family
		"""
# rule trim_viral_contig:
# 	input:
# 		depth = "viral_contig_bowtie2/{sample}.depth",
# 		fasta = "viral_contig/{sample}.fasta"
# 	output:
# 		fasta = "trimmed_viral_contig/{sample}.fasta",
# 		trimming_summary = "trimmed_viral_contig/{sample}_trimming.txt"
# 	params:
# 		trim_viral_contig = config['trim_viral_contig']
# 	shell:
# 		"""
# 		python {params.trim_viral_contig} {input.depth} {input.fasta} {output.fasta} {output.trimming_summary}
# 		"""


# rule summarise_virus:
# 	input:
# 		merged = "merged_blastn_refseq/{sample}.merged_blastn_refseq"
# 	output:
# 		summary = "summary/{sample}.summary"
# 	run:
		

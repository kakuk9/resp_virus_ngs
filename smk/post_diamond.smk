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

rule all:
	input:
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
		expand("all_blastn/{batch}_{cat2}.tsv", batch=config['batch'], cat2=["batchBlastnAnnotated","segmented","nonsegmented","unknownGenomeStructure"]),
		expand("reference_seq/{batch}_moreThan2Contigs.tsv", batch=config['batch'])

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
		list = "viral_contig/{sample}.txt",
		fasta =  "viral_contig/{sample}.fasta"
	params:
		contigs = f"{CONTIGS_DIR}/{{sample}}_filterContigs.fasta"
	shell:
		"""
		cat {input.m8} | csvtk cut -t -f qseqid | csvtk uniq | csvtk del-header -o {output.list}
		seqkit grep -f {output.list} {params.contigs} -o {output.fasta}
		"""


rule blastn_and_filter_annotate:
	input:
		viral_contigs="viral_contig/{sample}.fasta"
	output:
		blastn="raw_blastn/{sample}.blastn",
		blastn_tmp=temp("raw_blastn/{sample}.blastn.tmp"),
		blastn_filtered="blastn/{sample}.blastnannotated"
	params:
		blastn_filtered = config['blastn_filter_script'],
		bash=config['add_taxon_script'],
		prefix="{sample}"
	threads: 8
	shell:
		"""
		blastn -db nt -query {input.viral_contigs} -out {output.blastn} -num_threads {threads} -max_target_seqs 10 -outfmt "6 qseqid sseqid sscinames staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" 
		python {params.blastn_filtered} {output.blastn} {output.blastn_tmp} --colname qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
		bash {params.bash} -i {output.blastn_tmp} -f blastn -k staxid -e -c 4 -n qseqid,sseqid,sscinames,staxid,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,library_id,accNum
		mv {params.prefix}.blastnannotated {output.blastn_filtered}
		"""

rule combine_all_blastn:
	input:
		blastn=expand("blastn/{sample}.blastnannotated", sample=SAMPLES)
	output:
		blastn="all_blastn/{batch}_batchBlastnAnnotated.tsv",
		segmented = "all_blastn/{batch}_segmented.tsv",
		nonsegmented = "all_blastn/{batch}_nonsegmented.tsv",
		unknownGenomeStructure = "all_blastn/{batch}_unknownGenomeStructure.tsv"
	params:
		family = config['virus_family_segmentInfo'],
		genus = config['virus_genus_segmentInfo'],
		species = config['virus_species_segmentInfo']
	run:
		blastn = pd.DataFrame()
		for f in input.blastn:
			df = pd.read_csv(f,sep='\t')
			blastn = pd.concat([blastn,df if not df.empty else None])
		blastn = blastn.query('superkingdom == "Viruses"').reset_index(drop=True)
		virus_family_segment_orNot = pd.read_csv(params.family, sep='\t')
		virus_family_segment_orNot = dict(zip(virus_family_segment_orNot['family'], virus_family_segment_orNot['segmented']))
		virus_genus_segment_orNot = pd.read_csv(params.genus, sep='\t')
		virus_genus_segment_orNot = dict(zip(virus_genus_segment_orNot['genus'], virus_genus_segment_orNot['segmented']))
		virus_species_segment_orNot = pd.read_csv(params.species, sep='\t')
		virus_species_segment_orNot = dict(zip(virus_species_segment_orNot['species'], virus_species_segment_orNot['segmented']))
		def segment_or_not(SPECIES, GENUS, FAMILY):
			if SPECIES in virus_species_segment_orNot:
				return virus_species_segment_orNot[SPECIES]
			elif GENUS in virus_genus_segment_orNot:
				return virus_genus_segment_orNot[GENUS]
			elif FAMILY in virus_family_segment_orNot:
				return virus_family_segment_orNot[FAMILY]
			else:
				return "unknown"
		blastn['segmented'] = blastn.apply(lambda x: segment_or_not(x['species'], x['genus'], x['family']), axis=1)
		blastn.to_csv(output.blastn, mode='w',index=False,header=True,sep='\t')
		segmented = blastn.query('segmented == "yes"')
		segmented.to_csv(output.segmented, mode='w',index=False,header=True,sep='\t')
		nonsegmented = blastn.query('segmented == "no"')
		nonsegmented.to_csv(output.nonsegmented, mode='w',index=False,header=True,sep='\t')
		unknownGenomeStructure = blastn.query('segmented == "unknown"')
		unknownGenomeStructure.to_csv(output.unknownGenomeStructure, mode='w',index=False,header=True,sep='\t')

rule summarise_virus:
	input:
		blastn="all_blastn/{batch}_batchBlastnAnnotated.tsv"
	output:
		tsv = "reference_seq/{batch}_moreThan2Contigs.tsv"
	params:
		script = config['summarise_virus_scaffold_script'],
		batch = config['batch']
	shell:
		"""
		python {params.script} {input.blastn} -p {params.batch}
		"""

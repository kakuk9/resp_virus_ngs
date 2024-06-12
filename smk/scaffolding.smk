"""
"""

import os
import glob
import pandas as pd
import re

#summary_dir = config['summary_dir']
summary_dir = "/home3/2893911k/AH001/Batch1/NEW_NextSeqHigh/reference_seq"
#spades_directory = config['spades_directory']
spades_directory = "/home3/2893911k/AH001/Batch1/NEW_NextSeqHigh/metaspades"
parsing_scaffoldBuilder_script = "/home3/2893911k/kirsty_scripts/python/parsing_scaffoldBuilder.py"
batch_blastn = "/home3/2893911k/AH001/Batch1/NEW_NextSeqHigh/all_blastn/NEW_NextSeqHigh_batchBlastnAnnotated.tsv"

### Get sample name ###
def list_samples(dir):
    SCAFFOLDS =[]
    for file in os.listdir(dir):
        base=os.path.basename(file)
        if re.search(r"SCAFFOLD", base):
            sample = "_".join(base.split("_")[0:3])
            virus_species = "".join(base.split("_SCAFFOLD_")[1].split(".")[0])
            scaffold = sample+"_SCAFFOLD_"+virus_species
            SCAFFOLDS.append(scaffold)
    return(SCAFFOLDS)

SCAFFOLDS=list_samples(dir=summary_dir)
#print(SCAFFOLDS)

rule all:
    input:
        expand("scaffolding/extracted_contigs/{scaffold_id}.fasta", scaffold_id = SCAFFOLDS),
        expand("scaffolding/scaffoldBuilder/{scaffold_id}_Scaffold.fasta", scaffold_id = SCAFFOLDS),
        expand("scaffolding/parsing_scaffoldBuilder_output/{scaffold_id}_sequence_class.tsv", scaffold_id = SCAFFOLDS),
        "batch_scaffold_info.tsv"


rule extract_contigs:
    input:
        list_of_contigs = "reference_seq/{scaffold_id}.txt"
    output:
        fasta = "scaffolding/extracted_contigs/{scaffold_id}.fasta"
    params:
        contig_dir = f"{spades_directory}",
        sample = lambda wildcards: wildcards.scaffold_id.split("_SCAFFOLD_")[0]
    shell:
        """
        seqkit grep -f {input.list_of_contigs} {params.contig_dir}/{params.sample}/contigs.fasta > {output.fasta}
        """

rule scaffold_builder:
    input:
        fasta = "scaffolding/extracted_contigs/{scaffold_id}.fasta",
        ref = "reference_seq/{scaffold_id}.fasta"
    output:
        fasta = "scaffolding/scaffoldBuilder/{scaffold_id}_Scaffold.fasta"
    params:
        prefix = "{scaffold_id}",
        outputDir = "scaffolding/scaffoldBuilder/"
    run:
        import glob
        import os
        import subprocess
        subprocess.run(["mamba", "run", "-n", "scaffold_builder", "scaffold_builder.py", "-q", input.fasta, "-r", input.ref, "-p", params.prefix])
        subprocess.run(["mv", f"{params.prefix}_Scaffold.fasta",output.fasta])
        subprocess.run(["mv", f"{params.prefix}_overlap_alignment",params.outputDir])
        subprocess.run(["mv", f"{params.prefix}.coords",params.outputDir])
        subprocess.run(["mv", f"{params.prefix}_output.txt",params.outputDir])

rule parse_scaffold:
    input:
        expand("scaffolding/scaffoldBuilder/{scaffold_id}_Scaffold.fasta", scaffold_id = SCAFFOLDS)
    output:
        "scaffolding/parsing_scaffoldBuilder_output/{scaffold_id}_sequence_class.tsv"
    params:
        prefix = "scaffolding/scaffoldBuilder/{scaffold_id}",
        outputDir = "scaffolding/",
        script = parsing_scaffoldBuilder_script,
        blastn = batch_blastn
    shell:
        """
        python {params.script} -i {params.prefix} -o {params.outputDir} -b {params.blastn}
        """

rule batch_scaffold_info:
    input:
        expand("scaffolding/parsing_scaffoldBuilder_output/{scaffold_id}_sequence_class.tsv", scaffold_id = SCAFFOLDS)
    output:
        "batch_scaffold_info.tsv"
    params:
        blastn = batch_blastn
    run:
        blastn = pd.read_csv(params.blastn, sep="\t")
        blastn['unique_sequence_id'] = blastn['library_id'] + "_" + blastn['species'].str.replace(" ", "_") + "_" + blastn['qseqid']
        batch_df=pd.DataFrame()
        for f in input:
            df = pd.read_csv(f, sep="\t")
            batch_df = pd.concat([batch_df, df if not df.empty else None])
        combined_df = blastn.merge(batch_df, on=["unique_sequence_id", "species", "genus", "family", "kingdom", "superkingdom"], how="outer", validate="one_to_many")
        def get_class(class_value):
            if pd.isna(class_value):
                return "no_scaffolding_attempt"
            else:
                return class_value
        combined_df['class'] = combined_df['class'].apply(get_class)
        combined_df.drop_duplicates(inplace=True)
        combined_df.to_csv(output[0], sep="\t", index=False)

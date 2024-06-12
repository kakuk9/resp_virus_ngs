'''
Last updated on 17/05/2024
Written by Kirsty Kwok
Features:
1. Determine if there are more than 2 contigs for a species
2. Determine if the virus is segmented
3. Group virus by segmented or non-segmented
4. For non-segmented virus, determine if the virus can be scaffolded
5. Output the best reference for scaffolding
6. Output the scaffold_id_list
'''
import argparse
from Bio import Entrez
import os
import pandas as pd
import numpy as np
pd.options.mode.copy_on_write = True


# Useful functions for this script
def fasta_ambiguous(sequence):
    AMBIGUOUS_COUNT = 0
    for character in sequence:
        if not character.upper() in ["A", "T", "C", "G"]:
            AMBIGUOUS_COUNT += 1
    SEQ_LEN = len(sequence)
    if AMBIGUOUS_COUNT > 0:
        AMBIGUOUS_PERCENT = AMBIGUOUS_COUNT / SEQ_LEN * 100
    else:
        AMBIGUOUS_PERCENT = 0

    AMBIGUOUS_PERCENT = round(AMBIGUOUS_PERCENT, 2)
    return AMBIGUOUS_PERCENT

def fetch_fasta(accession_no, output_seq=False):
    Entrez.email = "ttkwok9@gmail.com"
    Entrez.api_key = "28d0e8086a4e32273db083b2a8182f213908"
    handle = Entrez.efetch(db="nucleotide", id=accession_no, rettype="fasta", retmode="text")
    record = handle.read()
    handle.close()
    lines=record.split("\n")
    header=lines[0].replace(">", "")
    sequence="".join(lines[1:])
    length = len(sequence)
    amb = fasta_ambiguous(sequence)
    if output_seq:
        return header,length,amb,sequence
    else:
        return header,length,amb

def SCAFFOLDING(qlen_sum, longest_qlen, longest_slen):
    if longest_slen - longest_qlen < 300:
        return("no_longContig")
    elif qlen_sum/longest_slen < 0.7:
        return("no_genomeCoverage<70%")
    else:
        return("yes")

def BEST_REF(INPUT_DICT):
    num_sseqid = len(INPUT_DICT)
    # If there is only one sseqid, return it
    if num_sseqid == 1:
        return(list(INPUT_DICT.keys())[0])
    else:
        #For more than 1 sseqid, we need to determine which reference is the best
        max_occurrence_key = max(INPUT_DICT, key=lambda x: INPUT_DICT[x]['occurrence'])
        max_slen_key = max(INPUT_DICT, key=lambda x: INPUT_DICT[x]['slen'])
        max_length_key = max(INPUT_DICT, key=lambda x: INPUT_DICT[x]['max_length'])
        highpident_key = max(INPUT_DICT, key=lambda x: INPUT_DICT[x]['highest_pident'])
        if max_occurrence_key == max_slen_key == max_length_key == highpident_key:
            best_ref = max_occurrence_key
        else:
            pident_threshold = 90
            # If longest reference with 0% ambiguous_percent and pident >= 90, return it
            if INPUT_DICT[max_slen_key]['ambiguous_percent'] == 0 and (INPUT_DICT[max_slen_key]['lowest_pident'] >= pident_threshold):
                best_ref = max_slen_key
            # Elif reference with the highest occurrence has 0% ambiguous_percent and pident >= 90, return it
            elif INPUT_DICT[max_occurrence_key]['ambiguous_percent'] == 0 and (INPUT_DICT[max_occurrence_key]['lowest_pident'] >= pident_threshold):
                best_ref = max_occurrence_key
            # Else return the reference with the lowest ambiguous_percent
            else:
                best_ref = min(INPUT_DICT, key=lambda x: INPUT_DICT[x]['ambiguous_percent'])
            return(best_ref)
        

# Main script
def main(args):
    
    if args.outputPrefix != "":
        output_prefix = args.outputPrefix + "_"
    else:
        output_prefix = ""
        
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    if not os.path.exists(args.outputReference):
        os.makedirs(args.outputReference)
   
    df = pd.read_csv(args.inputblastn, sep='\t')
    blastn = df.query('superkingdom == "Viruses"')
    blastn['sseqid_acc'] = blastn['sseqid'].str.split('|').str[3]

    if len(blastn) == 0:
        print("Empty dataframe/no viral sequences found.")
        exit(0)

    # This per_sample_summary_viral_species is for determining if there are more than 2 contigs for a species and if the virus is segmented
    per_sample_summary_viral_species = pd.merge(blastn.groupby('file_id')['species'].value_counts().reset_index(), 
                                                blastn[['species', 'kingdom', 'family', 'genus', 'segmented']].drop_duplicates(), 
                                                on='species', how='left')
    per_sample_summary_viral_species.to_csv(f'{args.output}/{output_prefix}virusSpeciesSummary.tsv', sep='\t', index=False)
    # This per_sample_summary_viral_family can be used for heatmap
    per_sample_summary_viral_family = blastn.astype({'family': 'category'}).groupby('file_id')['family'].value_counts().reset_index()
    per_sample_summary_viral_family.to_csv(f'{args.output}/{output_prefix}virusFamilySummary.tsv', sep='\t', index=False)

    from collections import defaultdict
    moreThan2Contigs = per_sample_summary_viral_species.query('count > 1').reset_index(drop=True).to_dict(orient='records')
        
    if len(moreThan2Contigs) == 0:
        print("No species with more than 2 contigs.")
        exit(0)
        
    # Iterate over moreThan2Contigs_dict
    for item in moreThan2Contigs:
        # Filter blastn based on file_id and species
        num_ambiguous_seq = 0
        sseqids, slens, pidents, lengths = [],[],[],[]
        ref_summary={}
        filtered_df = blastn[(blastn['file_id'] == item['file_id']) & (blastn['species'] == item['species'])]
        #if segmented/unknown, all sseqid will be outputed, if non-segmented, only sseqid with slen >= slen_cutoff will be outputed.
        segmented = item['segmented']
        # Concatenate qseqid and qlen with semicolon and add to the dictionary
        item['qseqid'] = ";".join([f"{row.qseqid}" for row in filtered_df.itertuples()])
        item['qlen_sum'] = int(sum([row.qlen for row in filtered_df.itertuples()]))
        item['longest_qlen'] = int(max([row.qlen for row in filtered_df.itertuples()]))
        item['longest_slen'] = int(max([row.slen for row in filtered_df.itertuples()]))
        sseqids.extend([f"{row.sseqid_acc}" for row in filtered_df.itertuples()])
        slens.extend([f"{row.slen}" for row in filtered_df.itertuples()])
        pidents.extend([f"{row.pident}" for row in filtered_df.itertuples()])
        lengths.extend([f"{row.length}" for row in filtered_df.itertuples()])
        slen_cutoff = round(max(list(map(int, slens)))*.85)

        for sseqid in set(sseqids):
            index_list = [i for i, n in enumerate(sseqids) if n == sseqid]
            slen=set(slens[i] for i in index_list)
            slen=int("".join(slen))
            if (slen >= slen_cutoff and segmented == "no") or (segmented!="no"):
                max_length,highest_pident=0,0
                occurrence = int(len(index_list))
                highest_pident = float(max([pidents[i] for i in index_list]))
                lowest_pident = float(min([pidents[i] for i in index_list]))
                max_length = int(max([lengths[i] for i in index_list]))
                header, length, ambiguous_percent =fetch_fasta(sseqid)
                if ambiguous_percent >=1 :
                    num_ambiguous_seq += 1
                ref_summary.update({sseqid:{"slen":slen,"occurrence":occurrence,"max_length":max_length,"highest_pident":highest_pident,"lowest_pident":lowest_pident, "ambiguous_percent":ambiguous_percent}})
        if num_ambiguous_seq < len(set(sseqids)) and segmented != "yes":
            ref_summary = {k:v for k,v in ref_summary.items() if v['ambiguous_percent'] < 1}
        #if num_ambiguous_seq
        item['ref_summary']=ref_summary

            
    df_moreThan2Contigs = pd.DataFrame(moreThan2Contigs)
    df_moreThan2Contigs = df_moreThan2Contigs[['file_id','species','segmented','count','kingdom','family','genus','qseqid','qlen_sum','longest_qlen','longest_slen','ref_summary']]
    df_nonSegmented = df_moreThan2Contigs.query('segmented == "no"').reset_index(drop=True)
    df_nonSegmented['scaffolding'] = df_nonSegmented.apply(lambda x: SCAFFOLDING(x['qlen_sum'], x['longest_qlen'], x['longest_slen']), axis=1)
    df_nonSegmented_scaffolding = df_nonSegmented.query('scaffolding == "yes"').reset_index(drop=True)
        
    df_nonSegmented_scaffolding['best_ref'] = df_nonSegmented_scaffolding['ref_summary'].apply(BEST_REF)
    sent_4_scaffolding = df_nonSegmented_scaffolding[['file_id', 'species', 'qseqid', 'best_ref']].copy()
    sent_4_scaffolding['scaffold_id'] = (sent_4_scaffolding['file_id'] + '_SCAFFOLD_' + sent_4_scaffolding['species']).str.replace(" ", "_")
    #sent_4_scaffolding.drop(['file_id', 'species'], axis=1, inplace=True)
    scaffold_id_best_ref = sent_4_scaffolding[['scaffold_id','best_ref']].copy()
    # header,length,amb,sequence
    scaffold_id_best_ref['header'], scaffold_id_best_ref['length'], scaffold_id_best_ref['ambiguous_percent'], scaffold_id_best_ref['sequence'] = zip(*scaffold_id_best_ref['best_ref'].apply(fetch_fasta, output_seq=True))
    for index, row in scaffold_id_best_ref.iterrows():
        with open(f'{args.outputReference}/{row["scaffold_id"]}.fasta', 'w') as f:
            f.write(f">{row['header']}\n{row['sequence']}")
    scaffold_id_best_ref.drop(['sequence'], axis=1, inplace=True)
    scaffold_id_best_ref.to_csv(f'{args.outputReference}/{output_prefix}nonSegmented_scaffold_id_bestRef.tsv', sep='\t', index=False)
    sent_4_scaffolding = dict(zip(sent_4_scaffolding['scaffold_id'], sent_4_scaffolding['qseqid']))
    for scaffold_id, qseqid in sent_4_scaffolding.items(): 
        with open(f'{args.outputReference}/{scaffold_id}.txt', 'w') as f:
            f.write("\n".join(qseqid.split(";")))
    # Output to files
    df_nonSegmented_scaffolding.to_csv(f'{args.outputReference}/{output_prefix}nonSegmented_4scaffolding.tsv', sep='\t', index=False)
    df_moreThan2Contigs.to_csv(f'{args.outputReference}/{output_prefix}moreThan2Contigs.tsv', sep='\t', index=False)
            

# Set up argparse and run the main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise virus data.')
    parser.add_argument('inputblastn', help='input blastn file')
    parser.add_argument('--outputPrefix', '-p', help='output prefix', default='')
    parser.add_argument('--output', '-o', help='output directory', default="summary")
    parser.add_argument('--outputReference', help='output reference directory', default="reference_seq")
    args = parser.parse_args()
    main(args)

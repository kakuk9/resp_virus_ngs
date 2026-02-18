
"""
Written by Kirsty Kwok, last updated on 02/04/2024
Usage: python blastn_filter.py <blastn>
Description: This python filters blastn output.
Steps:
1. Filter based on evalue and pident
2. Get sseqid accession number
3. Add a column to store the library ID
4. Rank the hits based on bitscore, pident, length and evalue
5. Keep the top hit for each query sequence
7. Use pytaxonkit to get the lineage info
8. Check if sscinames/lineage contains any words in the exclude_list
9. Exclude non-viral and hits that contain words in the exclude_list
10. Save the filtered DataFrame to the output file

Output:
    - 
Wishlist:
    - 
Bugs:
    -
"""

import argparse
import os
import pandas as pd
import re
import warnings
warnings.filterwarnings('ignore')

#Preset
HEADER= ['qseqid', 'sseqid', 'sscinames', 'staxid', 'staxids', 'pident', 'length', 'qlen', 'slen', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
COLUMN_TYPE={'qseqid': 'string', 'sseqid': 'string','sscinames': 'string', 'staxid': 'string', 'staxid': 'string',
                'pident': 'float', 'length': 'int', 'qlen': 'int', 'slen': 'int','mismatch': 'int', 'gapopen': 'int',
                'qstart': 'int', 'qend': 'int', 'sstart': 'int', 'send': 'int',
                'evalue' : 'float', 'bitscore': 'int'}


def main(args):

    sample_prefix = os.path.basename(args.InputBlastnResult).split('.')[0]

    global HEADER
    if args.colname != HEADER:
        HEADER=args.colname.split(',')

    # Only load exclude list if exclusion is not skipped
    if not args.skip_exclusion:
        exclude_list = open(args.ExcludeList).read().splitlines()
    else:
        exclude_list = []

    def check_sscinames(sscinames):
        if any(re.search(r'\b{}\b'.format(word), sscinames, re.IGNORECASE) for word in exclude_list):
            return "yes"
        else:
            return "no"

    essential_columns=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
    if not all(col in HEADER for col in essential_columns):
        raise ValueError(f"Column names must include {essential_columns}")

    blastn= pd.read_csv(args.InputBlastnResult, sep='\t', names=HEADER, keep_default_na=False)

    #blastn.astype(COLUMN_TYPE)
    
    original_num_qseqid = blastn['qseqid'].nunique()
    
    #qseqid sseqid sscinames staxid staxids pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
    #0      1      2         3      4       5      6      7    8    9        10      11     12   13     14   15    16

    # Filter the DataFrame based on the specified criteria (evalue and pident)
    if not args.skip_filtering:
        blastn = blastn[(blastn['evalue'] <= args.evalue) & (blastn['pident'] >= args.id)].reset_index(drop=True)
    
    if len(blastn) > 0:
        # Get sseqid accession number
        # This is removed because accession number will be extracted in add_taxon.sh
        #blastn = blastn.assign(acc_num=lambda blastn: blastn['sseqid'].str.split('|').str[3])
        # Here we add a column to store the library ID, didn't use sample ID because it's not unique  (e.g. DNA, RNA)
        blastn['library_id'] = sample_prefix
        #Sort the DataFrame based on the specified criteria
        blastn.sort_values(by=['qseqid', 'evalue', 'bitscore', 'length', 'pident', 'slen' ], ascending=[True, True, False, False, False, False], inplace=True)
        #Keep the top hit(s) for each query sequence
        blastn = blastn.groupby('qseqid').head(args.numTopHit2Keep).reset_index(drop=True)
        # Apply exclusion filtering only if not skipped
        if not args.skip_exclusion and 'sscinames' in blastn.columns:
            blastn['sscinames_exclude']= blastn['sscinames'].apply(check_sscinames)
            blastn = blastn.query('sscinames_exclude == "no"')
            blastn.drop(columns=['sscinames_exclude'], inplace=True)
        # Check if there is any '|' in sseqid
        if any("|" in s for s in blastn['sseqid']):
            blastn['accNum'] = blastn['sseqid'].str.split('|').str[3].str.split('.').str[0]
        else:
            blastn['accNum'] = blastn['sseqid'].str.split('.').str[0]

        blastn.to_csv(args.OuputBlastnResult, sep='\t', index=False)
        filtered_num_qseqid = blastn['qseqid'].nunique()
    else:
        HEADER.extend(['library_id', 'accNum'])
        empty_df = pd.DataFrame(columns=HEADER)
        #empty_df['library_id'] = sample_prefix
        empty_df.to_csv(args.OuputBlastnResult, sep='\t', index=False)
        filtered_num_qseqid = empty_df['qseqid'].nunique()
    
    num_qseqid_removed = original_num_qseqid - filtered_num_qseqid
    #print(f"Number of unique qseqid in input file: {original_num_qseqid}")
    #print(f"Number of unique qseqid removed after filtering: {num_qseqid_removed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script filters blastn results')
    parser.add_argument('InputBlastnResult',
                        help='blastn result prefix, e.g. /home3/user1/sample1.blastn')
    parser.add_argument('OuputBlastnResult',
                        help='blastn result prefix, e.g. /home3/user1/sample1.blastn')
    parser.add_argument('--ExcludeList',
                        help='List of words to be excluded from the sscinames/lineage', 
                        type=str, 
                        default='/home3/2893911k/kirsty_scripts/smk/respiratory_mNGS/exclude.txt')
    parser.add_argument('-n', '--numTopHit2Keep', type=int, help="Number of top hit(s) to be kept, default: 1", default=1)
    parser.add_argument('-e', '--evalue', type=float, help="E-value cut-off. default: 1e-10", default=1e-10)
    parser.add_argument('--colname', type=str, help=f"Column names, comma-separated, default:{HEADER}", default=HEADER)
    parser.add_argument('--id', type=float, help="Percentage identity cut-off. default: 70", default=70)
    parser.add_argument('--skip-filtering', action='store_true', help="Skip e-value and pident filtering, only keep top hits")
    parser.add_argument('--skip-exclusion', action='store_true', help="Skip exclusion list filtering")
    
    args = parser.parse_args()
    main(args)
"""
This script filters the Diamond m8 output file to keep only the top hit for each query sequence.
Usage: python diamond_m8_filter.py input_Diamond_m8 output_Diamond_m8 (optional: --topHit 1 --evalue 1e-10 --pident 60)
1. Filter the DataFrame based on the <=  evalue and >=  pident
2. sort_values(by=['qseqid', 'evalue', 'bitscore', 'length', 'pident' ], ascending=[True, True, False, False, False], inplace=True)
3. Group by qseqid and keep the topHit rows for each group
4. Save the filtered DataFrame to the output file
"""

import argparse
import pandas as pd


def main(args):
    input_m8 = pd.read_csv(args.input_Diamond_m8, 
                           sep='\t',
                           names=['qseqid', 'sseqid', 
                                  'pident', 'length', 'mismatch', 'gapopen', 
                                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    
    # Filter the DataFrame based on the specified criteria
    filtered_m8 = input_m8[(input_m8['evalue'] <= args.evalue) & (input_m8['pident'] >= args.id)].reset_index(drop=True)
    #input_m8.query('evalue <=  args.evalue and pident >= args.pident').reset_index(drop=True)
    
    # Sort the DataFrame based on the specified criteria
    filtered_m8.sort_values(by=['qseqid', 'evalue', 'bitscore', 'length', 'pident' ], ascending=[True, True, False, False, False], inplace=True)
    # Group by qseqid and keep the topHit rows for each group
    filtered_m8 = filtered_m8.groupby('qseqid').head(args.numTopHit2Keep)
    
    # Save the filtered DataFrame to the output file
    filtered_m8.to_csv(args.output_Diamond_m8, sep='\t', header=None, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Description of your script")
    parser.add_argument("input_Diamond_m8", help="Path to input file")
    parser.add_argument("output_Diamond_m8", help="Path to output file")
    parser.add_argument('-n', '--numTopHit2Keep', type=int, help="Number of top hit(s) to be kept, default: 1", default=1)
    parser.add_argument('-e', '--evalue', type=float, help="e-value cut-off", default=1e-10)
    parser.add_argument('--id', type=int, help="Percent identity cut-off", default=60)
    
    args = parser.parse_args()
    main(args)
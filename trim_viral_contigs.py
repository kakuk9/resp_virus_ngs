"""
Last update on 25/09/2024
Written by Kirsty Kwok
Description:
This script takes in contigs.fasta and its sequencing depth.
It will trim the contigs based on the depth information (it masks mid low coverage with N, and trimmed 5' and 3' low coverage region).
Usage: python trim_viral_contigs.py <input_depth> <input_fasta> [-d <min_depth>]
"""

import argparse
import os
import pyfastx
import pandas as pd
import subprocess

class NoContigsBelowThreshold(Exception):
    pass

def main (args):
    try: 
        with open(args.input_depth, 'r') as f:
            rows = [[row.strip().split('\t')[0], row.strip().split('\t')[1], row.strip().split('\t')[2]] 
                    for row in f.readlines() 
                    if int(row.strip().split('\t')[3]) < args.min_depth]       
        print(rows)
        if len(rows) == 0:
            raise NoContigsBelowThreshold

    except OSError:
        print('Error reading depth file')
    except NoContigsBelowThreshold:
        print('No contigs below threshold.')

    else:
        low_quality_contigs = set(row[0] for row in rows)
        combined_depths = []

        for contig in low_quality_contigs:
            poorly_covered_regions = [row for row in rows if row[0] == contig]
            if len(poorly_covered_regions) > 1:
                keep_going = True
                while keep_going:
                    keep_going = False
                    for i, row in enumerate(poorly_covered_regions):
                        if i < len(poorly_covered_regions) - 1:
                            try:
                                # If the end of the first region is the same as the start of the next region, combine them
                                if poorly_covered_regions[i][2] == poorly_covered_regions[i+1][1]:
                                    poorly_covered_regions[i][2] = poorly_covered_regions[i+1][2]
                                    poorly_covered_regions.pop(i+1)
                                    keep_going = True
                                    break
                            except IndexError:
                                pass
            combined_depths.extend(poorly_covered_regions)
            # Sorting the contigs by length (longest first)
            combined_depths = sorted(combined_depths, key=lambda x: int(x[0].split("_")[3]), reverse=True)
            with open('temp.depth', 'w') as f:
                for row in combined_depths:
                    f.write('\t'.join(row) + '\n')

        # Mask low quality regions with "N" in fasta file
        subprocess.run(f"bedtools maskfasta -fi {args.input_fasta} -bed temp.depth -fo temp.fasta", shell=True)
        processed_sequences = []
        # Convert "N" in 5' and 3' region to * and remove them
        for seq in pyfastx.Fasta('temp.fasta', build_index=True):
            sequence = [*str(seq.seq)]
            if seq.name in low_quality_contigs:
                for i in range(len(sequence)):
                    if sequence[i] == 'N':
                        sequence[i] = '*'
                    else:
                        break
                for i in range(len(sequence)-1, 0, -1):
                    if sequence[i] == 'N':
                        sequence[i] = '*'
                    else:
                        break
            name = seq.name
            if "N" in set(sequence) or "*" in set(sequence):
                name += "_trimmed"

            good_contig_length = len(sequence) - sequence.count('*') - sequence.count('N')
            if good_contig_length >= 50:
                processed_sequences.append(f">{name}")
                processed_sequences.append(''.join(sequence).replace('*', ''))
            
        fasta_prefix = os.path.splitext(args.input_fasta)[0].split('.fasta')[0]
        print(fasta_prefix)
        with open('temp2.fasta', 'w') as f:
            f.write('\n'.join(processed_sequences))
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim viral contigs')
    parser.add_argument('input_depth', help='Input depth file')
    parser.add_argument('input_fasta', help='Input fasta file')
    parser.add_argument('-d', '--min_depth', help='Minimum depth to keep, default=10x', default=10, type=int)
    args = parser.parse_args()
    main(args)
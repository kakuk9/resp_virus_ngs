import argparse
import os
import pandas as pd
from Bio import Entrez
from Bio import SeqIO

# Function to calculate percentage of ambiguous bases (non-ATCG)
def fasta_ambiguous(sequence):
    amb_count = sum(1 for base in sequence if base.upper() not in {"A", "T", "C", "G"})
    seq_len = len(sequence)
    return round((amb_count / seq_len * 100) if seq_len > 0 else 0, 2)

# Function to fetch FASTA from NCBI and compute length and ambiguous percent
def fetch_fasta(accession_no, email, api_key, output_seq=False):
    Entrez.email = email
    Entrez.api_key = api_key
    handle = Entrez.efetch(db="nucleotide", id=accession_no, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    genbank_title = record.description
    sequence = str(record.seq)
    length = len(sequence)
    amb_pct = fasta_ambiguous(sequence)
    if output_seq:
        return genbank_title, length, amb_pct, sequence
    return genbank_title, length, amb_pct

# Selects best reference, returning None if no "kept" contigs are found
def get_best_reference(df, email, api_key, amb_threshold, top_n):
    # Prioritize contigs that were kept after trimming
    kept_df = df[df['status'].str.contains("kept", na=False)].copy()

    # If no contigs were "kept", return None to signal that this group should be skipped.
    if kept_df.empty:
        return None

    # Find the row(s) with the longest FINAL (trimmed) contig length
    max_final_len = kept_df['final_length'].max()
    longest_contig_df = kept_df[kept_df['final_length'] == max_final_len]

    # Within the best contig(s), sort by reference length, then bitscore
    candidates = longest_contig_df.sort_values(['slen_b', 'bitscore_b'], ascending=[False, False])
    
    # Check top N candidates for ambiguity
    for _, row in candidates.head(top_n).iterrows():
        acc = row.get('sacc_b') or row.get('subject_acc') or row.get('accession')
        if not acc:
            continue
        try:
            _, _, amb_pct = fetch_fasta(acc, email, api_key)
            if amb_pct < amb_threshold:
                return row
        except Exception as e:
            print(f"Warning: Could not fetch {acc}. Error: {e}. Skipping.")
            continue
            
    # If no suitable candidate is found, return the top one
    if not candidates.empty:
        return candidates.iloc[0]
    else:
        return None

# CORRECTED: Writes a FASTA file with a unified filename, no if/else.
def write_fasta(sequence, header, accession, sample_prefix, genus, segment, output_dir):
    filename = f"{sample_prefix}--{genus}_{segment}.fasta"
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w') as f:
        f.write(f">{accession} {header}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')
    return filepath

# CORRECTED: Creates an empty FASTA file with a unified filename, no if/else.
def write_empty_fasta(sample_prefix, genus, segment, output_dir):
    filename = f"{sample_prefix}--{genus}_{segment}.fasta"
    filepath = os.path.join(output_dir, filename)
    open(filepath, 'w').close()
    return filepath

def main(args):
    try:
        blastn = pd.read_csv(args.input_merged_blastn, sep='\t', header=0)
        trim_summary = pd.read_csv(args.input_trim_summary, sep=r'\s+', header=0)
    except FileNotFoundError as e:
        raise ValueError(f"Input file not found: {e.filename}")
    if blastn.empty:
        open(args.output_file, 'w').close()
        exit(0)

    merged_df = pd.merge(blastn, trim_summary, left_on='qseqid', right_on='contig_name', how='left')
    merged_df['status'].fillna('kept', inplace=True)
    if 'final_length' not in merged_df.columns:
        merged_df['final_length'] = merged_df['qlen_b']
    else:
        merged_df['final_length'].fillna(merged_df['qlen_b'], inplace=True)

    def get_segment_key(row):
        is_segmented = (row.get('segmented_b') == 'yes') or (row.get('segmented_r') == 'yes')
        if is_segmented:
            segment_name = row.get('segment_name')
            if pd.notna(segment_name) and segment_name != 'N/A':
                return str(segment_name).replace(" ", "-")
            else:
                return "segmented-unknown"
        else:
            return "non-segmented"

    merged_df['segment_key'] = merged_df.apply(get_segment_key, axis=1)
    group_cols = ['genus_b', 'segment_key']
    os.makedirs(args.fasta_output_dir, exist_ok=True)

    results = []
    for group_keys, group_df in merged_df.groupby(group_cols):
        best = get_best_reference(group_df, args.entrez_email, args.entrez_api_key, args.amb_threshold, args.top_n)
        
        genus = group_keys[0]
        segment = group_keys[1]
        
        if best is None:
            print(f"Skipping group: {group_keys} -> All contigs were discarded or no valid reference found.")
            
            # Still create the empty FASTA file as a placeholder on the filesystem
            write_empty_fasta(args.sample_prefix, genus, segment, args.fasta_output_dir)
            
            # Add a placeholder row to the output TSV with the requested string
            results.append({
                'genus_b': genus, 
                'segment_key': segment,
                'qseqid': 'N/A',
                'sacc_b': 'N/A',
                'fasta_file': 'no_valid_reference_all_contigs_discarded',
                'status': 'all_contigs_discarded'
            })
            continue # Move to the next group

        acc = best.get('accNum_b') or best.get('sacc_b')
        print(f"Processing group: {group_keys} -> Best ref: {acc}")
        
        header, length, amb_pct, seq = fetch_fasta(acc, args.entrez_email, args.entrez_api_key, output_seq=True)
        fasta_path = write_fasta(seq, header, acc, args.sample_prefix, genus, segment, args.fasta_output_dir)
        
        best_info = best.to_dict()
        best_info.update({'fasta_file': fasta_path, 'ambiguous_pct': amb_pct})
        results.append(best_info)

    if results:
        best_df = pd.DataFrame(results)
        best_df.to_csv(args.output_file, sep='\t', index=False)
    else:
        open(args.output_file, 'w').close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select and fetch best reference per genus (and segment) from merged BLASTN results.")
    parser.add_argument("input_merged_blastn", type=str,help="Path to the input merged BLASTN TSV file.")
    parser.add_argument("input_trim_summary", type=str,help="Path to the input trimming summary TSV file.")
    parser.add_argument("-o", "--output_file", type=str, required=True,help="Path to write the selected best references metadata TSV.")
    parser.add_argument("--entrez_email", type=str, default = "k.kwok.1@research.gla.ac.uk", help="Email for NCBI Entrez.")
    parser.add_argument("--entrez_api_key", type=str, default = "28d0e8086a4e32273db083b2a8182f213908", help="API key for NCBI Entrez.")
    parser.add_argument("--amb_threshold", type=float, default=10.0,help="Maximum allowed percentage of ambiguous bases (default: 10%%).")
    parser.add_argument("--top_n", type=int, default=3,help="Number of top candidates to check for ambiguity (default: 3).")
    parser.add_argument("-p", "--sample_prefix", type=str, required=True,help="Prefix for naming output FASTA files (e.g., sample name).")
    parser.add_argument("--fasta_output_dir", type=str, required=True,help="Directory to write output FASTA files.")
    args = parser.parse_args()
    main(args)
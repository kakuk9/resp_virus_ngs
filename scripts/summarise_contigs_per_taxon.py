"""
Last update on 30/09/2025
Written by Kirsty Kwok
Description:
This script takes in contigs, their BLAST alignments, and read mapping data 
to produce a summary of viral statistics per taxonomic group and segment.

It now integrates contig trimming information and processes segmented viruses
by outputting one row per segment, with separate columns for total genome size
and individual segment size.
"""

import argparse
import subprocess
import pandas as pd
import pyfastx


def get_total_read_count(fastq_files: str) -> int:
    """Calculates the total number of reads from one or more FASTQ files."""
    try:
        return sum(len(pyfastx.Fastq(fq.strip())) for fq in fastq_files.split(','))
    except Exception:
        return 0

def parse_idxstats(idxstats_file: str) -> tuple:
    """
    Parses samtools idxstats file and returns:
    - total number of reads (sum of 3rd and 4th columns)
    - total number of viral reads (sum of 3rd column)
    - dictionary mapping contig name to read count (3rd column value)
    """
    total_reads = 0
    viral_reads = 0
    contig_read_counts = {}

    try:
        with open(idxstats_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 4:
                    continue

                contig_name = parts[0]
                mapped_reads = int(parts[2])
                unmapped_reads = int(parts[3])

                # Skip the special '*' line that contains unmapped reads
                if contig_name == '*':
                    continue

                total_reads += mapped_reads + unmapped_reads
                viral_reads += mapped_reads
                contig_read_counts[contig_name] = mapped_reads

        return total_reads, viral_reads, contig_read_counts
    except (FileNotFoundError, ValueError, IndexError):
        return 0, 0, {}

def run_subprocess(command: str) -> str:
    """Runs a shell command and returns its stripped stdout."""
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""

### UPDATED: Added your original function back to calculate TOTAL genome size.
def get_total_genome_size(col_name: str, uniq_value: str, refseq_tsv_path: str) -> float:
    """Estimates the median TOTAL genome size for a taxonomic group, summing segments if present."""
    try:
        refseq_df = pd.read_csv(refseq_tsv_path, sep='\t')
        col_name_cleaned = col_name.replace("_b", "")

        df_filtered = refseq_df[refseq_df[col_name_cleaned] == uniq_value].copy()

        if df_filtered.empty:
            return 0.0

        # Check if this is a segmented virus (has non-NA segments with multiple unique segments)
        has_segments = df_filtered['Segment'].notna().any()

        if has_segments:
            # Filter to only segmented entries
            df_segmented = df_filtered.dropna(subset=['Segment'])
            if not df_segmented.empty and df_segmented['Segment'].nunique() > 1:
                # Multi-segmented virus: sum median length of each segment
                return df_segmented.groupby('Segment')['Length'].median().sum()

        # Non-segmented virus or single segment: return median length of all entries
        return df_filtered['Length'].median()

    except (FileNotFoundError, KeyError):
        return 0.0

### UPDATED: This function is now dedicated to getting only an INDIVIDUAL segment's size.
def get_individual_segment_size(col_name: str, uniq_value: str, segment_name: str, refseq_tsv_path: str) -> float:
    """Estimates the median size for a SINGLE, specific segment of a taxonomic group."""
    # This function should not be called for non-segmented viruses.
    if segment_name == "non-segmented":
        return 0.0
    try:
        refseq_df = pd.read_csv(refseq_tsv_path, sep='\t')
        col_name_cleaned = col_name.replace("_b", "")
        
        df_filtered = refseq_df[refseq_df[col_name_cleaned] == uniq_value].copy()
        
        # Unsanitize the name to match the RefSeq file (e.g., "Segment-1" -> "Segment 1")
        original_segment_name = segment_name.replace('-', ' ').replace('_', ' ')
        segment_df = df_filtered[df_filtered['Segment'] == original_segment_name]
        
        return segment_df['Length'].median() if not segment_df.empty else 0.0

    except (FileNotFoundError, KeyError):
        return 0.0

def main(args):
    """Main execution function."""
    total_read_count = get_total_read_count(args.fq)
    _, total_viral_read_count, contig_read_counts = parse_idxstats(args.input_idxstats)

    trim_df = pd.DataFrame()
    if args.trimming_summary:
        try:
            trim_df = pd.read_csv(args.trimming_summary, sep=r'\s+', index_col='contig_name')
        except (FileNotFoundError, KeyError):
            pass # Fail silently if file/column is missing

    try:
        blast_df = pd.read_csv(args.input_merged_refseq_blastn, sep='\t')
        if blast_df.empty: blast_df = pd.DataFrame()
    except FileNotFoundError:
        blast_df = pd.DataFrame()

    uniq_groups = []
    if not blast_df.empty:
        def get_segment_key(row):
            is_segmented = (('segmented_r' in row.index and row['segmented_r'] == 'yes') or 
                           ('segmented_b' in row.index and row['segmented_b'] == 'yes'))

            if is_segmented and 'segment_name' in row.index and pd.notna(row['segment_name']):
                return str(row['segment_name']).replace(' ', '-').replace('_', '-')
            else:
                return "non-segmented"

        blast_df['segment_key'] = blast_df.apply(get_segment_key, axis=1)
        uniq_groups = blast_df[[args.column_name, 'segment_key']].drop_duplicates().values.tolist()

    if not uniq_groups:
        ### UPDATED: Added new columns to the "no results" output.
        df = pd.DataFrame({
            "library_id": [args.sample_id],
            "total_read_count": [total_read_count], "mapped_read_count": [0], "percent_mapped_total_read": [0],
            "total_viral_read_count": [total_viral_read_count], "rpkm": [0],
            "column_name": [args.column_name], "uniq_value": ["No_value_found"], "segment_name": ["N/A"],
            "unique_species_b": ["N/A"], "unique_species_r": ["N/A"], "unique_sscinames": ["N/A"],
            "expected_genome_size": [0], "expected_segment_size": [0],
            "contig_count": [0], "max_depth": [0], "min_depth": [0],
            "longest_contig_len": ["No_contig_found"], "shortest_contig_len": ["No_contig_found"],
            "contig_ids": ["No_contig_found"], "discarded_contig_ids": ["No_contig_found"]
        })
        df.to_csv(args.output_file, index=False, sep='\t')
        return

    rows = []
    for group, segment in uniq_groups:
        mask = (blast_df[args.column_name] == group) & (blast_df['segment_key'] == segment)
        all_contig_ids = blast_df.loc[mask, 'qseqid'].unique().tolist()

        # Collect unique species from _b and _r columns separately
        species_b_values = blast_df.loc[mask, 'species_b'].dropna().unique() if 'species_b' in blast_df.columns else []
        species_r_values = blast_df.loc[mask, 'species_r'].dropna().unique() if 'species_r' in blast_df.columns else []
        all_species_b = sorted(species_b_values)
        all_species_r = sorted(species_r_values)

        # Collect unique sscinames (column has no suffix)
        sscinames_values = blast_df.loc[mask, 'sscinames'].dropna().unique() if 'sscinames' in blast_df.columns else []
        all_sscinames = sorted(sscinames_values)

        kept_contig_ids = []
        discarded_contig_ids = []
        if not trim_df.empty:
            for contig_id in all_contig_ids:
                if contig_id in trim_df.index and "discarded" in trim_df.loc[contig_id, 'status']:
                    discarded_contig_ids.append(contig_id)
                else:
                    kept_contig_ids.append(contig_id)
        else:
            kept_contig_ids = all_contig_ids
        
        ### UPDATED: Calculate both total genome size and individual segment size.
        total_genome_size = get_total_genome_size(args.column_name.capitalize(), group, args.refseq_tsv)
        segment_size = get_individual_segment_size(args.column_name.capitalize(), group, segment, args.refseq_tsv)

        contig_lengths = [int(cid.split('_')[3]) for cid in kept_contig_ids]

        # Count reads from idxstats for kept contigs
        total_mapped_read_count = 0
        for contig in kept_contig_ids:
            total_mapped_read_count += contig_read_counts.get(contig, 0)

        # Parse depth file to get min/max depths
        with open(args.input_depth, 'r') as depth_file:
            depths = {int(float(line.strip().split('\t')[3])) for line in depth_file if line.strip().split('\t')[0] in kept_contig_ids}

        # Calculate simple percentage of total reads (not genome-size adjusted)
        percent_total_read = (total_mapped_read_count / total_read_count * 100) if total_read_count > 0 else 0

        # Calculate RPKM (Reads Per Kilobase per Million viral reads)
        # RPKM = (mapped_reads / (genome_size_kb)) / (viral_reads / 1_000_000)
        rpkm = (total_mapped_read_count / (total_genome_size / 1000)) / (total_viral_read_count / 1000000) if total_viral_read_count > 0 and total_genome_size > 0 else 0

        ### UPDATED: Added the new columns to the final output row.
        row_data = {
            "library_id": args.sample_id,
            "total_read_count": total_read_count,
            "mapped_read_count": total_mapped_read_count,
            "percent_mapped_total_read": percent_total_read,
            "total_viral_read_count": total_viral_read_count,
            "rpkm": rpkm,
            "column_name": args.column_name,
            "uniq_value": group,
            "segment_name": segment,
            "unique_species_b": "|".join(all_species_b),
            "unique_species_r": "|".join(all_species_r),
            "unique_sscinames": "|".join(all_sscinames),
            "expected_genome_size": total_genome_size,
            "expected_segment_size": segment_size,
            "contig_count": len(kept_contig_ids),
            "max_depth": max(depths) if depths else 0,
            "min_depth": min(depths) if depths else 0,
            "longest_contig_len": max(contig_lengths) if contig_lengths else 0,
            "shortest_contig_len": min(contig_lengths) if contig_lengths else 0,
            "contig_ids": "|".join(sorted(kept_contig_ids)),
            "discarded_contig_ids": "|".join(sorted(discarded_contig_ids))
        }
        rows.append(row_data)

    df = pd.DataFrame(rows)
    df.to_csv(args.output_file, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Summarize viral statistics from mapping and alignment data, with segment-aware reporting.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('input_merged_refseq_blastn', help='Input merged blastn file (TSV)')
    parser.add_argument('input_depth', help='Input depth file from final mapping')
    parser.add_argument('input_idxstats', help='Input idxstats file from final mapping')
    parser.add_argument('--trimming_summary', '-t', help='(Optional) TSV summary file from the contig trimming step.')
    parser.add_argument('--output_file', '-o', help='Output file name (TSV)', required=True)
    parser.add_argument('--sample_id', '-p', help='Sample ID', required=True)
    parser.add_argument('--column_name', '-c', help='Column name to group by (e.g., genus_b, family_b)', required=True)
    parser.add_argument('--fq', help = "Comma-separated list of input host-depleted fastq files.", required=True)
    parser.add_argument('--refseq_tsv', '-r', 
                        default="/home3/2893911k/kirsty_db/refseq_viral/new_complete/RefSeq_viral_20241023.tsv",
                        help="Path to the RefSeq viral metadata TSV.")
    
    args = parser.parse_args()
    main(args)
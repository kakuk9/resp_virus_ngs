"""
Last updated: 2025-10-14
This script is used to scaffold contigs using nucmer coords with flexible gap filling options.
The script stitches contigs together and fills uncovered regions with either:
  - N's (default, traditional scaffolding)
  - Reference sequence (creates hybrid reference for mapping/consensus)
  - Both (generates both output types)
The script will output the final scaffold(s), summary files, and alignment of overlapping regions.
Output files are named with _filln or _fillr suffix to indicate fill method.
Usage: python scaffold_nucmer_coords_stitch_with_ref.py <best_reference_fasta> <contig_fasta> [--fill {n,r,b}] [--coverage_threshold 90] [--output_prefix scaffold]
"""
import argparse
import pyfastx
import pandas as pd
import os
import subprocess
import sys
import collections

from nucmer_functions import downsize_contig_pool_range
from nucmer_functions import getSubseq_aln


def generate_scaffold(keep_contig_dict, fasta, contig_to_keep_fasta_obj, ref_seq, LEN_R, fill_method):
    """
    Generate scaffold with specified fill method.

    Args:
        keep_contig_dict: Dictionary of contigs to use for scaffolding
        fasta: Contig fasta object
        contig_to_keep_fasta_obj: Filtered contig fasta object
        ref_seq: Reference sequence (None if fill_method='n')
        LEN_R: Reference length
        fill_method: 'n' for N's, 'r' for reference sequence

    Returns:
        tuple: (final_seq, aln_check, overlapping_contig)
    """
    final_seq = ""
    aln_check = ""
    overlapping_contig = []

    for i in range(len(keep_contig_dict)):
        before_seq, seq, after_seq, overlap_cons = "", "", "", ""
        current_query = list(keep_contig_dict.keys())[i]
        gap = (list(keep_contig_dict.values())[i]["GAP2NEXT"])
        i_start = list(keep_contig_dict.values())[i]["S2"]
        i_S1 = list(keep_contig_dict.values())[i]["S1"]
        i_end = list(keep_contig_dict.values())[i]["E2"]
        i_ort = list(keep_contig_dict.values())[i]["ORT"]

        # sequence for current contig
        if i_ort == "+":
            i_seq = str(fasta[current_query])[i_start-1:i_end]
        elif i_ort == "-":
            i_seq = str(fasta[current_query].antisense)[i_end-1:i_start]
        else:
            i_seq=""
        start = 0
        end = len(i_seq)

        if i > 0:
            p = i - 1
            p_gap = (list(keep_contig_dict.values())[p]["GAP2NEXT"])
            if p_gap < 0:
                start = abs(p_gap)
        # alignment for overlapping regions
        if not i == len(keep_contig_dict)-1:
            n = i + 1
            next_query = list(keep_contig_dict.keys())[n]
            if gap < 0:
                n_start = list(keep_contig_dict.values())[n]["S2"]
                n_end = list(keep_contig_dict.values())[n]["E2"]
                n_ort = list(keep_contig_dict.values())[n]["ORT"]
                aln,overlap_cons = getSubseq_aln(current_query, next_query, i_start, i_end, i_ort, n_start, n_end, n_ort, gap, contig_to_keep_fasta_obj)
                overlapping_contig.extend([current_query, next_query])
                aln_check+=aln
        # final scaffold - fill uncovered regions with N's or reference sequence
        if i_S1 > 1 and i == 0:
            if fill_method == 'r':
                # Use reference sequence before first contig
                ref_start_pos = 0
                ref_end_pos = i_S1 - 1
                before_seq = ref_seq[ref_start_pos:ref_end_pos]
            else:
                # Use N's before first contig
                before_seq = "N" * (i_S1 - 1)
        if gap >= 0:
            if fill_method == 'r':
                # Use reference sequence for gap between contigs
                # Get reference coordinates for the gap
                i_E1 = list(keep_contig_dict.values())[i]["E1"]  # End position on reference for current contig
                if i < len(keep_contig_dict) - 1:
                    n_S1 = list(keep_contig_dict.values())[i+1]["S1"]  # Start position on reference for next contig
                    # Extract reference sequence for the gap region
                    after_seq = ref_seq[i_E1:n_S1-1]
                else:
                    # Last contig - no gap sequence needed yet
                    after_seq = ""
            else:
                # Use N's for gaps
                after_seq = "N" * gap
            seq = i_seq
        elif gap < 0:
            end = gap

        seq = i_seq[start:end]
        final_seq = final_seq + before_seq + seq + after_seq + overlap_cons

    # Add fill sequence after the last contig if it doesn't reach the end
    last_contig_E1 = list(keep_contig_dict.values())[-1]["E1"]
    if last_contig_E1 < LEN_R:
        if fill_method == 'r':
            # Use reference sequence from end of last contig to end of reference
            final_seq = final_seq + ref_seq[last_contig_E1:LEN_R]
        else:
            # Use N's from end of last contig to end of reference
            final_seq = final_seq + "N" * (LEN_R - last_contig_E1)

    return final_seq, aln_check, overlapping_contig


def main(args):
    # Set up output file names
    output_prefix = args.output_prefix
    delta_file = f"{output_prefix}.delta"
    coords_file = f"{output_prefix}.coords"
    contigs_to_keep_fasta = f"{output_prefix}_contigs_to_keep.fasta"
    tiling_file = f"{output_prefix}_tiling.tsv"
    aln_check_file = f"{output_prefix}_aln_check.txt"

    # Remove existing fasta index files to avoid conflicts
    for index_file in [f"{args.contig_fasta}.fxi", f"{args.best_reference_fasta}.fxi", f"{contigs_to_keep_fasta}.fxi"]:
        if os.path.exists(index_file):
            os.remove(index_file)

    # Run nucmer alignment
    print(f"Running nucmer alignment...")
    try:
        result = subprocess.run(
            f'nucmer --prefix={output_prefix} {args.best_reference_fasta} {args.contig_fasta}',
            shell=True, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running nucmer: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    # Generate coords file
    print(f"Generating coords file...")
    try:
        result = subprocess.run(
            f'show-coords -rcl {delta_file} > {coords_file}',
            shell=True, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running show-coords: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    # Load contig fasta file
    print(f"Loading contig fasta file...")
    try:
        fasta = pyfastx.Fasta(args.contig_fasta, build_index=True)
    except Exception as e:
        print(f"Error loading contig fasta: {e}", file=sys.stderr)
        sys.exit(1)

    # Load reference fasta file if needed for reference filling
    ref_seq = None
    if args.fill in ['r', 'b']:
        print(f"Loading reference fasta file for gap filling...")
        try:
            ref_fasta = pyfastx.Fasta(args.best_reference_fasta, build_index=True)
            # Get the first (and typically only) reference sequence
            ref_name = list(ref_fasta.keys())[0]
            ref_seq = str(ref_fasta[ref_name])
        except Exception as e:
            print(f"Error loading reference fasta: {e}", file=sys.stderr)
            sys.exit(1)

    # Load and filter nucmer coords file
    print(f"Filtering alignments by coverage threshold ({args.coverage_threshold}%)...")
    try:
        with open(coords_file, "r") as f:
            raw_coords = [
                x.strip().replace("|", "").split()
                for x in f.readlines()[5:]
                if len(x.strip().replace("|", "").split()) > 0 and
                   float(x.strip().replace("|", "").split()[-3]) >= args.coverage_threshold
            ]
    except FileNotFoundError:
        print(f"Error: Coords file {coords_file} not found", file=sys.stderr)
        sys.exit(1)
    except (ValueError, IndexError) as e:
        print(f"Error parsing coords file: {e}", file=sys.stderr)
        sys.exit(1)

    if not raw_coords:
        print(f"Warning: No alignments passed the coverage threshold of {args.coverage_threshold}%", file=sys.stderr)
        print(f"Creating empty output files for failed scaffolding...")

        # Create empty/minimal output files so workflow can continue
        # Determine output file paths
        if args.output_scaffold:
            final_scaffold_file = args.output_scaffold
        else:
            fill_suffix = f"_fill{args.fill}" if args.fill != 'b' else "_fillr"
            final_scaffold_file = f"{output_prefix}_final_scaffold{fill_suffix}.fasta"

        if args.output_tiling:
            tiling_file_out = args.output_tiling
        else:
            tiling_file_out = tiling_file

        if args.output_summary:
            summary_file_out = args.output_summary
        else:
            fill_suffix = f"_fill{args.fill}" if args.fill != 'b' else "_fillr"
            summary_file_out = f"{output_prefix}_summary{fill_suffix}.txt"

        # Write empty scaffold (just copy the input contigs as-is)
        with open(final_scaffold_file, 'w') as f_out, open(args.contig_fasta, 'r') as f_in:
            f_out.write(f_in.read())

        # Write empty tiling file
        with open(tiling_file_out, 'w') as f:
            f.write("CONTIG\tS1\tE1\tS2\tE2\tLEN1\tLEN2\t%COV\t%IDY\tORT\tGAP2NEXT\n")

        # Write summary indicating failure
        with open(summary_file_out, 'w') as f:
            f.write(f"Scaffolding failed: No alignments passed coverage threshold of {args.coverage_threshold}%\n")
            f.write(f"Input contigs copied to output without scaffolding\n")

        print(f"Empty outputs written. Exiting with success status.")
        sys.exit(0)

    coords = downsize_contig_pool_range(raw_coords)
    print(f"Found {len(coords)} contig alignments to scaffold")

    # Convert the processed 'coords' list into a dictionary for scaffolding
    print(f"Building contig dictionary for scaffolding...")
    keep_contig_dict = {}
    LEN_R = None  # Reference length - should be same for all alignments

    for item in coords:
        contig_name = item[12]  # The contig name is the 13th element
        S2 = int(item[2])
        E2 = int(item[3])

        # Get reference length from first alignment
        if LEN_R is None:
            LEN_R = int(item[7])

        if S2 < E2:
            orientation = "+"
        else:
            orientation = "-"

        keep_contig_dict[contig_name] = {
            "S1": int(item[0]),
            "E1": int(item[1]),
            "S2": S2,
            "E2": E2,
            "ORT": orientation,
            "GAP2NEXT": int(item[-1]),
        }

    # Write contigs to keep to a temporary fasta file
    print(f"Writing contigs to keep...")
    contig_to_keep = ""
    for x in keep_contig_dict:
        contig_to_keep += f">{x}\n{fasta[x]}\n"

    try:
        with open(contigs_to_keep_fasta, 'w') as f:
            f.write(contig_to_keep)
        contig_to_keep_fasta_obj = pyfastx.Fasta(contigs_to_keep_fasta, build_index=True)
    except Exception as e:
        print(f"Error writing contigs to keep: {e}", file=sys.stderr)
        sys.exit(1)

    # Determine which fill methods to use
    fill_methods = []
    if args.fill == 'b':
        fill_methods = ['n', 'r']
        print(f"Generating both scaffolds (N's and reference filling)...")
    else:
        fill_methods = [args.fill]
        fill_desc = "N's" if args.fill == 'n' else "reference sequence"
        print(f"Scaffolding contigs with {fill_desc}...")

    # Generate scaffold(s) for each fill method
    scaffolds = {}
    for fill_method in fill_methods:
        final_seq, aln_check, overlapping_contig = generate_scaffold(
            keep_contig_dict, fasta, contig_to_keep_fasta_obj, ref_seq, LEN_R, fill_method
        )
        scaffolds[fill_method] = {
            'seq': final_seq,
            'aln_check': aln_check,
            'overlapping_contig': overlapping_contig
        }

    # Output files
    print(f"Writing output files...")

    # Determine tiling file path (use explicit path if provided, otherwise use prefix)
    if args.output_tiling:
        tiling_file_final = args.output_tiling
    else:
        tiling_file_final = tiling_file

    # Write tiling file (shared for all fill methods)
    try:
        pd.DataFrame.from_dict(keep_contig_dict, orient='index').reset_index().rename(
            columns={"index": "CONTIG"}
        ).to_csv(tiling_file_final, sep="\t", index=False)
    except Exception as e:
        print(f"Error writing tiling file: {e}", file=sys.stderr)

    # Write scaffold files for each fill method
    for fill_method, scaffold_data in scaffolds.items():
        final_seq = scaffold_data['seq']
        aln_check = scaffold_data['aln_check']
        overlapping_contig = scaffold_data['overlapping_contig']

        # Set up output file names with fill method suffix
        fill_suffix = f"_fill{fill_method}"

        # Use explicit paths if provided, otherwise use prefix-based naming
        if args.output_scaffold and args.fill != 'b':
            # If explicit scaffold path provided and not generating both, use it directly
            final_scaffold_file = args.output_scaffold
        else:
            final_scaffold_file = f"{output_prefix}_final_scaffold{fill_suffix}.fasta"

        if args.output_summary and args.fill != 'b':
            # If explicit summary path provided and not generating both, use it directly
            summary_file = args.output_summary
        else:
            summary_file = f"{output_prefix}_summary{fill_suffix}.txt"

        aln_check_file_method = f"{output_prefix}_aln_check{fill_suffix}.txt"

        # Write final scaffold
        try:
            with open(final_scaffold_file, 'w') as f:
                f.write(f">final_scaffold\n{final_seq}\n")
        except Exception as e:
            print(f"Error writing final scaffold ({fill_method}): {e}", file=sys.stderr)
            sys.exit(1)

        # Write alignment check
        try:
            with open(aln_check_file_method, 'w') as f:
                f.write(aln_check)
        except Exception as e:
            print(f"Error writing alignment check ({fill_method}): {e}", file=sys.stderr)

        # Calculate composition
        comp = collections.defaultdict(int)
        for nt in final_seq:
            comp[nt] += 1

        # Collect positions with N(s)
        N_pos = set(i+1 for i in range(len(final_seq)) if final_seq[i] == 'N')

        # Collect N(s) position range
        if N_pos:
            N_pos_range = [(a,b)
                for l in [sorted([i for i in N_pos if {i-1,i+1} - N_pos])]
                for a,b in zip(l[::2],l[1::2]+[l[-1]]) if l]
        else:
            N_pos_range = []

        # Calculate summary statistics
        # For hybrid reference: calculate contig coverage vs reference fill
        total_contig_length = sum(abs(v["E2"] - v["S2"]) + 1 for v in keep_contig_dict.values())
        contig_coverage = round(total_contig_length / len(final_seq) * 100, 2) if len(final_seq) > 0 else 0
        gap_fill = round((len(final_seq) - total_contig_length) / len(final_seq) * 100, 2) if len(final_seq) > 0 else 0

        genome_coverage = round((len(final_seq)-comp["N"])/len(final_seq)*100, 2) if len(final_seq) > 0 else 0
        gc_content = round((comp["G"]+comp["C"])/len(final_seq)*100, 2) if len(final_seq) > 0 else 0
        n_content = round((comp["N"])/len(final_seq)*100, 2) if len(final_seq) > 0 else 0
        overlapping_contig = list(dict.fromkeys(overlapping_contig))

        # Write summary file
        fill_method_desc = "Reference sequence" if fill_method == 'r' else "N's"
        try:
            with open(summary_file, 'w') as f:
                f.write(f"Scaffold length: {len(final_seq)}\n")
                f.write(f"Fill method: {fill_method_desc}\n")
                f.write(f"Reference sequence: {args.best_reference_fasta}\n")
                f.write(f"Reference length: {LEN_R}\n")
                f.write(f"Contig coverage (% from contigs): {contig_coverage}%\n")
                if fill_method == 'r':
                    f.write(f"Reference fill (% from reference): {gap_fill}%\n")
                else:
                    f.write(f"Gap fill (% with N's): {gap_fill}%\n")
                f.write(f"Genome coverage: {genome_coverage}%\n")
                f.write(f"N content: {n_content}%\n")
                f.write(f"GC content: {gc_content}%\n")
                f.write(f"Composition: {dict(comp)}\n")
                f.write(f"Number of contig(s) used for scaffolding: {len(keep_contig_dict)}\n")
                f.write(f"Number of contig(s) overlapping: {len(overlapping_contig)}\n")
                f.write(f"Contig(s) used for scaffolding:\n{';'.join(list(keep_contig_dict.keys()))}\n")
                f.write(f"Contig(s) overlapping:\n{';'.join(overlapping_contig)}\n")
                if N_pos_range:
                    f.write("Ns position:\n{}\n".format("\n".join(f"{x}-{y}" for x, y in N_pos_range)))
                else:
                    f.write("Ns position: None\n")
            print(f"  [{fill_method_desc}] Scaffold: {final_scaffold_file}")
            print(f"  [{fill_method_desc}] Summary: {summary_file}")
            print(f"  [{fill_method_desc}] Contigs cover {contig_coverage}%, gaps filled with {gap_fill}%")
        except Exception as e:
            print(f"Error writing summary file ({fill_method}): {e}", file=sys.stderr)
            sys.exit(1)

    print(f"\nScaffolding complete!")
    print(f"Shared files written: {tiling_file}, {aln_check_file}")
    if args.fill == 'b':
        print(f"Generated both N-filled and reference-filled scaffolds")
    print(f"All output files written with prefix: {output_prefix}_*")

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Scaffold contigs using nucmer alignment with flexible gap filling options',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Traditional scaffolding with N's (default)
  python scaffold_nucmer_coords_stitch_with_ref.py reference.fasta contigs.fasta
  # Output: scaffold_final_scaffold_filln.fasta

  # Hybrid reference with reference sequence filling gaps
  python scaffold_nucmer_coords_stitch_with_ref.py reference.fasta contigs.fasta --fill r
  # Output: scaffold_final_scaffold_fillr.fasta

  # Generate both outputs
  python scaffold_nucmer_coords_stitch_with_ref.py reference.fasta contigs.fasta --fill b
  # Output: scaffold_final_scaffold_filln.fasta AND scaffold_final_scaffold_fillr.fasta

  # Custom settings
  python scaffold_nucmer_coords_stitch_with_ref.py reference.fasta contigs.fasta -f b -ct 95 -p my_scaffold

Fill options:
  -f n : Fill gaps with N's (traditional scaffolding, default)
  -f r : Fill gaps with reference sequence (creates hybrid reference for mapping/consensus)
  -f b : Generate both N-filled and reference-filled scaffolds

Output files:
  <prefix>_final_scaffold_filln.fasta  - Scaffold with N's filling gaps
  <prefix>_final_scaffold_fillr.fasta  - Scaffold with reference filling gaps
  <prefix>_summary_filln.txt           - Summary for N-filled scaffold
  <prefix>_summary_fillr.txt           - Summary for reference-filled scaffold
  <prefix>_tiling.tsv                  - Contig tiling information (shared)
  <prefix>_aln_check_fill*.txt         - Alignment check for overlaps
        """
    )
    parser.add_argument('best_reference_fasta', help='Path to the reference genome fasta file')
    parser.add_argument('contig_fasta', help='Path to the contigs fasta file')
    parser.add_argument('--coverage_threshold', '-ct', type=float, default=90.0,
                        help='Minimum coverage threshold for nucmer alignments (default: 90.0)')
    parser.add_argument('--output_prefix', '-p', type=str, default='scaffold',
                        help='Prefix for all output files (default: scaffold)')
    parser.add_argument('--fill', '-f', type=str, choices=['n', 'r', 'b'], default='n',
                        help='Fill method for uncovered regions: n=N\'s (default), r=reference sequence, b=both (generates both outputs)')

    # Add options for explicit output paths (overrides prefix-based naming)
    parser.add_argument('--output_scaffold', type=str, default=None,
                        help='Explicit path for output scaffold file (overrides prefix-based naming)')
    parser.add_argument('--output_tiling', type=str, default=None,
                        help='Explicit path for output tiling file (overrides prefix-based naming)')
    parser.add_argument('--output_summary', type=str, default=None,
                        help='Explicit path for output summary file (overrides prefix-based naming)')

    args = parser.parse_args()

    # If explicit output paths are provided but no output_prefix, derive prefix from output_scaffold path
    if args.output_scaffold and args.output_prefix == 'scaffold':
        # Extract directory and base filename from output_scaffold
        # e.g., "scaffolds/sample--virus_final_scaffold.fasta" -> "scaffolds/sample--virus"
        scaffold_path = args.output_scaffold
        scaffold_dir = os.path.dirname(scaffold_path)
        scaffold_base = os.path.basename(scaffold_path)
        # Remove _final_scaffold.fasta suffix to get the prefix
        if scaffold_base.endswith('_final_scaffold.fasta'):
            prefix_base = scaffold_base[:-len('_final_scaffold.fasta')]
            args.output_prefix = os.path.join(scaffold_dir, prefix_base) if scaffold_dir else prefix_base

    # Validate input files exist
    if not os.path.exists(args.best_reference_fasta):
        print(f"Error: Reference fasta file not found: {args.best_reference_fasta}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.contig_fasta):
        print(f"Error: Contig fasta file not found: {args.contig_fasta}", file=sys.stderr)
        sys.exit(1)

    main(args)
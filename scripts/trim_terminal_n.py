import argparse
import pyfastx
import sys
import os

def trim_terminal_n(sequence):
    """
    Trims N bases from both 5' and 3' ends of the sequence.
    Keeps internal N bases.
    """
    if not sequence:
        return ""

    # Handle case-insensitive matching
    seq_upper = sequence.upper()

    # Find first non-N position from 5' end
    start = 0
    for i, base in enumerate(seq_upper):
        if base != 'N':
            start = i
            break
    else:
        # All N's in sequence
        return ""

    # Find first non-N position from 3' end (searching backwards)
    end = len(seq_upper)
    for i in range(len(seq_upper) - 1, -1, -1):
        if seq_upper[i] != 'N':
            end = i + 1
            break

    # Return trimmed sequence from original (preserves case)
    return sequence[start:end]


def main(args):
    try:
        input = pyfastx.Fasta(args.input_fasta)

        if len(input) == 0:
            print("INFO: Input FASTA file is empty. Generating empty output file.")
            open(args.output_fasta, 'w').close()  # Create empty fasta
            print("---")
            print("INFO: Processing complete.")
            sys.exit(0)
        else:
            with open(args.output_fasta, 'w') as f_out:
                for seq in input:
                    seq_name = seq.name
                    seq_sequence = seq.seq
                    original_length = len(seq_sequence)

                    # Trim terminal N's from both ends
                    trimmed_sequence = trim_terminal_n(seq_sequence)
                    trimmed_length = len(trimmed_sequence)

                    # Only write sequences that have content after trimming
                    if trimmed_length > 0:
                        n_trimmed_5prime = original_length - trimmed_length - (original_length - seq_sequence.upper().rfind(trimmed_sequence[0]))
                        n_trimmed_3prime = original_length - trimmed_length - n_trimmed_5prime

                        # Calculate from actual positions
                        seq_upper = seq_sequence.upper()
                        start_pos = seq_upper.find(trimmed_sequence[0].upper())
                        n_trimmed_5prime = start_pos
                        n_trimmed_3prime = original_length - (start_pos + trimmed_length)

                        print(f"{seq_name}: Original={original_length}bp, Trimmed={trimmed_length}bp, "
                              f"5'-N removed={n_trimmed_5prime}, 3'-N removed={n_trimmed_3prime}")

                        f_out.write(f">{seq_name}\n{trimmed_sequence}\n")
                    else:
                        print(f"{seq_name}: All N's - sequence discarded (original={original_length}bp)")

            print("---")
            print("INFO: Processing complete.")

            # Cleanup: Remove pyfastx index files
            index_file = f"{args.input_fasta}.fxi"
            if os.path.exists(index_file):
                try:
                    os.remove(index_file)
                    print(f"INFO: Removed index file: {index_file}")
                except Exception as e:
                    print(f"WARNING: Could not remove {index_file}: {e}")


    except FileNotFoundError:
        print(f"ERROR: Input FASTA file not found at '{args.input_fasta}'", file=sys.stderr)
        sys.exit(1)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Trim N bases from 5' and 3' ends of sequences. Internal N's are preserved.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input_fasta', required=True, help='Input fasta file to be processed')
    parser.add_argument('-o', '--output_fasta', required=True, help='Path to the output fasta file')
    args = parser.parse_args()
    main(args)
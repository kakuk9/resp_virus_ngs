"""
Functions for nucmer alignment
"""
import pyfastx
import subprocess
from itertools import combinations

def trim_contigs_correct_orient(coords, fasta, output_fasta):
    """
    - trim contig based on nucmer coord
    - correct orientation of contig
    - add suffix to contig name based on trimming status
    """
    fasta = pyfastx.Fasta(fasta, build_index=True)
    temp_fasta_content = []
    for coord in coords:
        if int(coord[2]) < int(coord[3]):
            sequence = str(fasta[coord[-1]])[int(coord[2])-1:int(coord[3])]
            suffix = "_o"
        elif int(coord[2]) > int(coord[3]):
            sequence = str(fasta[coord[-1]].antisense)[int(coord[3])-1:int(coord[2])]
            suffix = "_rc"
        if len(sequence) != int(coord[-5]):
            suffix += "_untrimmed"
        else:
            suffix += "_trimmed"
        temp_fasta_content.append(f">{coord[-1]}{suffix}\n{sequence}")
    with open(output_fasta, "w") as f:
        f.write("\n".join(temp_fasta_content))

def run_nucmer(reference_fasta, query_fasta, out_prefix, coverage_threshold):
    """
    - run nucmer
    - run show-coords
    - filter nucmer coords based on coverage threshold
    - return coords
    """
    try:
        subprocess.run(
            f'nucmer {reference_fasta} {query_fasta} -p {out_prefix}',
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        subprocess.run(
            f'show-coords -rcl {out_prefix}.delta > {out_prefix}.coords',
            shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}\nReturn code: {e.returncode}\nOutput: {e.output}\nError: {e.stderr.decode() if e.stderr else ''}")
        return []

    try:
        with open(f'{out_prefix}.coords') as f:
            coords = [
                x.strip().replace("|", "").split()
                for x in f.readlines()[5:]
                if float(x.strip().replace("|", "").split()[-3]) >= coverage_threshold
            ]
        return coords
    except Exception as e:
        print(f"Error reading or processing coords file: {e}")
        return []

def downsize_contig_pool_range(coords):
    """
    - remove contigs that are completely overlapping with another contig
    - remove contigs that are partially overlapping with other contigs
    - calculate gap to next contig
    """
    to_be_removed = set()
    # completely overlapping with another contig
    for combo in combinations(coords, 2):
        q_range = set(range(int(combo[0][0]), int(combo[0][1])+1))
        s_range = set(range(int(combo[1][0]), int(combo[1][1])+1))
        if q_range.issubset(s_range):
            to_be_removed.add(combo[0][-1])
        elif s_range.issubset(q_range):
            to_be_removed.add(combo[1][-1])
    coords = [x for x in coords if x[-1] not in to_be_removed]

    partially_overlapping = set()
    repeat_partially_overlapping = True
    while repeat_partially_overlapping:
        for i in range(len(coords)):
            coords_len = len(coords)
            if i != 0 and i != len(coords)-1:
                neighbouring_contigs = set(list(range(int(coords[i-1][0]), int(coords[i-1][1])+1)) + list(range(int(coords[i+1][0]), int(coords[i+1][1])+1)))
                contig_to_check = set(list(range(int(coords[i][0]), int(coords[i][1])+1)))
                if contig_to_check.issubset(neighbouring_contigs):
                    partially_overlapping.add(coords[i][-1])
        coords = [x for x in coords if x[-1] not in partially_overlapping]
        if coords_len == len(coords):
            repeat_partially_overlapping = False

    ref_len = int(coords[0][7])
    for i in range(len(coords)):
        if i < len(coords)-1:
            GAP2NEXT = int(coords[i+1][0]) - int(coords[i][1])
        elif i == len(coords)-1:
            GAP2NEXT = ref_len - int(coords[i][1])
        coords[i].append(f"{int(GAP2NEXT)}")
    return coords

# Overlapping regions alignment
def getSubseq_aln(current_query, next_query, i_start, i_end, i_ort, n_start, n_end, n_ort, gap, fasta):
    if i_ort == "+":
        i_seq = str(fasta[current_query])[i_start-1:i_end][gap:]
    elif i_ort == "-":
        i_seq = str(fasta[current_query].antisense)[i_end-1:i_start][gap:]
    if n_ort == "+":
        n_seq = str(fasta[next_query])[n_start-1:n_end][:abs(gap)]
    elif n_ort == "-":
        n_seq = str(fasta[next_query].antisense)[n_end-1:n_start][:abs(gap)]

    from Bio import Align
    aligner = Align.PairwiseAligner(scoring="blastn", mode = 'global')
    alignments = aligner.align(i_seq, n_seq)
    pident = "{:.1f}".format(alignments[0].counts()[1]/sum(alignments[0].counts())*100)
    n_gaps, n_ident, n_mismatches = alignments[0].counts()
    aln1 = alignments[0][0, :]
    aln2 = alignments[0][1, :]    

    if pident == "100.0":
        overlap_consensus = aln1
    else:
        overlap_consensus = ""
        from Bio.Data import IUPACData
        d = {v:k for k,v in IUPACData.ambiguous_dna_values.items()}

        for pos in range(0, len(aln1)):
            if aln1[pos] == aln2[pos]:
                overlap_consensus += aln1[pos]
            else:
                nt="".join(sorted([aln1[pos], aln2[pos]]))
                if "-" not in nt:
                    # Check if this nucleotide combination exists in IUPAC dictionary
                    if nt in d:
                        am_nt = d[nt]
                        overlap_consensus += am_nt
                    else:
                        # If not in dictionary (e.g., same nucleotide sorted), just use the first one
                        overlap_consensus += aln1[pos]
                else:
                    # "-" will be sorted first, here we ignore gap but use lower case to flag the position
                    overlap_consensus += nt[1].lower()
    aln_output = f"GAP2NEXT: {gap}\n{current_query}\t{i_start}\t{i_end}\t{i_ort}\n{next_query}\t{n_start}\t{n_end}\t{n_ort}\nALN_SCORE: {alignments[0].score}\nPident: {pident}%\n{n_gaps} gaps, {n_ident} identities, {n_mismatches} mismatches\n{alignments[0].format()}\n"
    return aln_output, overlap_consensus
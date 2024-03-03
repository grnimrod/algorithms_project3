
import argparse

# Define the args for the command-line call
def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform global linear sequence alignment.")
    parser.add_argument("fasta_file", type=str, help="Path to input FASTA file containing sequences")
    parser.add_argument("-m", "--matrix", type=str, help="Path to scoring matrix file in Phylip-like format", required=True)
    parser.add_argument("-g", "--gap", type=int, help="Gap penalty", required=True)
    parser.add_argument('--hide-alignments', action='store_true', help='Hide aligned sequences')
    
    args = parser.parse_args()
    
    if args.matrix is None or args.gap is None:
        parser.error("-m <matrix_file> or -g <gapscore> is missing")

    return args


# read the fasta file containing the 2 sequences
# format looks like:
# >seq1
# acgtgtcaacgt 
# >seq2
# acgtcgtagcta 
def read_fasta_file(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:]
                sequence = ""
            else:
                sequence += line
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences


# Read the .txt file containing the substitution matrix
# uses the phylip-like format
def read_phylip_like_matrix(matrix_file):
    matrix = {}
    with open(matrix_file, 'r') as file:
        num_chars = int(file.readline().strip())
        characters = ['A', 'C', 'G', 'T']
        for char in characters:
            matrix[char] = {}
        for line in file:
            line = line.strip().split()
            char = line[0]
            scores = list(map(int, line[1:]))
            for i, score in enumerate(scores):
                matrix[char][characters[i]] = score
    return matrix


# Global alignment with linear gapcost
# takes 5 args, A and B are str representing sequences
# g is gap cost, sub_matrix the substitution_matrix
# and show_alignment defines if the alignment should be included

def global_linear(sequence1, sequence2, g, sub_matrix, hide_alignments=True):
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()
    
    # Get cost of match/mismatch
    def cost(x, y):
        return sub_matrix[x][y]
    
    # Recursively fill out the matrix
    def C(i, j):
        if T[i][j] is None:
            if i == 0 and j == 0: # Base case
                T[i][j] = 0

            elif i == 0: # Only seq A traversed
                T[i][j] = C(i, j - 1) + g
                
            elif j == 0: # Only seq B traversed
                T[i][j] = C(i - 1, j) + g

            else: # Look at 3 prior connected cells
                match = C(i-1, j-1) + cost(sequence1[i-1], sequence2[j-1]) # diagonal move
                insert = C(i, j-1) + g 
                deletion = C(i-1, j) + g 
                T[i][j] = min(match, insert, deletion)
        return T[i][j]

    n = len(sequence1)
    m = len(sequence2)
    T = [[None] * (m + 1) for _ in range(n + 1)] # make table T

    opt_cost = C(n, m) # optimal cost

    if not hide_alignments:
        return opt_cost, []
        
    # Now to find the possible alignments with backtracking - does not run in quadratic time, but only way i can think of to show all alignemnts
    def backtracking(n, m, aligned_seq1="", aligned_seq2=""):
        if n == 0 and m == 0: # Base case
            return [[aligned_seq1, aligned_seq2]]
    
        alignments = [] # initialize list to hold pairs of seqs
    
        if n > 0 and m > 0 and T[n][m] == T[n - 1][m - 1] + cost(sequence1[n - 1], sequence2[m - 1]): # DIAGONAL move
            alignments.extend(backtracking(n - 1, m - 1, sequence1[n - 1] + aligned_seq1, sequence2[m - 1] + aligned_seq2)) # recursively call backtracking 
        if n > 0 and T[n][m] == T[n - 1][m] + g: # LEFT move
            alignments.extend(backtracking(n - 1, m, sequence1[n - 1] + aligned_seq1, "-" + aligned_seq2))
        if m > 0 and T[n][m] == T[n][m - 1] + g: # VERTICAL move
            alignments.extend(backtracking(n, m - 1, "-" + aligned_seq1, sequence2[m - 1] + aligned_seq2))
    
        return alignments

    aligned_sequences = backtracking(n, m)
    
    return opt_cost, aligned_sequences


# write the output to a fasta file 
def write_alignment_to_fasta(aligned_sequences, output_file):
    with open(output_file, 'w') as file:
        if isinstance(aligned_sequences, int):  # If aligned_sequences is an integer (optimal cost)
            file.write(">Optimal cost: {}\n".format(aligned_sequences))
        else:  # If aligned_sequences is a list of aligned sequences
            optimal_cost, sequence_pairs = aligned_sequences
            file.write(">Optimal cost: {}\n".format(optimal_cost))
            for i, sequence_pair in enumerate(sequence_pairs):
                file.write(">Alignment {}\n".format(i + 1))
                for j, sequence in enumerate(sequence_pair):
                    #file.write(">seq{}\n".format(j + 1))
                    file.write(sequence + "\n")


def main():
    args = parse_arguments()
    fasta_file = args.fasta_file
    matrix_file = args.matrix
    gap_penalty = args.gap
    hide_alignments = args.hide_alignments

    sequences = read_fasta_file(fasta_file)
    scoring_matrix = read_phylip_like_matrix(matrix_file)

    aligned_sequences = {}
    for seq_id1 in sequences:
        aligned_sequences[seq_id1] = {}
        for seq_id2 in sequences:
            if seq_id1 != seq_id2:
                aligned_sequences[seq_id1][seq_id2] = global_linear(sequences[seq_id1], sequences[seq_id2], gap_penalty, scoring_matrix, hide_alignments=not hide_alignments)

    # Print the aligned sequences in FASTA alignment format
    for seq_id1 in aligned_sequences:
        for seq_id2 in aligned_sequences[seq_id1]:
            alignment = aligned_sequences[seq_id1][seq_id2]
            print(f"Alignment between {seq_id1} and {seq_id2}:")
            if isinstance(alignment, list):
                for aligned_seq in alignment:
                    print(aligned_seq)
            else:
                print("Optimal cost: {}".format(alignment))
    
    # Use the write function to make output file
    write_alignment_to_fasta(aligned_sequences, "linear_aligned_sequences.fasta")

if __name__ == "__main__":
    main()
import sys
from Bio import SeqIO
from itertools import combinations
from global_align_linear import *
from extend_msa import *

# Define function to read FASTA file
def read_fasta_file(file_path):
    """
    Reads a FASTA file and returns a dictionary of sequences.
    """
    sequences = {}
    with open(file_path, "r") as fasta_file:
        for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), start=1):
            sequences[f"seq{i}"] = str(record.seq)
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



# Define functions to generate pairwise combinations and run alignment script
def generate_pairwise_combinations(sequences):
    
    #Generates all possible pairwise combinations of sequences.
    
    pairwise_combinations = list(combinations(sequences, 2))
    return pairwise_combinations


def main():
    # Read sequences from input FASTA file
    fasta_file_path = sys.argv[1]
    sequences = read_fasta_file(fasta_file_path)
    #print(sequences)

    #read the scoring matrix in phylip fomat
    phylip_matrix = sys.argv[2]
    sub_matrix = read_phylip_like_matrix(phylip_matrix)
    #print(sub_matrix)

    # Generate pairwise combinations
    pairwise_combinations = generate_pairwise_combinations(sequences)
    #print(pairwise_combinations)

    # Create a dictionary to store sequences with initial values of 0
    sequence_scores = {seq_id: 0 for seq_id in sequences}
    #getteing the combinde alignmet score
    sp_approx_alignment_score=0
    #getting the list of alignments
    alignment=[]

    #finding the sequenst with the least distens to the other sequences
    for seqid1, seqid2 in pairwise_combinations:
        alignment_score = global_linear(sequences[seqid1], sequences[seqid2], 5, sub_matrix, hide_alignments=False) 
        sequence_scores[seqid2]+= alignment_score[0]
        sequence_scores[seqid1]+= alignment_score[0]
    
    # calculating the distens and getting the alignment to all other sequences  to the nereast one
    for i in sequences:
        if i!=min(sequence_scores):
            sp_alignment_score = global_linear(sequences[min(sequence_scores)], sequences[i], 5, sub_matrix, hide_alignments=True)    
            sp_approx_alignment_score+=sp_alignment_score[0] 
            converted_list = [list(chars) for chars in zip(*sp_alignment_score[1][0])]
            #print(converted_list)
            alignment.append(converted_list)
            #alignment[i]=converted_list

    #print(alignment)
    M=alignment[0]
    for i in range(1,len(alignment)):
        A=alignment[i]
        M=extend_msa(M, A)

    print(M)
    print(sp_approx_alignment_score)
            
    #extend_msa(M, A)


if __name__ == "__main__":
    main()
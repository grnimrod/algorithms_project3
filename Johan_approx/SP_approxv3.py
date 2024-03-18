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

def construct_alignment_strings(alignment_data):
    alignment_strings = {}
    num_sequences = len(alignment_data[0])  # Number of sequences

    for i in range(num_sequences):
        sequence_string = ""
        for row in alignment_data:
            sequence_string += row[i]
        alignment_strings[f"seq{i + 1}"] = sequence_string
    
    return alignment_strings

def write_alignment_to_fasta(alignment_strings, output_file):
    with open(output_file, 'w') as file:
        for sequence_id, alignment_string in alignment_strings.items():
            file.write(">" + sequence_id + "\n")
            file.write(alignment_string + "\n")

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
    calculate_avarage=0
    count=0

    #getting the list of alignments
    alignment=[]

    #finding the sequenst with the least distens to the other sequences
    for seqid1, seqid2 in pairwise_combinations:
        alignment_score = global_linear(sequences[seqid1], sequences[seqid2], 5, sub_matrix, hide_alignments=False) 
        sequence_scores[seqid2]+= alignment_score[0]
        sequence_scores[seqid1]+= alignment_score[0]
        calculate_avarage += alignment_score[0]
        count+=1

    


    # calculating the distens and getting the alignment to all other sequences  to the nereast one
    for i in sequences:
        if i!=min(sequence_scores):
            sp_alignment_score = global_linear(sequences[min(sequence_scores)], sequences[i], 5, sub_matrix, hide_alignments=True)    
            #print(min(sequence_scores),i,sp_alignment_score[0])
            converted_list = [list(chars) for chars in zip(*sp_alignment_score[1][0])]
            alignment.append(converted_list)
    
    # #print(alignment)
    M=alignment[0]
    print(alignment)
    for i in range(1,len(alignment)):
        A=alignment[i]
        M=extend_msa(M, A)
        #print(M)

    


    
    #not nessesary at the moment
    #output_file = "alignment.fasta"
    #write_alignment_to_fasta(construct_alignment_file,output_file)
    
    #format the data
    construct_alignment_file=construct_alignment_strings(M)

    #getting the multile alignment score

    cost = [[0, 5, 2, 5, 5],  # A
        [5, 0, 5, 2, 5],  # C
        [2, 5, 0, 5, 5],  # G
        [5, 2, 5, 0, 5],  # T
        [5, 5, 5, 5, 0]]  #-'

    dict_str2seq = {'a':0, 'c':1, 'g':2, 't':3, 'A':0, 'C':1, 'G':2, 'T':3, '-':4, 'N':0, 'R':0, 'S':0}

    def str2seq(s):
        try:
            seq = [dict_str2seq[c] for c in list(s)]
            return seq
        except KeyError as e:
            print("ERROR: Illegal character", e, "in input string.")
            sys.exit(1)
    
    
    def compute_sp_score(sequences):
    # Convert sequences dictionary to list of sequences
        seq_list = list(sequences.values())
    # Convert sequences to sequences of indices using str2seq function
        row = [str2seq(s) for s in seq_list]
    # Compute the score of each induced pairwise alignment
        score = 0
        for i in range(len(row)):
            for j in range(i + 1, len(row)):
                if len(row[i]) != len(row[j]):
                    print("ERROR: Rows", i, "and", j, "have different lengths.")
                    sys.exit(1)
                for c in range(len(row[i])):
                    score += cost[row[i][c]][row[j][c]]
        return score
    
    print(compute_sp_score(construct_alignment_file))




if __name__ == "__main__":
    main()
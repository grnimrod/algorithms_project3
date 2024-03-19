
# Import helper functions
from msa_sp_score_3k import *
from score_matrix import initiate_score_matrix
from read_fasta import read_fasta_file

import argparse
#################
# Function call example:
# (ctib) C:\Users\riber\Programming\masters\algorithms\project3>python src\sp_exact_3.py data\testdata_short.txt -g 5 -m data\score_matrix.txt
# 1st arg: sequence fasta file (path)
# 2nd arg: -g, <gap_score>
# 3rd arg: -m, <score_matrix> (path)
# 4th arg: --hide-alignments, hides the alignments hence skipping backtracking (NOT MANDOTORY)

# Output:
 
# Alignment:
# GTTCCGAAAGGCTAGCGCTAGGC-GCC-
# A-T--G-GAT-TT-AT-CTGCTC-TTCG
# --T--G-CATGCTGAAACTTCTCAACCA
#
# Optimal sum-of-pairs score: 198



#################
# Setup for arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform multiple sequence alignment of exactly 3 sequences")
    parser.add_argument("fasta_file", type=str, help="Path to input FASTA file containing sequences")
    parser.add_argument("-g", "--gap", type=int, help="Gap penalty", required=True)
    parser.add_argument("-m", "--matrix", type=str, help="Path to scoring matrix file", required=True)
    parser.add_argument('--hide-alignments', action='store_true', help='Hide aligned sequences')
    
    args = parser.parse_args()
    
    if args.gap is None or args.matrix is None:
        parser.error("-g <gapscore> or -m <matrix> is missing")

    return args




#####################################
############ sp_exact_3 #############
#####################################
def sp_exact_3(A, B, C, gap_cost, score_matrix, hide_alignments=False):

    # make 3d matrix filled with inf (its a minimyzation)
    n, m, l = len(A), len(B), len(C)
    dp = [[[float("inf")] * (l + 1) for _ in range(m + 1)] for _ in range(n + 1)]

    # first cell
    dp[0][0][0] = 0

    # Make edge-layers with gapcost
    for i in range(1, n+1):
        dp[i][0][0] = dp[i-1][0][0] + gap_cost
    for j in range(1, m+1):
        dp[0][j][0] = dp[0][j-1][0] + gap_cost
    for k in range(1, l+1):
        dp[0][0][k] = dp[0][0][k-1] + gap_cost 
    

    # Fill out 3D table
        # if states refereing to the score (see slide 18, pp4.1)
        # remember we are looking at 7 cells now since its 3d
    for i in range(0, n+1):
        for j in range(0, m+1):
            for k in range(0, l+1):
                v0, v1, v2, v3, v4, v5, v6, v7 = [float("inf")] * 8
                
                if i==0 and j==0 and k==0: 
                    v0 = 0

                if i>0 and j>0 and k>0:                 # matching all 3 seq
                    v1 = dp[i-1][j-1][k-1] + score_matrix[A[i-1]][B[j-1]] + score_matrix[A[i-1]][C[k-1]] + score_matrix[B[j-1]][C[k-1]]
                if i>0 and j>0 and k>=0:                # matching only 2 seq
                    v2 = dp[i-1][j-1][k] + score_matrix[A[i-1]][B[j-1]] + gap_cost + gap_cost
                if i>0 and j>=0 and k>0:                # matching only 2 seq
                    v3 = dp[i-1][j][k-1] + score_matrix[A[i-1]][C[k-1]] + gap_cost + gap_cost
                if i>=0 and j>0 and k>0:                # matching only 2 seq
                    v4 = dp[i][j-1][k-1] + score_matrix[B[j-1]][C[k-1]] + gap_cost + gap_cost
                if i>0 and j>=0 and k>=0:               # No matches
                    v5 = dp[i-1][j][k] + gap_cost + gap_cost
                if i>=0 and j>0 and k>=0:               # No matches
                    v6 = dp[i][j-1][k] + gap_cost + gap_cost
                if i>=0 and j>=0 and k>0:               # No matches
                    v7 = dp[i][j][k-1] + gap_cost + gap_cost
                
                dp[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)

    # return last cell or continue with backtracking
    score = dp[n][m][l]

    if hide_alignments == True:
        return score

    # Backtracking 
    def backtracking(A, B, C, dp, score_matrix, gap_cost):
        # Function for backtracking from the score return above.
        # checks if the current cells score corosponds to the 
        # scores of the any of the 7 prior cells it could 
        # potentially come from. 
        # As long as i, j, k arent 0 thenthey havent reach their corosponding edge
        ###################

        alignment = []
        i, j, k = n, m, l

        while i > 0 or j > 0 or k > 0:
            # All 3 match and is the same score going in diagonal
            if i>0 and j>0 and k>0 and dp[i][j][k] == dp[i-1][j-1][k-1] + score_matrix[A[i-1]][B[j-1]] + score_matrix[A[i-1]][C[k-1]] + score_matrix[B[j-1]][C[k-1]]:
                alignment.append([A[i-1], B[j-1], C[k-1]])
                i -= 1
                j -= 1
                k -= 1

            # The 3 combinations of "2 sequences aligned and 1 gap" 
            elif i>0 and j>0 and dp[i][j][k] == dp[i-1][j-1][k] + score_matrix[A[i-1]][B[j-1]] + gap_cost + gap_cost:
                alignment.append([A[i-1], B[j-1], "-"])
                i -= 1
                j -= 1
            elif i>0 and k>0 and dp[i][j][k] == dp[i-1][j][k-1] + score_matrix[A[i-1]][C[k-1]] + gap_cost + gap_cost:
                alignment.append([A[i-1], "-", C[k-1]])
                i -= 1
                k -= 1
            elif j>0 and k>0 and dp[i][j][k] == dp[i][j-1][k-1] + score_matrix[B[j-1]][C[k-1]] + gap_cost + gap_cost:
                alignment.append(["-", B[j-1], C[k-1]])
                j -= 1
                k -= 1

            # The 3 combinations of "2 gaps, so no aligned seqeuences"
            elif i>0 and dp[i][j][k] == dp[i-1][j][k] + gap_cost + gap_cost:
                alignment.append([A[i-1], "-", "-"])
                i -= 1
            elif j>0 and dp[i][j][k] == dp[i][j-1][k] + gap_cost + gap_cost:
                alignment.append(["-", B[j-1], "-"])
                j -= 1
            elif k>0 and dp[i][j][k] == dp[i][j][k-1] + gap_cost + gap_cost:
                alignment.append(["-", "-", C[k-1]])
                k -= 1

        alignment.reverse() # since backtracking goes backwards, everything should be reversed
        return alignment

    alignment_cols = backtracking(A, B, C, dp, score_matrix, gap_cost)
    alignment_rows = [[],[],[]] # initialize row list
    
    for a, b, c in alignment_cols: # get the letters from each col and insert as rows
        alignment_rows[0].append(a)
        alignment_rows[1].append(b)
        alignment_rows[2].append(c)

    joined_strings = [''.join(alignment) for alignment in alignment_rows] # join the strings
    
    return [joined_strings, score]


# Main collection of functions and argument parses
def main():
    args = parse_arguments()
    fasta_file = args.fasta_file
    gap_cost = args.gap
    matrix_file = args.matrix
    hide_alignments = args.hide_alignments

    # load fasta file with seqs
    sequence_dct = read_fasta_file(fasta_file)
    
    sequence_lst = [] # convert dct values to list
    iterator = iter(sequence_dct.values()) # Extracting the sequences without their corresponding headers and putting them into a list
    for _ in range(len(sequence_dct)):
        sequence_lst.append(next(iterator))

    # load score matrix
    score_matrix = initiate_score_matrix(matrix_file)
    
    seq1 = sequence_lst[0]
    seq2 = sequence_lst[1]
    seq3 = sequence_lst[2]
    
    aligned_sequences = sp_exact_3(seq1, seq2, seq3, gap_cost, score_matrix, hide_alignments = hide_alignments)

    # print a pretty output
    for i, obj in enumerate(sequence_dct):
        print(f">{obj}\n{aligned_sequences[0][i]}")
    print(f"\nOptimal sum-of-pairs score: {aligned_sequences[1]}")
        

if __name__ == "__main__":
    main()



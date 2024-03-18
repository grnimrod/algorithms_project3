from read_fasta import read_fasta_file
from score_matrix import initiate_score_matrix

sequences_dict = read_fasta_file("./data/testdata_short.txt")
sequences_list = []

# Extracting the sequences without their corresponding headers and putting them into a list
iterator = iter(sequences_dict.values())

for _ in range(len(sequences_dict)):
    sequences_list.append(next(iterator))

"""
As the sequences are put into a list reserving their original
order, it's more straightforward to use them as input parameters
for the alignment function in different testing scenarios
"""


def sp_exact_3(A, B, C, path_to_sc_matrix, gap_cost):
    sc_matrix = initiate_score_matrix(path_to_sc_matrix)

    """
    Create three dimensional dynamic programming table
    We fill it out with 'inf', as we are minimizing the score
    """
    n, m, l = len(A), len(B), len(C)
    dp = [[[float('inf')] * (l + 1) for _ in range(m + 1)] for _ in range(n + 1)]

    # Set dp[0][0][0] to 0
    dp[0][0][0] = 0

    # Initialize the first layer, row and column of the dp table with gap costs
    for i in range(1, n + 1):
        dp[i][0][0] = dp[i - 1][0][0] + gap_cost
    for j in range(1, m + 1):
        dp[0][j][0] = dp[0][j - 1][0] + gap_cost
    for k in range(1, l + 1):
        dp[0][0][k] = dp[0][0][k - 1] + gap_cost

    # Filling out the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            for k in range(1, l + 1):
                v1, v2, v3, v4, v5, v6, v7 = [0] * 7
                if i > 0 and j > 0 and k > 0:
                    v1 = dp[i - 1][j - 1][k - 1] + sc_matrix[A[i - 1]][B[j - 1]] + sc_matrix[A[i - 1]][C[k - 1]] + sc_matrix[B[j - 1]][C[k - 1]]
                if i > 0 and j > 0 and k >= 0:
                    v2 = dp[i - 1][j - 1][k] + sc_matrix[A[i - 1]][B[j - 1]] + gap_cost + gap_cost
                if i > 0 and j >= 0 and k > 0:
                    v3 = dp[i - 1][j][k - 1] + gap_cost + sc_matrix[A[i - 1]][C[k - 1]] + gap_cost
                if i >= 0 and j > 0 and k > 0:
                    v4 = dp[i][j - 1][k - 1] + gap_cost + gap_cost + sc_matrix[B[j - 1]][C[k - 1]]
                if i > 0 and j >= 0 and k >= 0:
                    v5 = dp[i - 1][j][k] + gap_cost + gap_cost
                if i >= 0 and j > 0 and k >= 0:
                    v6 = dp[i][j - 1][k] + gap_cost + gap_cost
                if i >= 0 and j >= 0 and k > 0:
                    v7 = dp[i][j][k - 1] + gap_cost + gap_cost
                dp[i][j][k] = min(v1, v2, v3, v4, v5, v6, v7)

    return dp[n][m][l]
    
    # Backtracking
    

print(sp_exact_3(sequences_list[0], sequences_list[1], sequences_list[2], "./data/score_matrix.txt", 5))
from score_matrix import initiate_score_matrix
from sp_exact_3 import sp_exact_3

seq1 = "AA"
seq2 = "AA"
seq3 = "AA"

score_matrix = initiate_score_matrix("./data/score_matrix.txt")

print(sp_exact_3(seq1, seq2, seq3, score_matrix, 5))
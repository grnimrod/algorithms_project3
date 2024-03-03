
# Let M be an multiple alignment (represented as a python a list of lists/columns).
#
# Let A be a pairwise alignment (represented as a python a list of lists/columns).
#
# The first row of M and A is the same string if gaps are removed
#
# Example for slide 11 in SP-MSA-Approx.pdf:

M = [['a','a','c'], ['-','t','t'], ['-','t','-'], ['c','c','c'], ['g','-','g'], ['t','t','a']]

A = [['a','a'], ['c','c'], ['g','g'], ['-','g'], ['t','t'], ['-','c']]

def print_alignment(A, highlight = -1):
    cols = len(A)
    if cols > 0:
        rows = len(A[-1])
        for i in range(rows):
            for j in range(cols):
                if j == highlight:
                    print(A[j][i].upper(), end=" ")
                else:
                    print(A[j][i], end=" ")
            print()
    else:
        print("Empty")
        
# Let MA be an MSA that is the extension of M with A cf. Gusfields 2-approximation algorithm.

MA = []
i = 0
j = 0

print("------------------------------------")
print("Start")
print("------------------------------------")  
print()
print("M:")
print_alignment(M)
print()
print("A:")
print_alignment(A)
print()

round = 0
while i < len(M) and j < len(A):

    print("------------------------------------")
    print("Round:", round)
    print("------------------------------------")    
    print()
    print("MA:")
    print_alignment(MA)
    print()
    print("i:", i)
    print("j:", j)
    print()
    print("M:")
    print_alignment(M,i)
    print()
    print("A:")
    print_alignment(A,j)
    print()
    round = round + 1

    # Invariant: (1) MA is a valid merge of all columns before column i in M
    # and all columns before column in A, and (2) the first row of M and A up
    # to (but not including) column i and j respectively is the same string
    # if gaps are removed.
    
    if M[i][0] == '-' and A[j][0] == '-':
        print("Apply case 1")
        print()
        # Case 1: The next column in MA is column i in M extended with the second symbol
        # in column j in A.
        M[i].append(A[j][1])
        MA.append(M[i])
        i = i + 1
        j = j + 1

    elif M[i][0] == '-' and A[j][0] != '-':
        print("Apply case 2")
        print()
        # Case 2: A[j][0] is a character, so the second symbol in column j in A, A[j][1],
        # must be in the column of MA that is the column in M where the first symbol corresponds
        # to A[j][0]. By the invariant, this column in M is the next column in M, where the first
        # symbol is a character, so we just moved forward in M until we find this column.
        M[i].append('-')
        MA.append(M[i])
        i = i + 1

    elif M[i][0] != '-' and A[j][0] == '-':
        print("Apply case 3")
        print()
        # Case 3: M[i][0] is a character, so column i in M must be in the column of MA that also
        # contains the second symbol from the column in A, where the first symbol is the character
        # corresponding to M[i][0]. By the invariant, this column in A is the next column in A,
        # where the first symbol is a character, so we just add columns from A to MA until we
        # find this column.
        c = ['-']*len(M[i])
        c.append(A[j][1])
        MA.append(c)
        j = j + 1

    elif M[i][0] != '-' and A[j][0] != '-':
        print("Apply case 4")
        print()
        # Case 4: By the invariant the characters M[i][0] and A[j][0] are at the same position
        # in the string spelled by the row of M and A if gaps are removed. The next column in
        # MA is thus column i in M extended with the second symbol in column j in A.
        M[i].append(A[j][1])
        MA.append(M[i])
        i = i + 1
        j = j + 1

if i < len(M):
    # add the remaining coloumns of M to MA
    while i < len(M):
        MA.append(M[i].append('-'))
        i = i + 1
        
if j < len(A):
    # add the remaining columns of A to MA
    k = len(MA[-1])
    while j < len(A):
        c = ['-']*(k-1)
        c.append(A[j][1])
        MA.append(c)
        j = j + 1

print("------------------------------------")
print("Done")
print("------------------------------------")
print()
print("MA:")
print_alignment(MA)
                      
    

            
        

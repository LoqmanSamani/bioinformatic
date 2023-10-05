import numpy as np
from datetime import datetime as t

def hamming_distance(seq1, seq2):
    """
    Calculate the Hamming distance between two sequences.

    Hamming distance is the number of substitutions (mismatches)
    between two sequences of the same length.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The Hamming distance between the two sequences.
    """
    num_sub = 0

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            num_sub += 1

    return num_sub

def edit_distance_recursive(seq1, seq2):
    """
    Calculate the edit distance between two sequences using a recursive approach.

    Edit distance (Levenshtein distance) is the minimum number of
    insertions, deletions, and substitutions required to transform
    one sequence into another.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The edit distance between the two sequences.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    if len(seq1) == 0:
        return len(seq2)
    
    elif len(seq2) == 0:
        return len(seq1)
    
    else:
        dis_hor = edit_distance_recursive(seq1[:-1], seq2) + 1
        dis_ver = edit_distance_recursive(seq1, seq2[:-1]) + 1
        if seq1[-1] == seq2[-1]:
            dis_diag = edit_distance_recursive(seq1[:-1], seq2[:-1])
        else:
            dis_diag = edit_distance_recursive(seq1[:-1], seq2[:-1]) + 1

    return min(dis_hor, dis_ver, dis_diag)

def edit_distance(seq1, seq2):
    """
    Calculate the edit distance between two sequences using dynamic programming.

    Edit distance (Levenshtein distance) is the minimum number of
    insertions, deletions, and substitutions required to transform
    one sequence into another.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The edit distance between the two sequences.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    d_matrix = np.zeros((len(seq1)+1, len(seq2)+1))

    # Initialize the first row and the first column of the matrix with increasing numbers!
    for i in range(len(seq1)+1):
        d_matrix[i, 0] = i

    for j in range(len(seq2)+1):
        d_matrix[0, j] = j

    for row in range(1, d_matrix.shape[0]):
        for col in range(1, d_matrix.shape[1]):
            dis_hor = d_matrix[row, col-1] + 1
            dis_ver = d_matrix[row-1, col] + 1

            if seq1[row-1] == seq2[col-1]:
                dis_diag = d_matrix[row-1, col-1]
            else:
                dis_diag = d_matrix[row-1, col-1] + 1

            d_matrix[row, col] = min(dis_hor, dis_ver, dis_diag)

    return d_matrix[-1, -1]





t_s = t.now()
# Calculate edit distance using dynamic programming (faster algorithm)
print(edit_distance('Shakcspeer', 'shac  ptea'))
t_e = t.now()

tr_s = t.now()
# Calculate edit distance using the slow recursive algorithm
print(edit_distance_recursive('Shakcspeer', 'shac  ptea'))
tr_e = t.now()

# Calculate the time taken by each algorithm
time1 = tr_e - tr_s  # Time taken by the recursive algorithm
time2 = t_e - t_s    # Time taken by the dynamic programming algorithm

print("Time taken by recursive algorithm:", time1)
print("Time taken by dynamic programming algorithm:", time2)

# Calculate the speedup factor (how many times faster the dynamic programming algorithm is)
speedup_factor = time1 / time2
print("Speedup factor:", speedup_factor)


"""
The outputs:

5.0
5

Time taken by recursive algorithm: 0:00:04.111366

Time taken by dynamic programming algorithm: 0:00:00.000114

Speedup factor: 36064.61403508772

you can see the dynamic programming algorithm is
significantly faster than the recursive algorithm!
"""


import numpy as np



def global_alignment_matrix(pattern, text, alphabet, penalty_scores, skip_score):
    """
    Calculate the global alignment matrix between a pattern and text using dynamic programming.

    Args:
        pattern (str): The pattern sequence.
        text (str): The text sequence.
        alphabet (list): The list of characters in the alphabet.
        penalty_scores (numpy.ndarray): The penalty scores for substitutions between characters in the alphabet.
        skip_score (int): The penalty score for skipping a character in the pattern or text.

    Returns:
        numpy.ndarray: The global alignment matrix.
    """



    # Convert pattern and text to uppercase for case-insensitive alignment
    pattern = pattern.upper()
    text = text.upper()

    # Create a matrix to store the alignment scores
    d_matrix = np.zeros((len(pattern) + 1, len(text) + 1))

    # Initialize the first column of the matrix with skip scores
    sk_s = skip_score
    for i in range(1, len(pattern) + 1):
        d_matrix[i, 0] = sk_s
        sk_s += skip_score

    # Initialize the first row of the matrix with skip scores
    sk_sc = skip_score
    for j in range(1, len(text) + 1):
        d_matrix[0, j] = sk_sc
        sk_sc += skip_score

    # Fill in the matrix using dynamic programming
    for row in range(1, d_matrix.shape[0]):
        for col in range(1, d_matrix.shape[1]):

            # Calculate the scores for horizontal, vertical, and diagonal moves
            dis_hor = d_matrix[row, col - 1] + penalty_scores[-1, alphabet.index(text[col - 1])]
            dis_ver = d_matrix[row - 1, col] + penalty_scores[alphabet.index(pattern[row - 1]), -1]

            if pattern[row - 1] == text[col - 1]:
                dis_diag = d_matrix[row - 1, col - 1]
            else:
                dis_diag = d_matrix[row - 1, col - 1] + penalty_scores[alphabet.index(pattern[row - 1]), alphabet.index(text[col - 1])]

            # Choose the minimum score among the three
            d_matrix[row, col] = min(dis_hor, dis_ver, dis_diag)

    return d_matrix

# Define the alphabet of bases
alphabet = ['A', 'G', 'C', 'T']

# Define the penalty scores for insertion, deletion, and substitution

"""
Since  the  ratio of  transition/transversion  substitution  is ~ 2.1,
we should  have different  scores  for  each of  these  two  different
substitutions.  The first  row of the  penalty  matrix  corresponds to
the first character in the alphabet (A) and the second row corresponds
to  the  second  character (G) and so on until the last row and column
corresponding to skip!
"""

penalty_scores = np.array([
    [0, 2, 4, 4, 8],
    [2, 0, 4, 4, 8],
    [4, 4, 0, 2, 8],
    [4, 4, 2, 0, 8],
    [8, 8, 8, 8, 8]
])


# Define the skip score
skip_score = 8



# Calculate the global alignment matrix
galign_matrix = global_alignment_matrix('GCTAGATGCTA', 'GCGTAGACGCTA', alphabet, penalty_scores, skip_score)


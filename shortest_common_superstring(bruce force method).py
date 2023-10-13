import itertools



def len_overlap(read1, read2, threshold):
    """
    Calculate the length of the overlap between two strings.

    Args:
        read1 (str): The first string.
        read2 (str): The second string.
        threshold (int): The minimum overlap length to consider.

    Returns:
        int: The length of the overlap between read1 and read2,
        or 0 if there's no overlap of at least threshold length.
    """
    start = 0

    while True:
        start = read1.find(read2[:threshold], start)

        if start == -1:
            return 0

        if read2.startswith(read1[start:]):
            return len(read1) - start

        start += 1

def shortest_common_superstring(strings):
    """
    Find the shortest common superstring
    that contains all input strings.

    Args:
        strings (list of str): A list of input strings.

    Returns:
        str: The shortest common superstring
        that combines all input strings.
    """
    scs = None

    for string in itertools.permutations(strings):
        sup = string[0]

        for i in range(len(strings) - 1):
            overlap = len_overlap(string[i], string[i + 1], threshold=3)
            sup += string[i + 1][overlap:]

        if scs is None or len(sup) < len(scs):
            scs = sup

    return scs





# Example usage:

#scs = shortest_common_superstring(['AGCTGGCTAGT', 'AGCTGATGCTGATG', 'ACGTCGTAGC'])


scs = set()

for _ in range(100):
    sc = shortest_common_superstring(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])

    scs.add(sc)




print(len(scs))










    

          

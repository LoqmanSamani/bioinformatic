import itertools
import random




class GreedySCS:





    def len_overlap(self, read1, read2, threshold):
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





    def shortest_common_superstring(self, strings):
        """
        Find the shortest common superstring that contains all input strings.

        Args:
            strings (list of str): A list of input strings.

        Returns:
            str: The shortest common superstring that combines all input strings.
        """
        scs = None

        for string in itertools.permutations(strings):
            sup = string[0]

            for i in range(len(strings) - 1):
                overlap = self.len_overlap(string[i], string[i + 1], threshold=1)
                sup += string[i + 1][overlap:]

            if scs is None or len(sup) < len(scs):
                scs = sup

        return scs
    






    def max_overlap(self, reads, threshold):
        """
        Find the pair of reads with the maximum overlap length.

        Args:
            reads (list of str): A list of input strings.
            threshold (int): The minimum overlap length to consider.

        Returns:
            tuple of (str, str, int): The two reads with the maximum overlap and the overlap length.
        """
        read1 = None
        read2 = None
        best_overlap = 0

        for r1, r2 in itertools.permutations(reads, 2):
            overlap = self.len_overlap(r1, r2, threshold=threshold)
            if overlap > best_overlap:
                read1 = r1
                read2 = r2
                best_overlap = overlap

        return read1, read2, best_overlap
    







    def greedy_scs(self, reads, threshold):
        """
        Build the shortest common superstring using a greedy algorithm.

        Args:
            reads (list of str): A list of input strings.
            threshold (int): The minimum overlap length to consider.

        Returns:
            str: The shortest common superstring constructed using the greedy algorithm.
        """
        read1, read2, best_overlap = self.max_overlap(reads, threshold)

        while best_overlap > 0:
            reads.remove(read1)
            reads.remove(read2)
            random.shuffle(reads)
            reads.append(read1 + read2[best_overlap:])
            read1, read2, best_overlap = self.max_overlap(reads, threshold)

        reads = ''.join(reads)

        return reads
    






# An Example

model = GreedySCS()

scs = model.shortest_common_superstring(['AAC', 'GGD', 'AAR','ATC', 'ACG', 'GHFTD', 'GFDRA'])

greedy_scs = model.greedy_scs(['AAC', 'GGD', 'AAR','ATC', 'ACG', 'GHFTD', 'GFDRA'], 2)

print(scs)

print(greedy_scs)






from itertools import permutations as per



class NaiveOverlap:

    def __init__(self, reads, threshold):

        """
        Initialize the NaiveOverlap class with a list of reads
        and a minimum required overlap length.

        Args:
            reads (list of str): A list of sequences (reads).
            threshold (int): The minimum length of the required overlap.
        """

        self.reads = self.upper(reads)  # Convert all reads to uppercase.
        self.threshold = threshold  # Store the threshold for overlap length.



    def upper(self, reads):
        """
        Convert all reads to uppercase.

        Args:
            reads (list of str): A list of sequences (reads).

        Returns:
            list of str: A list of sequences with all 
            characters converted to uppercase.
        """
        return [read.upper() for read in reads]
    



    def len_overlap(self, read1, read2, threshold):

        """
        Find the length of the overlap between two sequences.

        Args:
            read1 (str): The first sequence.
            read2 (str): The second sequence.
            threshold (int): The minimum length of the required overlap.

        Returns:
            int: The length of the overlap between read1 and read2,
            or 0 if no overlap meets the threshold.
        """

        start = 0  # Initialize the starting position for searching the overlap.


        while True:
            # Search for the first occurrence of a substring of read1 in the beginning of read2,
            # where the length of the substring is determined by the threshold.

            start = read1.find(read2[:threshold], start)


            if start == -1:
                # If no match is found, return 0, indicating no overlap meeting the threshold.

                return 0
            

            if read2.startswith(read1[start:]):
                # If read2 starts with the substring of read1 found at the current start position,
                # it means an overlap is found.
                # Calculate and return the length of the overlap by subtracting the start position
                # from the length of read1.

                return len(read1) - start
            

            start += 1  # Move to the next position in read1 to continue searching for overlaps.




    def naive_overlap(self):

        """
        Find and return all pairs of reads that have an
        overlap meeting the specified threshold.

        Returns:
            dict: A dictionary where keys are pairs of reads that overlap,
            and values are the lengths of the overlaps.
        """

        overlaps = {}


        for read1, read2 in per(self.reads, 2):

            overlap_len = self.len_overlap(read1, read2, threshold=self.threshold)


            if overlap_len > 0:

                overlaps[(read1, read2)] = overlap_len


        return overlaps
    


# Example usage:

reads = ["TTACGT", "GGTACCGT", "ACGTAGACGCTA"]

threshold = 3

naive_overlap_finder = NaiveOverlap(reads, threshold)

overlap_results = naive_overlap_finder.naive_overlap()

print(overlap_results)  

# Output: {('TTACGT', 'ACGTAGACGCTA'): 4}





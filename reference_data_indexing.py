

import random





class Indexing:

    
    def __init__(self, data, len_k_mer):
        """
        Initialize the Indexing class with the reference data and the length of k-mer.

        Args:
        - data (str): The reference sequence to be indexed.
        - len_k_mer (int): The length of k-mer for indexing and searching.

        This constructor creates a hash table for quick pattern matching.

        """
        self.data = data
        self.len_k_mer = len_k_mer
        
        # Create a hash table for efficient pattern matching
        self.hash_table = self.make_hash_table(self.data, self.len_k_mer)







        
        
    def make_hash_table(self, data, k):
        """
        Create a hash table to index substrings of the reference sequence.

        Args:
        - data (str): The reference sequence.
        - k (int): The length of k-mer for indexing.

        Returns:
        - list: A sorted list of tuples containing k-mers and their positions in the reference sequence.

        This method iterates through the reference sequence and stores each k-mer and its position(s) in the hash table.

        """
        hash_table = {}
        
        for i in range(len(data) - k + 1):
            
            if data[i: i+k] not in hash_table:
                
                hash_table[data[i: i+k]] = [i]
                
            else:
                hash_table[data[i: i+k]].append(i)
                
        # Convert the hash table into a sorted list for easy access
        hash_list = [(key, val) for key, val in hash_table.items()]
        sorted_hash_list = sorted(hash_list)
        
        return sorted_hash_list








    
    def query(self, pattern):
        """
        Search for exact matches of a pattern in the reference sequence.

        Args:
        - pattern (str): The pattern to search for in the reference sequence.

        Returns:
        - list: A list of tuples containing matched k-mers and their positions in the reference sequence.

        This method searches for exact matches of the input pattern in the reference sequence using the hash table.

        """
        hits = []

        if len(pattern) >= self.len_k_mer:
            
            for i in range(len(pattern) - self.len_k_mer + 1):

                for j in range(len(self.hash_table)):

                    if pattern[i: i + self.len_k_mer] in self.hash_table[j][0]:

                        hits.append((pattern[i: i + self.len_k_mer], self.hash_table[j][1]))
        else:
            raise ValueError("There is a problem, it may be the length of your pattern. "
                             "(the length of the pattern must be equal or shorter than the length of k_mer)")

        return hits







# Generate a random reference sequence
data = ''.join([random.choice('ATCG') for _ in range(1000000)])


# Create an Indexing object with a k-mer length of 10
model = Indexing(data, 10)


# Generate a random pattern to search for
pattern = ''.join([random.choice('ATCG') for _ in range(15)])


# Query the model to find matches for the pattern
matches = model.query(pattern)



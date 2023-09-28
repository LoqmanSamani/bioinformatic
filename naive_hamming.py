


import random

class NaiveHamming:
    
    def __init__(self):
        pass
    
    # Reverse the pattern
    def reverse_pattern(self, pattern):
        """
        Reverse a given DNA pattern.
        
        Parameters:
        - pattern (str): The DNA pattern to be reversed.
        
        Returns:
        - str: The reversed DNA pattern.
        """
        return pattern[::-1]
    
    def naive_hamming(self, data, pattern, num_mismatches):
        """
        Find approximate matches of a pattern and its reverse complement 
        in a DNA sequence (data).

        Parameters:
        - data (str): A DNA sequence.
        - pattern (str): The pattern to search for.
        - num_mismatches (int): The maximum number of allowed mismatches.

        Returns:
        - list: List of starting indices of approximate matches.
        """
        hits = []  # List to store starting indices of matches
        
        r_pattern = self.reverse_pattern(pattern) # Reverse the pattern using the reverse_pattern method
        
        for i in range(len(data) - len(pattern) + 1):
        
            hd = num_mismatches  # Number of allowed mismatches

            is_match = True

            for j in range(len(pattern)):
            
                if data[i + j] != pattern[j]:
                
                    hd -= 1  # Reduce the number of unmatched characters if there is a mismatch
                
                    if hd == 0:
                    
                        is_match = False
                    
                        break

            if is_match:
            
                hits.append(i)

        if pattern != r_pattern:  # Check if the pattern and its reverse are not the same
        
            for i in range(len(data) - len(r_pattern) + 1):
            
                hd = num_mismatches
            
                is_match = True

                for j in range(len(r_pattern)):
                
                    if data[i + j] != r_pattern[j]:
                    
                        hd -= 1
                    
                        if hd == 0:
                        
                            is_match = False
                        
                            break

                if is_match:
                
                    hits.append(i)
                
        return hits

# Generate a random DNA sequence
data =''.join([random.choice('AATGGGGCGC') for _ in range(10000)])

# Define the pattern to search for
pattern = "ATCGCT"

# Create an instance of the NaiveHamming class
model = NaiveHamming()

# Define the maximum number of allowed mismatches
mismatch = 3

# Find approximate matches of the pattern in the DNA sequence
hits = model.naive_hamming(data, pattern, mismatch)



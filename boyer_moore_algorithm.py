

import random





class BoyerMoorePreProcessing:
    """
    BoyerMoorePreProcessing: This class handles the preprocessing
    steps required for the Boyer-Moore algorithm, including
    creating tables for the bad character and good suffix rules.

    This class includes methods for computing the bad character
    rule, good suffix rule, and other required tables.
    """


    def __init__(self, pattern, alphabet='ATCG'):
        self.pattern = pattern
        self.alphabet = alphabet
        self.alpha_map = dict([(alphabet[i], i) for i in range(len(alphabet))])
        self.bad_char = self.dense_bad_char_tab(self.pattern, self.alpha_map)
        self.lp, self.big_l, self.small_l_prime = self.good_suffix_table(self.pattern)





    def bad_character_rule(self, i, c):

        assert c in self.alpha_map
        assert i < len(self.bad_char)
        ci = self.alpha_map[c]
        return i - (self.bad_char[i][ci] - 1)





    def good_suffix_rule(self, i):

        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]





    def match_skip(self):

        return len(self.small_l_prime) - self.small_l_prime[1]






    def dense_bad_char_tab(self, pattern, alpha_map):

        tab = []
        nxt = [0] * len(alpha_map)
        for i in range(0, len(pattern)):
            c = pattern[i]
            assert c in alpha_map
            tab.append(nxt[:])
            nxt[alpha_map[c]] = i + 1
        return tab




    def good_suffix_table(self, pattern):

        n = self.n_array(pattern)
        lp = self.big_l_prime_array(pattern, n)
        return lp, self.big_l_array(pattern, lp), self.small_l_prime_array(n)






    def n_array(self, s):

        return self.z_array(s[::-1])[::-1]





    def z_array(self, s):

        assert len(s) > 1
        z = [len(s)] + [0] * (len(s) - 1)

        for i in range(1, len(s)):
            if s[i] == s[i - 1]:
                z[1] += 1
            else:
                break

        r, l = 0, 0
        if z[1] > 0:
            r, l = z[1], 1

        for k in range(2, len(s)):
            assert z[k] == 0
            if k > r:
                # Case 1
                for i in range(k, len(s)):
                    if s[i] == s[i - k]:
                        z[k] += 1
                    else:
                        break
                r, l = k + z[k] - 1, k
            else:
                # Case 2
                # Calculate length of beta
                nbeta = r - k + 1
                zkp = z[k - l]
                if nbeta > zkp:
                    # Case 2a: zkp wins
                    z[k] = zkp
                else:
                    # Case 2b: Compare characters just past r
                    nmatch = 0
                    for i in range(r + 1, len(s)):
                        if s[i] == s[i - k]:
                            nmatch += 1
                        else:
                            break
                    l, r = k, r + nmatch
                    z[k] = r - k + 1
        return z





    def big_l_prime_array(self, pattern, n):

        lp = [0] * len(pattern)
        for j in range(len(pattern) - 1):
            i = len(pattern) - n[j]
            if i < len(pattern):
                lp[i] = j + 1
        return lp






    def big_l_array(self, pattern, lp):

        l = [0] * len(pattern)
        l[1] = lp[1]
        for i in range(2, len(pattern)):
            l[i] = max(l[i - 1], lp[i])
        return l

    def small_l_prime_array(self, n):

        small_lp = [0] * len(n)
        for i in range(len(n)):
            if n[i] == i + 1:  # prefix matching a suffix
                small_lp[len(n) - i - 1] = i + 1
        for i in range(len(n) - 2, -1, -1):  # "smear" them out to the left
            if small_lp[i] == 0:
                small_lp[i] = small_lp[i + 1]
        return small_lp












class BoyerMooreAlgorithm:

        """
        BoyerMooreAlgorithm: This class implements
        the Boyer-Moore algorithm for exact pattern matching in a sequence.

        It uses the preprocessed pattern from BoyerMoorePreProcessing
        to efficiently find exact matches in the sequence.
        """
    
    def __init__(self, data):

        """
        Initializes the BoyerMooreAlgorithm class with the reference sequence data.
        """



        self.data = data
        
        
    def boyer_moore(self, pattern, p_bm, sequence):
        """
        Searches for exact matches of a given pattern within a reference sequence
        using the Boyer-Moore algorithm.

        Args:
            pattern (str): The pattern to search for in the reference sequence.
            p_bm (BmPreProcessing): The preprocessed pattern for Boyer-Moore algorithm.
            sequence (str): The reference sequence to search within.

        Returns:
            list: A list of starting positions in the sequence where the pattern is found.
        """
        i = 0
        num_iter = 0
        occurrences = []
        while i < len(sequence) - len(pattern) + 1:
            
            shift = 1  # Number of positions to shift the search window
            
            mismatched = False
            
            for j in range(len(pattern)-1, -1, -1):
                # Iterate through the pattern from end to start
                
                if pattern[j] != sequence[i + j]:
                   
                    skip_bp = p_bm.bad_character_rule(j, sequence[i + j])
                    # Calculate the shift based on bad character rule
                    skip_gs = p_bm.good_suffix_rule(j)
                    # Calculate the shift based on good suffix rule
                    
                    shift = max(skip_bp, skip_gs)  # Choose the maximum shift value
                    
                    mismatched = True
                    break
                    
            # If there is no mismatch, an occurrence is found
            if not mismatched:
                occurrences.append(i)
                skip_gs = p_bm.match_skip()
                shift = max(skip_gs, shift)
                    
            i += shift  # Shift the search window
            
        return occurrences

    
    


# Generate a random reference sequence (data)
data = ''.join([random.choice('ATCG') for _ in range(10000)])

# Define the pattern to search for in the reference sequence
pattern = 'TTAGCTT'

# Create an instance of the BoyerMooreAlgorithm class with the reference sequence (data)
model = BoyerMooreAlgorithm(data)

# Perform preprocessing on the pattern to create a preprocessed pattern (p_bm)
p_bm = BoyerMoorePreProcessing(pattern)

# Use the Boyer-Moore algorithm to find occurrences of the pattern in the reference sequence
hits = model.boyer_moore(pattern, p_bm, data)

# The 'hits' variable now contains a list of starting positions where the pattern is found
# in the reference sequence.

# Print the positions where the pattern was found
print("Occurrences of the pattern in the reference sequence:")
for pos in hits:
    print(f"Pattern found at position {pos}")


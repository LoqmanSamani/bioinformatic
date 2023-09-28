"""
The  NaiveAlgorithm  class is a  Python  implementation  for basic  DNA  sequence  analysis  tasks.
It is designed to handle input files in FASTA or FASTQ format and provides several functionalities,
including  reading  sequences from input files, calculating  base content statistics, finding exact
matches of  a  given DNA  pattern and  its reverse  complement  within a concatenated sequence, and 
generating  synthetic  DNA reads from a  reference genome. The class utilizes the Biopython library
for reading sequences and offers flexibility for various sequence analysis operations.
"""


from Bio import SeqIO  # Import SeqIO from Bio library
import collections as co
import random





class NaiveExactMachingAlgorithm:
   

    def __init__(self, file_path):
        
        """
        Initialize the NaiveAlgorithm class.

        Parameters:
        - file_path (str): The path to the input FASTA or FASTQ file.
        """
        
        self.file = file_path
        self.file_format_checker()  # Determine the file format
        self.sequences = self.file_reader()  # Dictionary to store sequences
        self.file_format = None

        # Dictionary for finding the complement of a DNA sequence
        self.complement_dict = {
            'A': 'T',
            'G': 'C',
            'T': 'A',
            'C': 'G'
        }
        
        
        
        

    def file_format_checker(self):
        
        """
        Check the format of the input file (FASTA or FASTQ).

        Returns:
        - str: The format of the input file ('FASTA', 'FASTQ', or 'UNKNOWN').
        """
        
        with open(self.file, 'r') as read_file:
            first_line = read_file.readline().rstrip()  # The first line of the file
            

            if first_line.startswith('>'):
                self.file_format = 'FASTA'
                
            elif first_line.startswith('@'):
                self.file_format = 'FASTQ'
                
            else:
                self.file_format = 'UNKNOWN'
                

            if self.file_format not in ['FASTA', 'FASTQ']:
                raise ValueError("This app can only handle FASTA or FASTQ files. The format of the input file is not recognized.")

                
                
                
                
                
                
    def file_reader(self):
        
        """
        Read sequences from the input file.

        Returns:
        - dict: A dictionary where keys are sequence IDs and values are tuples (length, sequence).
        """
        
        sequences = {}

        
        with open(self.file, 'r') as handle:
            
            if self.file_format == 'FASTA':
                
                for seq in SeqIO.parse(handle, 'fasta'):
                    sequences[seq.id] = [len(seq.seq), seq.seq]
                    
            elif self.file_format == 'FASTQ':
                
                for seq in SeqIO.parse(handle, 'fastq'):
                    sequences[seq.id] = [len(seq.seq), seq.seq]
                    
            else:
                raise ValueError("Please make sure that the format of the input file is FASTA or FASTQ")

                
        return sequences
    
    
    
    
    

    def base_content(self):
        
        """
        Calculate the number and percentage of each base in sequences.

        Parameters:
        - sequences (dict): Dictionary of sequences.

        Returns:
        - dict: A dictionary with base counts and percentages.
        """
        
        count = co.Counter()  # Dictionary to count bases

        
        for sequence in self.sequences.values():
            count.update(sequence[1])

            
        percentages = {}
        
        for key, val in count.items():
            
            if val > 0:  # Avoid zero division
                percentages[key] = val / sum(count.values())

                
        return count, percentages
    
    
    
    
    

    def concatenate_sequences(self):
        
        """
        Concatenate sequences into a single string.

        Parameters:
        - sequences (dict): Dictionary of sequences.

        Returns:
        - str: Concatenated sequence.
        """
        
        concatenated_sequences = ''
        
        for val in self.sequences.values():
            concatenated_sequences += val[1]
            
            
        return concatenated_sequences
    
    
    
    
    

    def reverse_complement(self, pattern):
        
        """
        Compute the reverse complement of a pattern DNA.

        Parameters:
        - concatenated_sequences (str): Input DNA sequence.

        Returns:
        - str: Reverse complement of the input sequence.
        """
        
        pattern = pattern.upper()
        
        r_pattern = ''

        for char in pattern:
            r_pattern = self.complement_dict[char] + r_pattern

            
        return r_pattern
    
    
    
    
    
    

    def naive(self, concatenated_sequences, pattern, r_pattern):
        
        """
        Find exact matches of a pattern and its reverse complement 
        in a concatenated sequence.

        Parameters:
        - concatenated_sequences (str): Concatenated DNA sequence.
        - pattern (str): Pattern to search for.
        - r_pattern (str): reverse complement of the pattern sequence

        Returns:
        - list: List of starting indices of exact matches.
        """
        
        occurrences = []  # List to store starting indices of matches
        
        #r_sequence = self.reverse_complement(concatenated_sequences)

        for i in range(len(concatenated_sequences) - len(pattern) + 1):
            is_match = True

            
            for j in range(len(pattern)):
                
                if concatenated_sequences[i + j] != pattern[j]:
                    is_match = False
                    break

                    
            if is_match:
                occurrences.append(i)
                
        if pattern != r_pattern: # check if pattern and its reverse are not the same
            
            for i in range(len(concatenated_sequences) - len(r_pattern) + 1):
                
                is_match = True

                for j in range(len(pattern)):
                
                    if concatenated_sequences[i + j] != r_pattern[j]:
                        is_match = False
                        break

                    
                if is_match:
                    occurrences.append(i)
            
        return occurrences
    
    
    
    
    

    def generate_reads(self, concatenated_sequences, num_reads, len_reads):
        
        """
        Generate synthetic DNA reads from a given genome.

        Parameters:
        - concatenated_sequences (str): Concatenated DNA sequence.
        - num_reads (int): The number of reads to generate.
        - len_reads (int): The desired length of each read.

        Returns:
        - list: A list containing the generated DNA reads.
        """
        
        reads = []  # Initialize an empty list to store generated reads

        for i in range(num_reads):
            
            random_start = random.randint(0, len(concatenated_sequences) - len_reads)
            generated_read = concatenated_sequences[random_start:random_start + len_reads]
            reads.append(str(generated_read))

            
        return reads
    
    
    
    
    
    
    
    
# Usage example:

path = 'address_of_a_file'
pattern = 'AGCTGGGTCANN'


model = NaiveExactMachingAlgorithm(path)

base_count, base_percentages = model.base_content()

concatenated_seq = model.concatenate_sequences()

r_pattern = model.reverse_complement(pattern)

naive_matches = model.naive(concatenated_seq, pattern, r_pattern)

generated_reads = model.generate_reads(concatenated_seq, num_reads=200, len_reads=20)


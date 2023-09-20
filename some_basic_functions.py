#!/usr/bin/env python
# coding: utf-8

# In[64]:


from Bio.Data import CodonTable
import numpy as np 


class SequenceAnalyser:
    
    def __init__(self):
        
        # Initialize the class with dictionaries and bases for DNA sequence analysis
        self.bases = ['A', 'G', 'C', 'T']
        
        # Dictionary for DNA to RNA translation
        self.translate_dict = {
            'A':'U',
            'G':'C',
            'T':'A',
            'C':'G'
        }
        
        # Dictionary for finding the complement of a DNA sequence
        self.complement_dict = {
            'A':'T',
            'G':'C',
            'T':'A',
            'C':'G'
        }
        
        # Dictionary to store codons and their corresponding amino acids.
        self.codon_dict = {}  
        
        # Initialize the genetic code for translation
        genetic_code = CodonTable.unambiguous_dna_by_id[1]
        for codon, amino_acid in genetic_code.forward_table.items():
            self.codon_dict[codon] = amino_acid
            
            
            
            
            
        
    def seq_cleaner(self, seq):
        """
        Clean a DNA sequence by removing non-standard
        characters and return the cleaned sequence.
        """
        
        cleaned_seq = ''.join([char for char in seq if char in self.bases])

        return cleaned_seq
    
    
    
    
        
        
    def gc_percentage(self, seq):
        """
        Calculate and return the GC percentage of a DNA sequence.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        g_counter = cleaned_seq.count('G')
        c_counter = cleaned_seq.count('C')
        gc_percentage = ((g_counter + c_counter) / len(cleaned_seq)) * 100
        
        return round(gc_percentage, 2)
    
    
    
    
    
    def sec_transcriptor(self, seq):
        """
        Transcribe a DNA sequence into RNA and return the transcribed sequence.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        translated_seq = ''.join([self.translate_dict[char] for char in cleaned_seq])
        
        return translated_seq
    
    
    
    
    
    def seq_translator(self, seq):
        """
        Translate a DNA sequence into a polypeptide (protein) sequence.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        codons = [cleaned_seq[i:i+3] for i in range(len(cleaned_seq))]
        polypeptide = ''.join([self.codon_dict[codon] if codon in self.codon_dict else '?' for codon in codons ])
        
        return polypeptide
    
    
    
    
    
    
    def seq_complement(self, seq):
        """
        Find the complement of a DNA sequence and return the complemented sequence.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        seq_complement = ''.join([self.complement_dict[base] for base in cleaned_seq])
        
        return seq_complement
    
    
    
    
    
    def reverse_complement(self, seq):
        """
        Find the reverse complement of a DNA sequence and return it.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        reverse_complement = ''.join([self.complement_dict[base] for base in cleaned_seq[::-1]])
        
        return reverse_complement
    
    
    
    
    
    def length_seq(self, seq):
        """
        Calculate and return the length of a DNA sequence.
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        
        return len(cleaned_seq)
    
    
    
    
    def numcleotide_frequency(self, seq):
        """
        Function to calculate the nucleotide frequencies in a DNA sequence
        """
        
        frequencies = {}
        cleaned_seq = self.seq_cleaner(seq)
        
        a_count = cleaned_seq.count('A')
        g_count = cleaned_seq.count('G')
        c_count = cleaned_seq.count('C')
        t_count = cleaned_seq.count('T')
        
        frequencies['A'] = a_count / len(cleaned_seq)
        frequencies['G'] = g_count / len(cleaned_seq)
        frequencies['C'] = c_count / len(cleaned_seq)
        frequencies['T'] = t_count / len(cleaned_seq)
        
        return frequencies
    
    
    
    def gc_skew(self, seq):
        """
        Function to calculate the GC skew of a DNA sequence
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        
        g_counter = cleaned_seq.count('G')
        c_counter = cleaned_seq.count('C')

        gc_skew = (g_counter - c_counter) / (g_counter + c_counter)
        
        return gc_skew
    
    
    
    
    def transcription(self, seq):
        """
        Function to transcribe a DNA sequence into an RNA sequence
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        transcribed_seq = cleaned_seq.replace('T', 'U')
    
        return transcribed_seq
    
    
    
    
    def reverse_transcription(self, seq):
        """
        Function to reverse transcribe an RNA sequence into a DNA sequence
        """
        
        cleaned_seq = self.seq_cleaner(seq)
        reverse_transcribed_seq = cleaned_seq.replace('U', 'T')
        
        return reverse_transcribed_seq
    
    
    
    
    
    def codon_count(self, seq):
        """
        Function to count the occurrences of each codon in a DNA sequence
        """
        
        codon_count = {}
        
        cleaned_seq = self.seq_cleaner(seq)
        
        codons = [cleaned_seq[i:i+3] for i in range(len(cleaned_seq))]
        
        for codon in codons:
            if codon not in codon_count:
                codon_count[codon] = 1
            else:
                codon_count[codon] += 1
                
        return codon_count
    
    
    
        
        
        
        

# Generate a random DNA sequence
random_sequence = ''.join(np.random.choice(['A', 'G', 'C', 'T'], 100))

# Create an instance of SequenceAnalyser
sequence_analyser = SequenceAnalyser()

# Clean the DNA sequence
cleaned_sequence = sequence_analyser.seq_cleaner(random_sequence)
print("Cleaned Sequence:", cleaned_sequence)

# Calculate the GC percentage of the sequence
gc_percentage = sequence_analyser.gc_percentage(cleaned_sequence)
print("GC Percentage:", gc_percentage, "%")

# Transcribe the DNA sequence to RNA
transcribed_sequence = sequence_analyser.sec_transcriptor(cleaned_sequence)
print("Transcribed Sequence (RNA):", transcribed_sequence)

# Translate the RNA sequence to a polypeptide sequence
polypeptide_sequence = sequence_analyser.seq_translator(transcribed_sequence)
print("Polypeptide Sequence:", polypeptide_sequence)

# Find the complement of the DNA sequence
sequence_complement = sequence_analyser.seq_complement(cleaned_sequence)
print("Sequence Complement:", sequence_complement)

# Find the reverse complement of the DNA sequence
reverse_complement_sequence = sequence_analyser.reverse_complement(cleaned_sequence)
print("Reverse Complement:", reverse_complement_sequence)

# Calculate the length of the DNA sequence
sequence_length = sequence_analyser.length_seq(cleaned_sequence)
print("Sequence Length:", sequence_length)

# Calculate the nucleotide frequency of the DNA sequence
nucleotide_frequency = sequence_analyser.numcleotide_frequency(cleaned_sequence)
print("Nucleotide Frequency:", nucleotide_frequency)

# Calculate the GC skew of the DNA sequence
gc_skew = sequence_analyser.gc_skew(cleaned_sequence)
print("GC Skew:", gc_skew)

# Transcribe the DNA sequence to RNA
rna_sequence = sequence_analyser.transcription(cleaned_sequence)
print("Transcribed Sequence (RNA):", rna_sequence)

# Reverse transcribe the RNA sequence to DNA
reverse_transcribed_sequence = sequence_analyser.reverse_transcription(rna_sequence)
print("Reverse Transcribed Sequence (DNA):", reverse_transcribed_sequence)

# Count the occurrences of each codon in the DNA sequence
codon_counts = sequence_analyser.codon_count(cleaned_sequence)
print("Codon Counts:", codon_counts)

              


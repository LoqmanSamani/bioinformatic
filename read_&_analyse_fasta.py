


from Bio import SeqIO as s


class FastaReader:
    
    def __init__(self, file):
        
        self.file = file
        self.sequences = self.fasta_reader()
        
        
    def fasta_reader(self):
        
        sequences = {}
        
        with open(self.file, 'r') as handle:
            
            for seq in s.parse(handle, 'fasta'):
                
                sequences[seq.id] = [len(seq.seq), seq.seq]
                
        self.sequences = sequences
        
        
        
    def num_sequences(self):
        
        return len(self.sequences)




    def max_min_seqs(self):
        
        longest_seq = None
        shortest_seq = None

        for key, val in self.sequences.items():
            
            if longest_seq is None or val[0] > longest_seq[1][0]:
                
                longest_seq = (key, val)
                
            if shortest_seq is None or val[0] < shortest_seq[1][0]:
                
                shortest_seq = (key, val)
                

        return longest_seq, shortest_seq
    
    
    
    
    
    def orf_finder(self):
        
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        
        orfs_in_each_seq = {}
        
        for key, val in self.sequences.items():
            
            orfs = []
            
            forward = (key, val[1])
            reverse = (key, val[1][::-1])
            
            for frame in range(3):
        
                forward_frame = forward[1][frame:]
                f_codons = [forward[i:i+3] for i in range(0, len(forward_frame), 3)]
            
                f_start_index = None
                f_stop_index = None
            
                for index, codon in enumerate(f_codons):
                
                    if codon == start_codon:
                        f_start_index = index
                    
                    elif codon in stop_codons:
                        f_stop_index = index
                    
                        if f_start_index is not None:
                            orfs.append(forward_frame[f_start_index:f_stop_index+1])
                    
                            f_start_index = None
                            f_stop_index = None
                            
                            
                reverse_frame = reverse[1][frame:]
                r_codons = [reverse[i:i+3] for i in range(0, len(reverse_frame), 3)]
                
                r_start_index = None
                r_stop_index = None  
            
                for index, codon in enumerate(r_codons):
                
                    if codon == start_codon:
                        r_start_index = index
                    
                    elif codon in stop_codons:
                        r_stop_index = index
                    
                        if r_start_index is not None:
                            orfs.append(reverse_frame[r_start_index:r_stop_index+1])
                    
                            r_start_index = None
                            r_stop_index = None
                        
    
            # save each orf with its length in a tuple
            orfs = [(len(orf)*3, ''.join(orf)) for orf in orfs]
        
            orfs_in_each_seq[key] = orfs
    
        return orfs_in_each_seq
    
    


    def repeat_finder(self, len_repeat):
        
        repeat_in_seq = {}
        
        for key, val in self.sequences.items():
            
            repeats = {}
            max_repeat = ('', 0)
            
            for i in range(len(val[1]) - len_repeat + 1):
                
                repeat = key[1][i:i + len_repeat]
                repeats[repeat] = repeats.get(repeat, 0) + 1

                if repeats[repeat] > max_repeat[1]:
                    max_repeat = (repeat, repeats[repeat])
                    
                    
            repeat_in_seq[key] = (max_repeat, repeats)
            

        return repeat_in_seq


        
        
# Instantiate the FastaReader class with a FASTA file path
fasta_reader = FastaReader('example.fasta')


# Call the fasta_reader method to read and parse the FASTA file
fasta_reader.fasta_reader()


# Get the number of sequences in the file
num_seqs = fasta_reader.num_sequences()
print(f"There are {num_seqs} sequences in the file.")



# Find the longest and shortest sequences
longest_seq, shortest_seq = fasta_reader.max_min_seqs()
print(f"The longest sequence is '{longest_seq[0]}' with a length of {longest_seq[1][0]} bases.")
print(f"The shortest sequence is '{shortest_seq[0]}' with a length of {shortest_seq[1][0]} bases.")



# Find ORFs in the sequences
orf_results = fasta_reader.orf_finder()
for seq_id, orfs in orf_results.items():
    print(f"ORFs in sequence '{seq_id}':")
    for orf in orfs:
        length, sequence = orf
        print(f"Length: {length}, Sequence: {sequence}")
        
        
        

# Find repeats in the sequences
len_repeat = 4  # Define the length of the repeat you want to search for
repeat_results = fasta_reader.repeat_finder(len_repeat)
for seq_id, (max_repeat, repeats) in repeat_results.items():
    print(f"Repeats in sequence '{seq_id}':")
    print(f"Most frequent repeat: '{max_repeat[0]}' with a count of {max_repeat[1]}")
    print("All repeats:")
    for repeat, count in repeats.items():
        print(f"Repeat: '{repeat}', Count: {count}")
    
                        


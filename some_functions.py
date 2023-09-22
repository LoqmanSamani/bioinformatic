

import random

def longest_common_prefix(seq1, seq2):
    
    common_prefix = ''
    i = 0
    while i < len(seq1) and i < len(seq2):
        if seq1[i] == seq2[i]:
            common_prefix += seq1[i]
            i +=1
            
        else:
            break
            
    if len(common_prefix) == 0:
        
        return "There is no common residue in seq2 & seq1."
    
    return common_prefix, i



seq1 = ''.join([random.choice('ACTG') for _ in range(1000)])
seq2 = ''.join([random.choice('ACTG') for _ in range(1000)])

seq = longest_common_prefix(seq1, seq2)







def match_sequences(seq1, seq2):
    
    if len(seq1) == len(seq2) and seq1 == seq2:
        
        return 'The input sequences are the same.'
    
    elif len(seq1) == len(seq2):
        
        return 'The input sequences are the same length but not identical.'
    else:
        return 'The input sequences are not the same length'
    
    
seq1 = ''.join([random.choice('ACTG') for _ in range(1000)])
seq2 = ''.join([random.choice('ACTG') for _ in range(1000)])

match_sequences = match_sequences(seq1, seq2)








def read_fastq(file_path):
    
    seqs = []
    scores = []
    
    reading_sequence = False 
    
    with open(file_path, 'r') as file:
        
        for line in file:
            
            if line.startswith('@') and not reading_sequence:
                
                reading_sequence = True
                
            elif reading_sequence:
                
                seqs.append(line.rstrip())
                
                reading_sequence = False
                
            elif line.startswith('+'):
                
                reading_sequence = True
                
            elif reading_sequence:
                
                scores.append(line.rstrip())
    
    seqs_and_scores = list(zip(seqs, scores))
    
    
    return seqs_and_scores    


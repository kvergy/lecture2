#!/usr/bin/env python
# coding: utf-8

# In[8]:


import sys
import os

dict_dna = {'A': 135, 'T' : 126, 'C' : 111, 'G' : 151}
dict_2chain = {'A': 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
dict_rna = {'A': 135, 'U' : 112, 'C' : 111, 'G' : 151}
class origin:
    dicti = dict_dna
    def __init__(self, seq):
        self.chain = seq
        self.len = len(seq)
        
    def weight(self):
        weight = 0
        for residue in self.chain:
            weight += self.dicti[residue]
        return weight  
   
            
class Class_dna(origin):
    dicti = dict_dna

        
    def transcription(self):
        rna = ''
        for residue in self.chain:
            if residue == 'T':
                rna += 'U'
            else:
                rna += residue
        return Class_rna(rna)    
    def second_chain(self):
        chain2 = ''
        for residue in self.chain:
            chain2 += dict_2chain[residue]
        return Class_dna(chain2)
class Class_rna(origin):
    dicti = dict_rna
         
    def translation(self):
        translation = { 
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',                  
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W', 
            } 
        protein = ''
        rna = self.chain
        i = 0
        for i in range(0,len(rna)-3,3):
            protein += translation[rna[i:i + 3]]
        return Class_protein(protein) 

class Class_protein(origin):
    dicti = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

     

def find_orf(seq):
    end = ['TAG', 'TAA', 'TGA']
    i = 0
    dna = ''
    orf = ''
    j = 0
    k = 0
    for k in range(3):
        orf = ''
        for i in range(k, len(seq)-3):
            if seq[i:i + 3] == 'ATG':
                j = i + 3
                while j < len(seq) - 3 and seq[j:j + 3] not in end:
                    orf += seq[j:j+3]
                    j += 3
                if len(orf) > len(dna) and len(orf) > 120:
                    dna = orf

    
    return Class_dna(dna)

#filename = stdin
seq=''

for path in sys.stdin:
    with open(path.strip()) as inf:
        file = inf.read().split('\n')
        for row in file[1:]:
            seq += row
    chain1 = Class_dna(seq)
    chain2 = chain1.second_chain()
    orf1 = find_orf(chain1.chain)
    orf2 = find_orf(chain2.chain)
    dna = orf2
    if orf1.len >= orf2.len:
        dna = orf1

    rna = dna.transcription()
    protein = rna.translation()
    print("Matrix DNA:{}".format(dna.chain))
    print("Matrix DNA weight:{}".format(dna.weight()))
    print("RNA:{}".format(rna.chain))
    print("RNA weight:{}".format(rna.weight()))
    print("Protein:{}".format(protein.chain))
    print("Protein weight:{}".format(protein.weight()))













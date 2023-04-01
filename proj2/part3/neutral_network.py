from itertools import combinations, permutations
from scipy.special import binom
import numpy as np


table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }




class BasicNeutralNetwork:
    def __init__(self, n_mutations) -> None:
        self.mutations = {}
        self.n_mutations = n_mutations
        self._build()
        
    def _build(self):
        for codon in table:
            syn = []
            nonsyn = []
            for cmp_codon in table:
                if sum([(codon[i] != cmp_codon[i]) for i in range(3)]) == self.n_mutations:
                    if table[codon] == table[cmp_codon]:
                        syn.append(cmp_codon)
                    else:
                        nonsyn.append(cmp_codon)
            self.mutations[codon] = (syn, nonsyn)
            
    def get_mutation(self, codon):
        return self.mutations[codon]



class NeutralNetworkCounter:
    def __init__(self, seq):
        self.seq = seq
        self.size = len(seq)
        
    def build(self, n_mutations):
        self.neighborhood = [BasicNeutralNetwork(i+1) for i in range(n_mutations)]
        self.syn = []
        for n in self.neighborhood:
            self.syn.append(np.array([len(n.get_mutation(self.seq[i:i+3])[0]) for i in range(0, len(self.seq), 3)]))
        
        syn_cnt, nonsyn_cnt = 0, 0
        syn_cnt += np.sum(self.syn[n_mutations-1])
        if n_mutations == 2:
            syn_cnt += (np.power(np.sum(self.syn[0]), 2) - np.sum(np.power(self.syn[0], 2))) / 2
        elif n_mutations == 3:
            syn_cnt += np.sum(self.syn[0]) * np.sum(self.syn[1]) - np.dot(self.syn[0], self.syn[1])
            syn_cnt += (np.power(np.sum(self.syn[0]), 3) - 3*np.dot(self.syn[0], np.power(self.syn[0], 2)) + 2*np.sum(np.power(self.syn[0], 3))) / 6
        
        nonsyn_cnt = binom(self.size, n_mutations) * 3**n_mutations - syn_cnt
            
        return int(syn_cnt), int(nonsyn_cnt)
                
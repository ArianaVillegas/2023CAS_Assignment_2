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
        self.neighborhood = [BasicNeutralNetwork(i+1) for i in range(2)]
        self.syn = []
        for n in self.neighborhood:
            self.syn.append(np.array([len(n.get_mutation(self.seq[i:i+3])[0]) for i in range(0, len(self.seq), 3)]))
            
        uniq = [[1, sum(self.syn[0])], [1, sum(self.syn[1])]]
        for idx in range(2):
            for i in range(2, 1 + n_mutations//(idx+1)):
                cur = uniq[idx][i-1] * uniq[idx][1] 
                for k in range(i-1):
                    cur -= (-1)**k * np.sum(np.power(self.syn[idx], 2+k)) * uniq[idx][i-2-k]
                cur //= i
                uniq[idx].append(cur)
                
        syn_cnt = 0
        for i in range(0, n_mutations + 1, 2):
            cur_cnt = uniq[0][n_mutations - i] * uniq[1][i//2] 
            syn_cnt += cur_cnt
        
        nonsyn_cnt = binom(self.size, n_mutations) * 3**n_mutations - syn_cnt
            
        return int(syn_cnt), int(nonsyn_cnt)


def main(inpath, n):
    f = open(inpath, "r")
    seq = ''.join(f.readlines()).replace('\n', '') + 'A'
    f.close()
    
    neutral_network = NeutralNetworkCounter(seq)
    syn, nonsyn = neutral_network.build(n)
    
    print('Synonymous:', syn)
    print('Non-synonymous', nonsyn)
    print('Total', syn + nonsyn)
    

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
            "Get the approximate number of synonymous, non-synonymous, and the total number of mutations")
    parser.add_argument("--file", type=str, default="part1/seq.txt",
            help="The relative path to the antigen escape calculator data")
    parser.add_argument("--n", type=int, default=1,
            help="Number of neutral network")
    args = parser.parse_args()
    main(args.file, args.n)

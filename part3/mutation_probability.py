# Wrtten in collaboration with Christopher Leap

from neutral_network import table, BasicNeutralNetwork
import numpy as np
import pandas as pd

# Change number from 1, 2, or 3, for respective n-point strict mutation; leave as 1 for non-strict n-point mutation
base_muts = BasicNeutralNetwork(1)
base_muts._build()

#for key, value in base_muts.mutations.items():
    #print(codon+":", muts)

idx2aa = list(set(table.values()))
idx2aa.sort()
#print(idx2aa, len(idx2aa))

aa2idx = {aa: i for i, aa in enumerate(idx2aa)}
#print(aa2idx)

count_table = np.zeros((len(idx2aa), len(idx2aa)))

for condonA, muts in base_muts.mutations.items():
    amino = table[condonA]
    row = aa2idx[amino]
    count_table[row][row] += len(muts[0])

    for codonB in muts[1]:
        amino = table[codonB]
        col = aa2idx[amino]
        count_table[row][col] += 1

for i, row in enumerate(count_table):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{int(val):2d} ", end="")
    print()

trans_tablenp = count_table/np.sum(count_table, axis=0)

for i, row in enumerate(trans_tablenp):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

# Uncomment for non-strict mutation
'''
trans_table2p = np.dot(trans_tablenp, trans_tablenp)

for i, row in enumerate(trans_table2p):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

trans_table3p = np.dot(trans_table2p, trans_tablenp)

for i, row in enumerate(trans_table3p):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()
'''

inpath = "../part2/aamut_fitness_all.csv"
gene = "S"
df = pd.read_csv(inpath)
df["keep"] = df.apply(lambda row: row["gene"] == gene, axis=1)
df = df[df["keep"]]

#print(df)
#print(df["clade_founder_aa"].value_counts())
values = df["clade_founder_aa"].value_counts()

gene_count = np.zeros(len(idx2aa))

for i, amino in enumerate(idx2aa):
    if amino not in values:
        continue
    else:
        gene_count[i] = values[amino]
print(gene_count)

gene_dist = gene_count/np.sum(gene_count)
print(gene_dist)

gene_transnp = np.dot(gene_dist, trans_tablenp)
print(gene_transnp)

# For strict 1- 2- or 3-point mutation, change the argument number on line 17 respectively
expected_strict_gene_transnp = gene_transnp*np.sum(gene_count)
print(expected_strict_gene_transnp)

# Uncomment for non-strict mutation
'''
nonstrict_gene_trans2p = np.dot(gene_dist, trans_table2p)
print(nonstrict_gene_trans2p)

nonstrict_gene_trans3p = np.dot(gene_dist, trans_table3p)
print(nonstrict_gene_trans3p)
'''

synnonsyn_count_table = np.zeros((len(idx2aa), 2))

for condon, muts in base_muts.mutations.items():
    amino = table[condon]
    row = aa2idx[amino]
    synnonsyn_count_table[row][0] += len(muts[0])
    synnonsyn_count_table[row][1] += len(muts[1])

for i, row in enumerate(synnonsyn_count_table):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{int(val):2d} ", end="")
    print()

synonsyn_trans_tablenp = (synnonsyn_count_table.T/np.sum(synnonsyn_count_table.T, axis=0)).T

for i, row in enumerate(synonsyn_trans_tablenp):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

synonsyn_gene_trans_table1p = np.dot(gene_dist, synonsyn_trans_tablenp)
print(synonsyn_gene_trans_table1p)

# For 2- or 3-point mutation, change the argument number on line 17 to 1 or 2, respectively
synonsyn_gene_trans_tablenp = np.dot(gene_transnp, synonsyn_trans_tablenp)
print(synonsyn_gene_trans_tablenp)

#temp2aa = np.dot(gene_dist, trans_tablenp)
#temp2syn = np.dot(gene_dist, trans_table1p2)

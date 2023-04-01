from neutral_network import table, BasicNeutralNetwork
import numpy as np


filename = 'seq.txt'

f = open(filename, "r")
seq = ''.join(f.readlines()).replace('\n', '') + 'A'
f.close()

print(len(seq))


base_muts = BasicNeutralNetwork(1)
base_muts._build()

#for codon, muts in base_muts.mutations.items():
    #print(codon+":", muts)


idx2aa = list(set(table.values()))
idx2aa.sort()
print(idx2aa, len(idx2aa))

aa2idx = {aa: i for i, aa in enumerate(idx2aa)}
print(aa2idx)

count_table = np.zeros((len(idx2aa), len(idx2aa)))

for condonA, mutsA in base_muts.mutations.items():
    amino = table[condonA]
    row = aa2idx[amino]
    count_table[row][row] += len(mutsA[0])

    for codonB in mutsA[1]:
        amino = table[codonB]
        col = aa2idx[amino]
        count_table[row][col] += 1

for i, row in enumerate(count_table):
    print(f"{i:2d} ", end="")
    for val in row:
        print(f"{int(val):2d} ", end="")
    print()

trans_table1p = count_table/np.sum(count_table, axis=0)

for i, row in enumerate(trans_table1p):
    print(f"{i:2d} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

'''
trans_table2p = np.dot(trans_table, trans_table)

for i, row in enumerate(trans_table2p):
    print(f"{i:2d} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

trans_table3p = np.dot(trans_table2p, trans_table)

for i, row in enumerate(trans_table3p):
    print(f"{i:2d} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()
'''

#OLD CODE
'''
syn_count = 0
nonsyn_count = 0
temp = []
temp2 = []
temp3 = []
checked = []

for amino, acid in base_muts.mutations.items():
    if amino in checked:
        continue
    syn_count += len(base_muts.mutations[amino][0])
    nonsyn_count += len(base_muts.mutations[amino][1])
    for amino2 in base_muts.mutations.keys():
        if amino2 == amino:
            continue
        if len(acid[0]) != 0:
            if len(acid[0]) == 3:
                if acid[0][0] == amino2 or acid[0][1] == amino2 or acid[0][2] == amino2:
                    syn_count += 3
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
            elif len(acid[0]) == 2:
                if acid[0][0] == amino2 or acid[0][1] == amino2:
                    syn_count += 2
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
            elif len(acid[0]) == 1:
                if acid[0] == amino2:
                    syn_count += 1
                    nonsyn_count += len(acid[1])
                    checked.append(amino2)
        else:
            continue
    print(syn_count)
    temp.append(syn_count)
    temp2.append(nonsyn_count+syn_count)
    temp3.append([syn_count, nonsyn_count+syn_count])
    if len(temp) != 0:
        if syn_count != temp[-1]:
            temp.append(syn_count)
    else:
        temp.append(syn_count)

    #syn_count = 0
    #nonsyn_count = 0

for amino in table.values():
    for amino2 in table.values():
        if amino == amino2:
            syn_count += 1
        else:
            nonsyn_count += 1
    print(syn_count, nonsyn_count)
    syn_count = 0
    nonsyn_count = 0
'''



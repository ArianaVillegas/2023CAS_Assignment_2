from neutral_network import NeutralNetworkCounter, BasicNeutralNetwork
import numpy as np


filename = 'seq.txt'

f = open(filename, "r")
seq = ''.join(f.readlines()).replace('\n', '') + 'A'
f.close()

print(len(seq))

# seq = 'CGCGATACATGAATC'
basic = BasicNeutralNetwork(1)
count = np.array([[len(basic.get_mutation(codon)[0]), len(basic.get_mutation(codon)[1])] for codon in basic.mutations])[:-1]
print(basic.mutations)
# print(np.sum(count, axis=0))
# print(len(count))
# print(138/(429+138))

# neutral_network = NeutralNetworkCounter(seq)
# syn, nonsyn = neutral_network.build(2)
# print(syn, nonsyn)
# print(syn + nonsyn)
# print(120325256403408 * 39872)
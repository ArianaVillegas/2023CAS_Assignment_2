from neutral_network import NeutralNetworkCounter, BasicNeutralNetwork


filename = 'seq.txt'

f = open(filename, "r")
seq = ''.join(f.readlines()).replace('\n', '') + 'A'
f.close()

print(len(seq))

# basic_neutral_network = BasicNeutralNetwork()
# basic_neutral_network.build(2)

neutral_network = NeutralNetworkCounter(seq)
syn, nonsyn = neutral_network.build(1)
print(syn, nonsyn)
print(syn + nonsyn)
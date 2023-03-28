# Define the reference sequence
filename = 'seq.txt'

f = open(filename, "r")
seq = ''.join(f.readlines()).replace('\n', '') + 'A'
f.close()

print(len(seq))

# Define a function to generate all possible single-step mutants
def generate_single_step_mutants(seq):
    # List of possible nucleotides
    nucleotides = ["A", "C", "G", "T"]
    
    # Loop through each position in the sequence
    for i in range(len(seq)):
        # Loop through each possible nucleotide substitution
        for nuc in nucleotides:
            # Generate the mutated sequence
            mutant_seq = seq[:i] + nuc + seq[i+1:]
            yield mutant_seq

# Define a function to generate all possible n-step mutants
def generate_n_step_mutants(seq, n):
    # Start with the single-step mutants of the reference sequence
    mutants = set(generate_single_step_mutants(seq))
    
    # Generate all possible n-step mutants
    for i in range(n-1):
        new_mutants = set()
        for mutant in mutants:
            new_mutants.update(generate_single_step_mutants(mutant))
        mutants.update(new_mutants)
    
    # Remove the reference sequence from the set of mutants
    mutants.remove(seq)
    
    return mutants

# Generate all possible 3-step mutants of the reference sequence
mutants = generate_n_step_mutants(seq, 3)

# Print the number of mutants
print(len(mutants))

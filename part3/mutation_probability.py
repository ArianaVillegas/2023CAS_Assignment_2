# Collaboration by Emmanuel Ohiri Christopher Leap

from neutral_network import table, BasicNeutralNetwork
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TABLE_DIR = "../tables"
STRICT = False
count = 1


idx2aa = list(set(table.values()))
idx2aa.sort()
aa2idx = {aa: i for i, aa in enumerate(idx2aa)}

def write_trans_table(path, table, dtype, caption, label):
    with open(path, "w") as file:
        file.write("\\begin{table*}[]\n\t\centering\n\t\\begin{tabular}")
        file.write("{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c}\n")

        # Write the headers
        file.write("\t\t")
        for i in range(len(table)):
            if idx2aa[i] == "_":
                file.write(" & \\_")
            else:
                file.write(f" & {idx2aa[i]}")
        file.write(" \\\\\n")

        for i, row in enumerate(table):
            file.write("\t\t\\hline\n")
            if idx2aa[i] == "_":
                file.write("\t\t\\_")
            else:
                file.write(f"\t\t{idx2aa[i]}")
            for val in row:
                if dtype == int:
                    file.write(f" & {int(val):2d}")
                elif dtype == float:
                    file.write(f" & {val:.2f}")
            file.write(" \\\\\n")

        file.write("\t\\end{tabular}\n")
        file.write("\t\\caption{" + caption + "}\n")
        file.write("\t\\label{tab:" + label + "}\n")
        file.write("\\end{table*}\n")

def make_trans_table(num_hops):
    # Change number from 1, 2, or 3, for respective n-point strict mutation; leave as 1 for non-strict n-point mutation
    base_muts = BasicNeutralNetwork(num_hops)
    base_muts._build()

    count_table = np.zeros((len(idx2aa), len(idx2aa)))

    for condonA, muts in base_muts.mutations.items():
        amino = table[condonA]
        row = aa2idx[amino]
        count_table[row][row] += len(muts[0])

        for codonB in muts[1]:
            amino = table[codonB]
            col = aa2idx[amino]
            count_table[row][col] += 1

    write_trans_table(f"{TABLE_DIR}/count_table_{num_hops}_pt.tex",
            count_table,
            int,
            caption=f"A table of the number of {num_hops} single-point mutations that each amino acid can go through to get from one amino acid to another. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
            label=f"count_table_{num_hops}_pt")

    trans_table = count_table/np.sum(count_table, axis=0)
    return trans_table

# Single-point transition table
trans_tablenp = make_trans_table(1)
write_trans_table(f"{TABLE_DIR}/trans_tablenp.tex",
        trans_tablenp,
        float,
        caption="A transition table for 1 single-point mutation. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
        label="trans_tablenp")

for i, row in enumerate(trans_tablenp):
    print(f"{idx2aa[i]} ", end="")
    for val in row:
        print(f"{val:.2f} ", end="")
    print()

# Non-strict transition tables
trans_table2p = np.dot(trans_tablenp, trans_tablenp)
write_trans_table(f"{TABLE_DIR}/trans_table2p_non_strict.tex",
        trans_table2p,
        float,
        caption="A transition table for 2 non-strict single-point mutations. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
        label="trans_table2p_non_strict")

trans_table3p = np.dot(trans_table2p, trans_tablenp)
write_trans_table(f"{TABLE_DIR}/trans_table3p_non_strict.tex",
        trans_table3p,
        float,
        caption="A transition table for 3 non-strict single-point mutations. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
        label="trans_table3p_non_strict")

# Strict transition tables
trans_table2p = make_trans_table(2)
write_trans_table(f"{TABLE_DIR}/trans_table2p_strict.tex",
        trans_table2p,
        float,
        caption="A transition table for 2 strict single-point mutations. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
        label="trans_table2p_strict")

trans_table3p = make_trans_table(3)
write_trans_table(f"{TABLE_DIR}/trans_table3p_strict.tex",
        trans_table3p,
        float,
        caption="A transition table for 3 non-strict single-point mutations. The row corresponds to the starting amino acid and the column corresponds to the final amino acid.",
        label="trans_table3p_strict")

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

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 6

def plot_dist(dist, title, path):
    fig, axs = plt.subplots(figsize=(3.5,3.5))
    axs.bar(idx2aa, gene_dist, color="white", edgecolor="black")
    axs.set_xlabel("Amino acid")
    axs.set_ylabel("Percentage found in genome")
    axs.set_title(title)
    plt.savefig(path, bbox_inches="tight")

plot_dist(
    gene_dist,
    "Original Genome",
    "../figures/aa_dist_original.pdf")

gene_transnp = np.dot(gene_dist, trans_tablenp)
print(gene_transnp)

plot_dist(
    gene_transnp,
    "1 Single-point Mutation",
    "../figures/aa_dist_1p.pdf")

# Generate synonymous vs. non-synonymous
expected_strict_gene_transnp = gene_transnp*np.sum(gene_count)

gene_trans2p = np.dot(gene_dist, trans_table2p)
gene_trans3p = np.dot(gene_dist, trans_table3p)

synnonsyn_count_table = np.zeros((len(idx2aa), 2))

base_muts = BasicNeutralNetwork(1)
base_muts._build()
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
synonsyn_gene_trans_table2p = np.dot(gene_trans2p, synonsyn_trans_tablenp)
synonsyn_gene_trans_table3p = np.dot(gene_trans3p, synonsyn_trans_tablenp)

def plot_pie(synonsyn, title, path):
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.pie(synonsyn, labels=["Synonymous", "Non-synonymous"], autopct='%1.4f%%')
    ax.legend()
    ax.set_title(title)
    plt.savefig(path, bbox_inches="tight")

plot_pie(synonsyn_gene_trans_table1p, "After 1 Single-point Mutation",
    "../figures/synonsyn_1p.pdf")
plot_pie(synonsyn_gene_trans_table2p, "After 2 Single-point Mutations",
    "../figures/synonsyn_2p.pdf")
plot_pie(synonsyn_gene_trans_table3p, "After 3 Single-point Mutations",
    "../figures/synonsyn_3p.pdf")

def write_synnonsyn_table(path, table, dtype, caption, label):
    with open(path, "w") as file:
        file.write("\\begin{table}[]\n\t\centering\n\t\\begin{tabular}")
        file.write("{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c}\n")

        # Write the headers
        file.write("\t\t& Synonymous & Non-synonymous \\\\\n")

        for i, row in enumerate(table):
            file.write("\t\t\\hline\n")
            if idx2aa[i] == "_":
                file.write("\t\t\\_")
            else:
                file.write(f"\t\t{idx2aa[i]}")
            for val in row:
                if dtype == int:
                    file.write(f" & {int(val):2d}")
                elif dtype == float:
                    file.write(f" & {val:.2f}")
            file.write(" \\\\\n")

        file.write("\t\\end{tabular}\n")
        file.write("\t\\caption{" + caption + "}\n")
        file.write("\t\\label{tab:" + label + "}\n")
        file.write("\\end{table}\n")
write_synnonsyn_table(f"{TABLE_DIR}/synnonsyn_table.tex",
        synonsyn_trans_tablenp,
        float,
        caption="The ratio of synonymous versus non-synonymous mutations for each amino acid.",
        label="synnonsyn")

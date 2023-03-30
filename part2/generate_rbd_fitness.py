import pandas as pd
import re

START_SITE = 331
END_SITE   = 531
NUM_SITES = END_SITE - START_SITE + 1
AMINO_ACIDS = "RKHDEQNSTYWFAILMVGPC*"

def main(inpath, outpath, gene):
    df = pd.read_csv(inpath)
    df["keep"] = df.apply(lambda row:
            row["aa_site"] >= START_SITE and
            row["aa_site"] <= END_SITE and
            row["gene"] == gene,
        axis=1)
    df = df[df["keep"]]
    df["original"] = df["clade_founder_aa"]
    df["site"] = df["aa_site"]
    df = df.pivot(index=["site", "original"], columns="mutant_aa", values="delta_fitness")
    df.to_csv(outpath)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
            "Saves all of the fitness data relevant to RBD from the given fitness file")
    parser.add_argument("--file", type=str, default="part2/aamut_fitness_all.csv",
            help="The relative path to the amino acid mutation-fitness data")
    parser.add_argument("--out", type=str, default="part2/aamut_fitness_rbd.csv",
            help="The relative path to save the converted table to")
    parser.add_argument("--gene", type=str, default="S",
            help="Which gene to look at (default S)")
    args = parser.parse_args()
    main(args.file, args.out, args.gene)

import pandas as pd

START_SITE = 331
END_SITE   = 531
NUM_SITES = END_SITE - START_SITE + 1
AMINO_ACIDS = "RKHDEQNSTYWFAILMVGPC*"

def main(path):
    df = pd.read_csv(path)

    data = {
        "site": list(range(START_SITE, END_SITE + 1)),
        "aa": [""]*NUM_SITES
    }
    for aa in AMINO_ACIDS:
        data[f"{aa} mutation fitness"] = [0.0]*NUM_SITES

    MUT_I = 2
    FIT_I = 3
    for mutation in df["mutation"]: #(_, _, mutation, fitness) in df:
        print(mutation)
        site = int(mutation[1:-2])
        aa = mutation[0]

        if site < START_SITE or site > END_SITE:
            continue
        data["aa"][site] = aa

    df2 = pd.DataFrame(data)
    print(df2)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="part2/Neher_aa_fitness.csv",
            help="The relative path to the Neher amino acid fitness data")
    args = parser.parse_args()
    main(args.file)

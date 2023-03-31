import pandas as pd
import numpy as np


def main(inpath, n, filter):
    df = pd.read_csv(inpath)
    df['synonymous'] = df['escaped'] <= filter
    
    print(sum(df['count']))
    print(sum(df[df['synonymous'] == True]['count']))
    print(sum(df[df['synonymous'] == False]['count']))
    
    subset = df[df['synonymous'] == True]['count']
    uniq = [1, sum(subset)]

    for i in range(2, 1 + n):
        cur = uniq[i-1] * uniq[1] 
        for k in range(i-1):
            cur -= (-1)**k * np.sum(np.power(subset, 2+k)) * uniq[i-2-k]
        cur //= i
        uniq.append(cur)
    
    total = sum(df['count'])**n
    syn = uniq[-1]
    
    print('Synonymous:', syn)
    print('Non-synonymous', total - syn)
    print('Total', total)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
            "Get the approximate number of synonymous (antigenically neutral), non-synonymous, and the total number of mutations")
    parser.add_argument("--file", type=str, default="part1/escape_calculator.csv",
            help="The relative path to the antigen escape table")
    parser.add_argument("--n", type=int, default=1,
            help="Number of neutral network")
    parser.add_argument("--filter", type=float, default=0.15,
            help="Bloom filter threshold")
    args = parser.parse_args()
    main(args.file, args.n, args.filter)
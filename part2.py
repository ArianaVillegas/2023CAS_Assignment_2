

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--breadth", type=int, default=5,
            help="The breadth with which to search the mutation space")
    parser.add_argument("--max-dist", type=int, default=1,
            help="The maximum number of mutations to explore past the neutral network")

    args = parser.parse_args()

    print(f"Searching the mutation graph with a breadth of {args.breadth} and max distance {args.max_dist}")

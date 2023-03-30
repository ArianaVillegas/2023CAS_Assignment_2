import random

# Constants
START_SITE = 331
END_SITE   = 531
NEUTRAL_THRESHOLD = 0.5 # Any fitness change less than this is marked neutral
AMINO_ACIDS = "RKHDEQNSTYWFAILMVGPC*"

def is_neutral(mutations, fitness):
    if abs(fitness) <= NEUTRAL_THRESHOLD:
        return True
    return False

def get_fitness(mutations):
    return random.gauss(0.0, 1.0) * 5.0 / 3.0

class RBD:
    def __init__(self, mutations=None):
        if mutations is None:
            self.mutations = []
            self.is_neutral = True
            self.fitness = 0.0
        else:
            self.mutations = [mut for mut in mutations]
            self.fitness = get_fitness(mutations)
            self.is_neutral = is_neutral(mutations, self.fitness)

    def __repr__(self):
        res = "\n\tRBD ("
        res += f"fitness: {self.fitness: 4.2f}, "
        res += f"is_neutral: {str(self.is_neutral):5s}, "
        res += "mutations:"
        for site, aa in self.mutations:
            res += f" ({site:3d}, {aa})"
        res += ")"
        return res

    def get_mutated(self, site=None, amino_acid=None):
        if site is None:
            site = random.randint(START_SITE, END_SITE)
        if amino_acid is None:
            amino_acid = random.choice(AMINO_ACIDS)
        res = RBD(self.mutations)
        res.mutations.append((site, amino_acid))
        #TODO: sort mutations by site and remove mutations with same a.a.
        #TODO: sample mutations from empirical mutation distribution
        return res

def dfs_helper(node, explored, dist, depth, breadth, max_dist, max_depth):
    res = []
    if dist >= max_dist or depth >= max_depth:
        return res

    for i in range(breadth):
        child = node.get_mutated()
        while child in explored or child in res: #TODO: Recalculate until new
            child = node.get_mutated()
        res.append(child)
        tmp_dist = dist
        if child.is_neutral:
            tmp_dist = 0
        else:
            tmp_dist += 1
        tmp = dfs_helper(child, explored + res, tmp_dist, depth+1,
                breadth, max_dist, max_depth)
        res.extend(tmp)
    return res

def dfs(breadth, max_dist, max_depth):
    """ Depth-first search for the RBD mutation network

    Performs a depth-first search of the mutation network. The maximum depth
    for the search is achieved when a node is ``max_dist`` mutations outside
    of the neutral network, or we have explored ``max_depth`` mutations away
    from the root.

    Args:
      breadth: integer. The number of mutations samples from each node in the
        network
      max_dist: integer. The maximum number of mutations outside the neutral
        network to try. If a genome is ``max_dist`` mutations from the neutral
        network, it becomes a leaf node of the dfs tree.
      max_depth: integer. The maximum number of hops to explore from the root.

    Returns:
      explored: list(RBD)
        A list of the mutated RBD genomes explored during the dfs.
    """
    root = RBD()
    explored = [root]
    dist = 0
    depth = 0
    return explored + dfs_helper(root, explored, dist, depth,
            breadth, max_dist, max_depth)

def main(breadth, max_dist, max_depth):

    print(f"Searching the mutation graph with a breadth of {breadth} and max distance {max_dist}")
    res = dfs(breadth, max_dist, max_depth)
    print(f"Result: {res}")

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--breadth", type=int, default=5,
            help="The breadth with which to search the mutation space")
    parser.add_argument("--max-dist", type=int, default=1,
            help="The maximum number of mutations to explore past the neutral network")
    parser.add_argument("--max-depth", type=int, default=5,
            help="The maximum number of mutations to explore from the root")
    parser.add_argument("--seed", type=int, default=42,
            help="The seed used for the pseudo-random number generator")
    parser.add_argument("--neutral-threshold", type=int, default=0.5,
            help="The threshold for what fitness change is considered neutral")
    args = parser.parse_args()

    random.seed(args.seed)

    main(breadth=args.breadth, max_dist=args.max_dist, max_depth=args.max_depth)

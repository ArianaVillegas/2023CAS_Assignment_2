from graph_tool.all import Graph, graph_draw
from graph_tool.draw import radial_tree_layout
import random
import pandas as pd
import matplotlib
import math
import os

from pathlib import Path
import sys
sys.path.append("part2/SARS2_RBD_Ab_escape_maps")
import bindingcalculator as bc

# Global variables
START_SITE = 331
END_SITE   = 531
NEUTRAL_CRITERION = "fitness" # The criterion to determine neutrality
FITNESS_THRESHOLD = 0.5 # Any fitness change less than this is marked neutral
ESCAPE_THRESHOLD = 0.96 # Any % antibody binding above this is marked neutral
AMINO_ACIDS = "RKHDEQNSTYWFAILMVGPC*"
MUT_DATA_PATH = "part2/aamut_fitness_rbd.csv"
OUT_DIR = "figures/"
COLOR_GOOD = (121 / 255.0, 159 / 255.0, 203 / 255.0, 1.0)
COLOR_BAD =  (249 / 255.0, 102 / 255.0,  94 / 255.0, 1.0)
COLOR_ROOT = (0.0, 0.0, 0.0, 1.0)

g_mut_data = None
g_v_fitness_color = None
g_v_fitness_str = None
g_bc = None

def is_neutral(mutations, fitness, escape):
    if NEUTRAL_CRITERION == "fitness":
        return abs(fitness) <= FITNESS_THRESHOLD # delta fitness <= threshold
    if NEUTRAL_CRITERION == "escape":
        return escape >= ESCAPE_THRESHOLD # % of antibodies >= threshold
    return False

def is_valid_site(site):
    return (site in g_mut_data.index and site in g_bc.sites)

def is_valid_mutation(site, amino_acid):
    if (g_mut_data.loc[[site]]["original"] == amino_acid).any():
        return False
    val = g_mut_data.loc[[site]][amino_acid]
    if len(val) > 1:
        val = max(val)
    if math.isnan(val):
        return False
    return True

def get_fitness(mutations):
    global g_mut_data

    fitness = 0.0
    if g_mut_data is None:
        g_mut_data = pd.read_csv(MUT_DATA_PATH, index_col="site")

    for site, aa in mutations:
        tmp = g_mut_data.loc[[site]][aa].values[0]
        fitness += tmp

    return fitness

def get_escape(mutations):
    sites = [site for site, _ in mutations if site in g_bc.sites]
    res = g_bc.binding_retained(sites)
    return res

class RBD:
    def __init__(self, graph, mutations=None):
        self.graph = graph
        self.vertex = graph.add_vertex()
        if mutations is None:
            self.mutations = []
            self.is_neutral = True
            self.fitness = 0.0
            self.escape = 1.0
        else:
            self.mutations = [mut for mut in mutations]
            self.fitness = get_fitness(mutations)
            self.escape = get_escape(mutations)
            self.is_neutral = is_neutral(mutations, self.fitness, self.escape)
        if self.is_neutral:
            g_v_fitness_color[self.vertex] = (1, 1, 1, 1)
        elif self.fitness < 0:
            g_v_fitness_color[self.vertex] = COLOR_BAD
        elif self.fitness > 0:
            g_v_fitness_color[self.vertex] = COLOR_GOOD
        g_v_fitness_str[self.vertex] = f"{self.fitness:4.2f}"

    def __repr__(self):
        res = "\n\tRBD ("
        res += f"fitness: {self.fitness: 4.2f}, "
        res += f"escape: {self.escape: .2f}, "
        res += f"is_neutral: {str(self.is_neutral):5s}, "
        res += "mutations:"
        for site, aa in self.mutations:
            res += f" ({site:3d}, {aa})"
        res += ")"
        return res

    def get_mutated(self, site=None, amino_acid=None):
        while site is None or not is_valid_site(site):
            site = random.randint(START_SITE, END_SITE)
        while amino_acid is None or not is_valid_mutation(site, amino_acid):
            amino_acid = random.choice(AMINO_ACIDS)
        res = RBD(self.graph, self.mutations + [(site, amino_acid)])
        self.graph.add_edge(self.vertex, res.vertex)
#         res.mutations.append((site, amino_acid))
#         res.fitness = get_fitness(res.mutations)
#         res.is_neutral = is_neutral(res.mutations, res.fitness)
        #TODO: sort mutations by site and remove mutations with same a.a.
        #TODO: sample mutations from empirical mutation distribution
        return res

def dfs_helper(node, explored, dist, depth, breadth, max_dist, max_depth):
    res = []
    if dist >= max_dist or depth >= max_depth:
        return res

    for i in range(breadth):
        child = node.get_mutated()
        if child in explored or child in res:
            # TODO given the size of the search space, the probability of this
            # occuring is so low, I am not handling this case
            print("WHOAH")
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
    global g_v_fitness_color
    global g_v_fitness_str

    graph = Graph(directed=False)
    g_v_fitness_color = graph.new_vertex_property("vector<double>")
    g_v_fitness_str = graph.new_vertex_property("string")

    root = RBD(graph)
    g_v_fitness_color[root.vertex] = COLOR_ROOT
    explored = [root]
    dist = 0
    depth = 0

    # Generate search tree
    res = explored + dfs_helper(root, explored, dist, depth,
            breadth, max_dist, max_depth)

    graph_draw(
            graph,
#             pos=radial_tree_layout(graph, 0),
            vertex_fill_color=g_v_fitness_color,
            vcmap=matplotlib.cm.bwr,
#             vertex_text=g_v_fitness_str, #graph.vertex_index,
#             font_size=8,
            output=f"{OUT_DIR}/neutral_network_{NEUTRAL_CRITERION}_b{breadth}_d{max_depth}.pdf")
    return res

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
    parser.add_argument("--neutral-criterion", type=str, default="fitness",
            help="The criterion we use to determine neutrality of mutations. Can be either 'fitness' or 'escape'.")
    parser.add_argument("--fitness-threshold", type=float, default=1.0,
            help="The threshold for what fitness change is considered neutral for the fitness criterion.")
    parser.add_argument("--escape-threshold", type=float, default=0.96,
            help="The threshold for fraction of antibody binding is considered neutral for the escape criterion.")
    parser.add_argument("--mut-csv", type=str, default="part2/aamut_fitness_rbd.csv",
            help="The relative path to the rbd-mutation-fitness data")
    parser.add_argument("--out-dir", type=str, default="figures/",
            help="The directory to save figures into")
    args = parser.parse_args()

    random.seed(args.seed)
    FITNESS_THRESHOLD = args.fitness_threshold
    ESCAPE_THRESHOLD = args.escape_threshold
    NEUTRAL_CRITERION = args.neutral_criterion
    MUT_DATA_PATH = args.mut_csv
    OUT_DIR = args.out_dir

    g_mut_data = pd.read_csv(MUT_DATA_PATH, index_col="site")
    g_bc = bc.BindingCalculator()

    main(breadth=args.breadth, max_dist=args.max_dist, max_depth=args.max_depth)

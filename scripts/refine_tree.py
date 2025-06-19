import argparse
from Bio import Phylo


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Midpoint-root and ladderize a Newick tree"
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree file")
    parser.add_argument(
        "--out_tree", required=True, help="Output refined tree filename"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Read the input tree
    tree = Phylo.read(args.tree, "newick")

    # Midpoint root and ladderize
    tree.root_at_midpoint()
    tree.ladderize()

    # Write the refined tree
    Phylo.write(tree, args.out_tree, "newick")

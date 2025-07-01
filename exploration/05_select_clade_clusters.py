# %%
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo

from utils_04 import parse_mash_triangle_output, get_xy_positions


def get_top_MLST(df, top_n=30):
    """Get the top N MLST values from the dataframe."""
    return df["MLST"].value_counts().head(top_n).index.tolist()


def draw_MLST(tree, df, top_mlst, ax, color_label):
    positions = get_xy_positions(tree)
    cmap = plt.get_cmap("tab20")
    mlst_colors = {
        mlst: cmap(i % 20) for i, mlst in enumerate(top_mlst)
    }  # Assign a color to each MLST
    for clade in tree.get_terminals():
        mlst = df.loc[clade.name, color_label]
        if mlst in top_mlst:
            x, y = positions[clade.name]
            ax.plot(x, y, ".", color=mlst_colors[mlst])
    # add legend
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker=".",
            color="w",
            label=str(mlst),
            markerfacecolor=mlst_colors[mlst],
            markersize=10,
        )
        for mlst in top_mlst
    ]
    ax.legend(
        handles=handles,
        title=color_label,
        loc="lower left",
        frameon=False,
    )


# %%


def load_data(mash_file, tree_file, metadata_file):
    """Load and process mash distance data, phylogenetic tree, and metadata.

    Args:
        mash_file: Path to mash triangle output file
        tree_file: Path to phylogenetic tree file in newick format
        metadata_file: Path to combined metadata CSV file

    Returns:
        tuple: (tree, df, info) where tree is the phylogenetic tree,
               df is the reordered mash distance matrix, and info is the metadata
    """
    # Load mash distance data
    df = parse_mash_triangle_output(mash_file)

    # Load phylogenetic tree
    tree = Phylo.read(tree_file, "newick")
    order = [leaf.name for leaf in tree.get_terminals()]

    # Reorder DataFrame according to tree order
    df = df.reindex(index=order, columns=order)

    # Load metadata
    info = pd.read_csv(
        metadata_file,
        index_col=["chromosome_acc"],
        parse_dates=["Assembly Release Date"],
        dtype={"MLST": str},
    )

    return tree, df, info


species = "hpylori"

mash_fname = f"../results/{species}/mash_triangle.tsv"
tree_fname = f"../results/{species}/attotree.nwk"
info_fname = f"../results/{species}/combined_metadata.csv"
tree, df, info = load_data(mash_fname, tree_fname, info_fname)

# distance distribution plot
plt.figure(figsize=(10, 6))
sns.ecdfplot(df.values.flatten(), color="k")
# overlay for MLST
top_mlst = get_top_MLST(info, top_n=20)
for mlst in top_mlst:
    mask = info["MLST"] == mlst
    sns.ecdfplot(
        df.loc[info[mask].index, info[mask].index].values.flatten(),
        label=mlst,
    )
plt.xscale("log")
plt.title("Mash Distance Distribution")
plt.xlabel("Mash Distance")
plt.ylabel("Density")
plt.legend(title="MLST", loc="upper left")
plt.tight_layout()
plt.show()

# %%

# for every large ST (>10 isolates)
threshold = 10
mlst_counts = info["MLST"].value_counts()
large_mlst = mlst_counts[mlst_counts > threshold].index.tolist()
# take the common ancestor clade for each large MLST
mlst_clades = {}
for mlst in large_mlst:
    # get all isolates for this MLST
    MLST_isolates = info[info["MLST"] == mlst].index.tolist()
    # get the common ancestor clade
    clade = tree.common_ancestor(MLST_isolates)
    mlst_clades[mlst] = clade
print(f"Found {len(mlst_clades)} MLST clades with more than {threshold} isolates.")
# %%


# define a function to color all nodes downstream of a node
def color_downstream_nodes(node, color):
    """Color all nodes downstream of a given node."""
    node.color = color
    for child in node.clades:
        color_downstream_nodes(child, color)


def score(n_MLST, n_extra):
    """Calculate a score based on the number of MLST isolates and extra isolates."""
    alpha = 2  # weight for MLST isolates
    score = n_MLST * alpha - n_extra
    return score


# for each clade, assign to each internal clade node:
# - what fraction of all MLST isolates are descendants of this node
# - what fraction of all descendats of this node are MLST isolates
for mlst, clade in mlst_clades.items():
    # get all isolates for this MLST
    MLST_isolates = info[info["MLST"] == mlst].index.tolist()
    # assign to each internal node the fraction of MLST isolates
    nodes_df = []
    for node in clade.get_nonterminals():
        # get all isolates in the descendants
        descendants = [leaf.name for leaf in node.get_terminals()]
        # fraction of MLST isolates in the descendants
        n_leaves = len(descendants)
        n_MLST = len(set(descendants) & set(MLST_isolates))
        fraction_MLST_leaves = n_MLST / n_leaves
        fraction_MLST_in_descendants = n_MLST / len(MLST_isolates)
        # distance from the clade root
        dist = clade.distance(node)
        nodes_df.append(
            {
                "node": node,
                "fraction_MLST_leaves": fraction_MLST_leaves,
                "fraction_MLST_in_descendants": fraction_MLST_in_descendants,
                "clade_root_dist": dist,
                "n_MLST_included": n_MLST,
                "n_non_MLST_included": n_leaves - n_MLST,
                "n_leaves": n_leaves,
                "score": score(n_MLST, n_leaves - n_MLST),
            }
        )
    # convert to DataFrame
    nodes_df = pd.DataFrame(nodes_df)

    # select isolate with the highest score
    best_node = nodes_df.loc[nodes_df["score"].idxmax(), "node"]

    sns.scatterplot(
        data=nodes_df,
        x="fraction_MLST_leaves",
        y="fraction_MLST_in_descendants",
        hue="clade_root_dist",
        palette="rainbow",
        legend="brief",
    )
    # mark the best node
    plt.scatter(
        nodes_df.loc[nodes_df["node"] == best_node, "fraction_MLST_leaves"],
        nodes_df.loc[nodes_df["node"] == best_node, "fraction_MLST_in_descendants"],
        color="red",
        marker="x",
    )
    plt.xlabel("Fraction of captured isolates that have correct ST")
    plt.ylabel("Fraction of ST isolates that are captured")
    plt.title(f"MLST {mlst} Clade Analysis")
    plt.tight_layout()
    plt.show()

    # color all nodes downstream blue
    color_downstream_nodes(clade, "silver")
    color_downstream_nodes(best_node, "black")

    # draw clade
    fig, ax = plt.subplots(figsize=(10, 10))
    Phylo.draw(
        clade,
        do_show=False,
        label_func=lambda x: "",
        axes=ax,
    )
    # draw MLST isolates
    top_mlst = get_top_MLST(info, top_n=20)
    draw_MLST(clade, info, top_mlst, ax, "MLST")
    ax.set_title(f"Phylogenetic Tree for MLST {mlst}")
    plt.tight_layout()
    plt.show()


# %%

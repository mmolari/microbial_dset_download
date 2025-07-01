from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns
import argparse
from collections import defaultdict
import json
from utils import parse_mash_triangle_output, get_xy_positions, color_downstream_nodes


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Select clade clusters based on MLST analysis"
    )
    parser.add_argument(
        "--mash-file", required=True, help="Path to mash triangle output file"
    )
    parser.add_argument(
        "--tree-file",
        required=True,
        help="Path to phylogenetic tree file in newick format",
    )
    parser.add_argument(
        "--metadata-file", required=True, help="Path to combined metadata CSV file"
    )
    parser.add_argument(
        "--output-figs", required=True, help="Directory for output figures"
    )
    parser.add_argument(
        "--output-clusters", required=True, help="Path for output clusters file"
    )
    parser.add_argument("--species", required=True, help="Species name")
    parser.add_argument(
        "--thr-size", required=True, type=int, help="Threshold size parameter"
    )
    parser.add_argument(
        "--false-positive-penalty",
        required=True,
        type=float,
        help="Penalty for false positives (non-selected ST) in cluster",
    )
    return parser.parse_args()


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


def score(n_MLST, n_extra, false_positive_penalty):
    """Calculate a score based on the number of MLST isolates and extra isolates."""
    return n_MLST - false_positive_penalty * n_extra


def score_internal_nodes(root, all_st_isolates, false_positive_penalty):
    """Refine the root clade based on the number of MLST isolates and extra isolates."""
    n_all_st_isolates = len(all_st_isolates)
    nodes_df = []
    for node in root.get_nonterminals():
        all_descendants = [leaf.name for leaf in node.get_terminals()]
        n_leaves = len(all_descendants)
        n_st = len(set(all_st_isolates) & set(all_descendants))
        n_extra = n_leaves - n_st
        dist = root.distance(node)
        nodes_df.append(
            {
                "node": node,
                "n_leaves": n_leaves,
                "n_st": n_st,
                "n_extra": n_extra,
                "score": score(n_st, n_extra, false_positive_penalty),
                "recall": n_st / n_all_st_isolates,
                "precision": n_st / n_leaves,
                "distance": dist,
            }
        )
    nodes_df = pd.DataFrame(nodes_df)
    # Sort by score, then by number of leaves, then by distance
    nodes_df = nodes_df.sort_values(
        by=["score", "n_leaves", "distance"], ascending=[False, False, True]
    ).reset_index(drop=True)

    return nodes_df


def precision_recall_plot(nodes_df, ST, species, fname):
    """Create a precision-recall plot for the internal nodes."""
    sns.scatterplot(
        data=nodes_df,
        x="precision",
        y="recall",
        hue="distance",
        palette="rainbow",
        legend="brief",
    )
    # mark the best node
    best_node = nodes_df.iloc[0]
    plt.scatter(
        best_node["precision"],
        best_node["recall"],
        color="red",
        marker="x",
    )
    plt.xlabel("Precision (ST isolates selected / all selected)")
    plt.ylabel("Recall (ST isolates selected / all ST)")
    plt.title(f"precision recall plot - {species} - ST {ST}")
    plt.legend(title="Clade Root Distance")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()


def mark_MLST(root, info_df, ax):
    iso_names = [leaf.name for leaf in root.get_terminals()]
    df = info_df.loc[iso_names]

    # assign colors to MLSTs
    clade_mlst = df["MLST"].value_counts().index.tolist()
    mlst_colors = defaultdict(lambda: "black")  # default color for others
    for i, mlst in enumerate(clade_mlst[:20]):
        mlst_colors[mlst] = plt.get_cmap("tab20")(i % 20)

    positions = get_xy_positions(root)
    for clade in root.get_terminals():
        mlst = df.loc[clade.name, "MLST"]
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
        for mlst in clade_mlst[:20]
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        frameon=False,
    )


def plot_selected_clade_in_tree(root, best_node, info, species, ST, svfile):
    # color all nodes downstream blue
    color_downstream_nodes(root, "silver")
    color_downstream_nodes(best_node, "black")

    # draw clade
    fig, ax = plt.subplots(figsize=(8, 10))
    Phylo.draw(
        root,
        do_show=False,
        label_func=lambda x: "",
        axes=ax,
    )
    # draw MLST isolates
    mark_MLST(root, info, ax)
    ax.set_title(f"{species} - ST {ST} clade selection")
    plt.tight_layout()
    plt.savefig(svfile, dpi=300)
    plt.close(fig)


def display_results(nodes_df, clade_root, info_df, svfld, ST, species):
    # display ROC curve
    precision_recall_plot(nodes_df, ST, species, svfld / f"ST_{ST}_ROC.png")
    # and sub-tree
    best_node = nodes_df.iloc[0]["node"]
    svfile = svfld / f"ST_{ST}_clade_selection.png"
    plot_selected_clade_in_tree(clade_root, best_node, info_df, species, ST, svfile)


def save_results(clade_root, info, mash_dist, svfld, ST, species):
    # save list of isolate accession numbers
    fname = svfld / f"ST_{ST}_isolates.txt"
    iso_names = [leaf.name for leaf in clade_root.get_terminals()]
    np.savetxt(fname, np.array(iso_names), fmt="%s")

    # save metadata for the selected clade
    fname = svfld / f"ST_{ST}_isolates_metadata.csv"
    df = info.loc[iso_names]
    df.to_csv(fname, index_label="chromosome_acc")

    # save general information
    fname = svfld / f"ST_{ST}_info.json"
    all_MLST_counts = info["MLST"].value_counts()
    clade_MLST_counts = df["MLST"].value_counts()
    md = mash_dist.loc[iso_names, iso_names]
    # mean of upper-triangular matrix
    mean_mash_dist = md.values[np.triu_indices(len(md), k=1)].mean()
    std_mash_dist = md.values[np.triu_indices(len(md), k=1)].std()

    info_dict = {
        "species": species,
        "ST": ST,
        "n_isolates": len(iso_names),
        "n_ST_isolates": int(clade_MLST_counts[ST]),
        "all_ST": [
            (int(ct), int(all_MLST_counts[st])) for st, ct in clade_MLST_counts.items()
        ],
        "mean_mash_dist": mean_mash_dist,
        "std_mash_dist": std_mash_dist,
    }
    # save as JSON
    with open(fname, "w") as f:
        json.dump(info_dict, f, indent=4)


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()

    # Load data
    tree, df, info = load_data(args.mash_file, args.tree_file, args.metadata_file)

    # create output directories if they do not exist
    figs_fld = Path(args.output_figs)
    figs_fld.mkdir(parents=True, exist_ok=True)
    clusters_fld = Path(args.output_clusters)
    clusters_fld.mkdir(parents=True, exist_ok=True)

    for ST, count in info["MLST"].value_counts().items():
        # only consider STs with enough isolates
        if count < args.thr_size:
            continue

        print(f"{args.species} - Processing ST {ST} with {count} isolates")

        # get list of isolates for this ST
        st_isolates = info[info["MLST"] == ST].index.tolist()
        # get clade root for this ST
        st_root = tree.common_ancestor(st_isolates)

        # score internal nodes
        nodes_df = score_internal_nodes(
            st_root, st_isolates, args.false_positive_penalty
        )

        # pick best node based on score
        best_node = nodes_df.iloc[0].to_dict()
        # if less than threshold size, skip
        if best_node["n_leaves"] < args.thr_size:
            print(
                f"Skipping ST {ST} after scoring: only {best_node['n_leaves']} leaves found"
            )
            continue

        # create result figures
        display_results(nodes_df, st_root, info, figs_fld, ST, args.species)

        # save results
        save_results(st_root, info, df, clusters_fld, ST, args.species)

# %%
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo

from sklearn.cluster import DBSCAN


# %%


def parse_mash_triangle_output(file_path: str) -> pd.DataFrame:
    """
    Parse a lower-triangular matrix file produced by `mash triangle` into a pandas DataFrame.
    """
    with open(file_path, "r") as f:
        # First line: number of sequences
        n = int(f.readline().split()[-1])

        # Preallocate an n×n float array filled with zeros
        mat = np.zeros((n, n), dtype=float)
        ids = []

        for i, line in enumerate(f):
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            fp, *vals = parts
            ids.append(Path(fp).stem)

            if vals:
                # parse into 1d array of length i
                row_vals = np.fromiter(
                    (float(x) for x in vals), dtype=float, count=len(vals)
                )
                # fill lower triangle [i, 0:i] and mirror into [0:i, i]
                mat[i, : len(row_vals)] = row_vals
                mat[: len(row_vals), i] = row_vals

    # wrap in DataFrame
    return pd.DataFrame(mat, index=ids, columns=ids)


def get_xy_positions(tree):
    """Get x, y positions for each clade in the tree."""
    positions = {}
    depths = tree.depths()
    for i, clade in enumerate(tree.get_terminals()):
        positions[clade.name] = (depths[clade], i + 1)
    return positions


def get_top_MLST(df, top_n=30):
    """Get the top N MLST values from the dataframe."""
    df_mlst = df[["MLST"]].copy()
    return df_mlst["MLST"].value_counts().head(top_n).index.tolist()


def draw_MLST(tree, df, top_mlst, ax, positions):
    """Draw MLST labels on the tree."""
    cmap = plt.get_cmap("tab20")
    mlst_colors = {
        mlst: cmap(i % 20) for i, mlst in enumerate(top_mlst)
    }  # Assign a color to each MLST
    for clade in tree.get_terminals():
        mlst = df.loc[clade.name, "MLST"]
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
        title="MLST",
        loc="lower left",
        frameon=False,
    )


def plot_tree(tree, df, ax):
    top_MLST = get_top_MLST(df, top_n=20)
    positions = get_xy_positions(tree)

    Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
    draw_MLST(tree, df, top_MLST, ax, positions)


def cluster_distance_matrix(
    distance_matrix: pd.DataFrame, eps: float, min_samples: int
) -> np.ndarray:
    """
    Cluster a distance matrix using DBSCAN.

    Args:
        distance_matrix: Square DataFrame with pairwise distances
        eps: Maximum distance between samples in the same neighborhood
        min_samples: Minimum number of samples in a neighborhood for a core point

    Returns:
        Array of cluster labels (-1 for noise points)
    """
    db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    return db.fit_predict(distance_matrix.values)


# %%

species = "saureus"

fname = f"../results/{species}/mash_triangle.tsv"
df = parse_mash_triangle_output(fname)

# get order from the tree
fname = f"../results/{species}/attotree.nwk"
tree = Phylo.read(fname, "newick")
order = [leaf.name for leaf in tree.get_terminals()]

# reorder DataFrame
df = df.reindex(index=order, columns=order)

# get metadata
fname = f"../results/{species}/combined_metadata.csv"
info = pd.read_csv(
    fname,
    index_col=["chromosome_acc"],
    parse_dates=["Assembly Release Date"],
    dtype={"MLST": str},
)

# get the 95th percentile threshold for the heatmap
high_threshold = np.quantile(df.values, 0.95)


# Create subplots with tree on the left and heatmap on the right
fig, (ax1, ax2) = plt.subplots(
    1,
    2,
    figsize=(16, 8),
    gridspec_kw={"width_ratios": [1, 2]},
    sharey=True,
)

# Plot tree on the left
plot_tree(tree, info, ax1)

# Plot heatmap on the right
sns.heatmap(
    df,
    ax=ax2,
    cmap="viridis",
    square=True,
    cbar_kws={"label": "Mash distance"},
    vmin=0,
    vmax=high_threshold,
    xticklabels=False,
    yticklabels=False,
)
ax2.set_title("Mash Distances Heatmap")

plt.tight_layout()
plt.show()


# %%
def plot_mash_distance_analysis(df, info, top_n=20):
    """
    Plot histogram and cumulative distributions of mash distances.

    Args:
        df: Mash distance matrix (DataFrame)
        info: Metadata DataFrame with MLST information
        top_n: Number of top MLST types to include in cumulative plot
    """
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Left subplot: histogram of mash distances
    sns.histplot(
        df.values.flatten(),
        stat="density",
        ax=ax1,
    )
    ax1.set_title("Distribution of Mash Distances")
    ax1.set_xlabel("Mash Distance")
    ax1.set_ylabel("Density")

    # Right subplot: cumulative distributions of distances for top sequence types
    top_mlst = get_top_MLST(info, top_n=top_n)
    cmap = plt.get_cmap("tab20")

    for i, mlst in enumerate(top_mlst):
        # get samples with this MLST
        mlst_samples = info[info["MLST"] == mlst].index
        # filter to samples present in distance matrix
        mlst_samples = [s for s in mlst_samples if s in df.index]

        if len(mlst_samples) > 1:
            # get pairwise distances within this MLST group
            mlst_df = df.loc[mlst_samples, mlst_samples]
            # get upper triangle (excluding diagonal)
            mask = np.triu(np.ones_like(mlst_df, dtype=bool), k=1)
            distances = mlst_df.values[mask]

            sns.ecdfplot(
                distances,
                label=f"MLST {mlst}",
                color=cmap(i % 20),
                stat="proportion",
                ax=ax2,
            )

    # freeze x-limits
    xlim = ax2.get_xlim()

    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    sns.ecdfplot(
        df.values[mask],
        label="All Samples",
        color="black",
        stat="proportion",
        linestyle="--",
        ax=ax2,
    )

    ax2.set_xlim(xlim)
    ax2.set_xlabel("Mash Distance")
    ax2.set_ylabel("Cumulative Probability")
    ax2.set_title("Cumulative Distribution of Mash Distances by Sequence Type")
    ax2.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    plt.show()


# Call the function
plot_mash_distance_analysis(df, info)


# %%
# cluster labels
labels = cluster_distance_matrix(df, eps=0.002, min_samples=10)

# create a dataframe with cluster labels and MLST
cluster_df = pd.DataFrame(
    {
        "MLST": info["MLST"],
        "Cluster": labels,
    },
    index=info.index,
)

# Get MLST counts and assign -1 to those with less than 10 elements
mlst_counts = cluster_df["MLST"].value_counts()
small_mlst = mlst_counts[mlst_counts < 10].index
cluster_df.loc[cluster_df["MLST"].isin(small_mlst), "MLST"] = "rest"

# visualize confusion matrix
confusion_matrix = pd.crosstab(
    cluster_df["MLST"],
    cluster_df["Cluster"],
    # margins=True,
    # normalize="columns",
)
plt.figure(figsize=(12, 8))
sns.heatmap(
    confusion_matrix,
    annot=True,
    cmap="Blues",
    cbar_kws={"label": "Count"},
)
plt.title("Confusion Matrix of MLST vs Cluster Labels")
plt.xlabel("Cluster Labels")
plt.ylabel("MLST")
plt.tight_layout()
plt.show()

# %%

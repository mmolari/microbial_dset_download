# %%
import numpy as np
import pandas as pd
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo

from utils_04 import parse_mash_triangle_output, get_xy_positions

# %%


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

species = "hpylori"

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

# distance distribution plot
plt.figure(figsize=(10, 6))
sns.ecdfplot(df.values.flatten(), color="k")
# overlay for MLST
top_mlst = get_top_MLST(info, top_n=20)
for mlst in top_mlst:
    sns.ecdfplot(
        df.loc[
            info[info["MLST"] == mlst].index, info[info["MLST"] == mlst].index
        ].values.flatten(),
        label=mlst,
    )
plt.xscale("log")
plt.title("Mash Distance Distribution")
plt.xlabel("Mash Distance")
plt.ylabel("Density")
plt.legend(title="MLST")
plt.tight_layout()
plt.show()

# %%
from sklearn.cluster import HDBSCAN


def cluster_distance_matrix(
    distance_matrix: pd.DataFrame, eps: float, min_samples: int
) -> np.ndarray:
    """
    Cluster a distance matrix using HDBSCAN.
    """
    db = HDBSCAN(
        cluster_selection_epsilon=eps,
        min_cluster_size=min_samples,
        metric="precomputed",
    )
    return db.fit_predict(distance_matrix.values)


labels = cluster_distance_matrix(df, eps=0.005, min_samples=10)

# create a dataframe with cluster labels and MLST
cluster_df = pd.DataFrame(
    {
        "MLST": info.loc[df.index, "MLST"].values,
        "Cluster": labels,
    },
    index=df.index,
)

# Get MLST counts and assign -1 to those with less than 10 elements
mlst_counts = cluster_df["MLST"].value_counts()
small_mlst = mlst_counts[mlst_counts < 10].index
cluster_df.loc[cluster_df["MLST"].isin(small_mlst), "MLST"] = "small"
# set NAN to "rest"
cluster_df["MLST"].fillna("NA", inplace=True)

# visualize confusion matrix
confusion_matrix = pd.crosstab(
    cluster_df["MLST"],
    cluster_df["Cluster"],
    # margins=True,
    normalize="index",  # Normalize by row to get proportions
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

fig, axs = plt.subplots(1, 2, figsize=(12, 20))

ax = axs[0]
# plot tree
Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
top_MLST = get_top_MLST(info, top_n=20)
draw_MLST(tree, cluster_df, top_MLST, ax, "MLST")
ax.set_title("Phylogenetic Tree with MLST")

ax = axs[1]
# plot tree
Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
# plot cluster labels
top_clusters = cluster_df["Cluster"].value_counts().head(20).index.tolist()
draw_MLST(tree, cluster_df, top_clusters, ax, "Cluster")
ax.set_title("Phylogenetic Tree with Cluster Labels")
plt.tight_layout()
plt.show()

# %%

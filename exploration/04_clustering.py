# %%
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo
from sklearn.cluster import HDBSCAN, DBSCAN


from utils_04 import get_xy_positions, load_data
from sklearn.metrics import (
    mutual_info_score,
    normalized_mutual_info_score,
    adjusted_rand_score,
)

# %%


def get_top_MLST(df, top_n=30):
    """Get the top N MLST values from the dataframe."""
    return df["MLST"].value_counts().head(top_n).index.tolist()


def draw_MLST(tree, series, ax):
    positions = get_xy_positions(tree)
    mask = series == -1
    mask |= series.isna()
    s_count = series[~mask].value_counts()
    cmap = plt.get_cmap("tab20")
    colors = defaultdict(lambda: "black")
    for i, (s, ct) in enumerate(s_count[:20].items()):
        colors[s] = cmap(i % 20)

    for clade in tree.get_terminals():
        s = series.loc[clade.name]
        # skip nan and -1
        if pd.isna(s) or s == -1:
            continue
        x, y = positions[clade.name]
        ax.plot(x, y, ".", color=colors[s])
    # add legend
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker=".",
            color="w",
            label=str(s),
            markerfacecolor=colors[s],
            markersize=10,
        )
        for s in s_count.index[:20]
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        frameon=False,
    )


def cluster_distance_matrix_hdbscan(
    distance_matrix: pd.DataFrame, eps: float, min_samples: int
) -> pd.Series:
    """
    Cluster a distance matrix using HDBSCAN.
    """
    db = HDBSCAN(
        cluster_selection_epsilon=eps,
        min_cluster_size=min_samples,
        metric="precomputed",
    )
    labels = db.fit_predict(distance_matrix.values)
    # create a pandas Series with the labels
    labels_series = pd.Series(labels, index=distance_matrix.index)
    return labels_series


# def cluster_distance_matrix_dbscan(
#     distance_matrix: pd.DataFrame, eps: float, min_samples: int
# ) -> pd.Series:
#     """
#     Cluster a distance matrix using DBSCAN.
#     """
#     db = DBSCAN(
#         eps=eps,
#         min_samples=min_samples,
#         metric="precomputed",
#     )
#     labels = db.fit_predict(distance_matrix.values)
#     # create a pandas Series with the labels
#     labels_series = pd.Series(labels, index=distance_matrix.index)
#     return labels_series


def evaluate_clustering(eps: float) -> tuple[dict, pd.Series]:
    """
    Evaluate clustering performance for a given epsilon value.

    Args:
        eps: Epsilon parameter for HDBSCAN clustering

    Returns:
        tuple: (clustering_results_dict, labels_series)
    """
    labels_series = cluster_distance_matrix_hdbscan(mash_df, eps=eps, min_samples=10)
    cluster_df = pd.concat(
        {
            "MLST": mlst_conf_series,
            "Cluster": labels_series,
        },
        axis=1,
    )

    # Calculate clustering metrics
    mi = mutual_info_score(cluster_df["MLST"], cluster_df["Cluster"])
    nmi = normalized_mutual_info_score(cluster_df["MLST"], cluster_df["Cluster"])
    ari = adjusted_rand_score(cluster_df["MLST"], cluster_df["Cluster"])

    clustering_results = {
        "eps": eps,
        "mi": mi,
        "nmi": nmi,
        "ari": ari,
        "num_clusters": len(set(labels_series) - {-1}),
        "num_mlst": len(set(mlst_conf_series) - {"small", "NA"}),
        "num_samples": len(labels_series),
        "num_samples_in_clusters": len(labels_series[labels_series != -1]),
        "num_samples_in_mlst": len(
            mlst_conf_series[~mlst_conf_series.isin(["small", "NA"])]
        ),
    }

    return clustering_results, cluster_df


def plot_clustering_results(clustering_results: pd.DataFrame, best_eps: float):
    """
    Plot clustering results for different epsilon values.
    """
    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    metric_labels = {
        "mi": "Mutual Information",
        "nmi": "Normalized Mutual Information",
        "ari": "Adjusted Rand Index",
        "num_clusters": "Number of Clusters",
        "num_mlst": "Number of (large) MLST Types",
        "num_samples": "Total Samples",
        "num_samples_in_clusters": "Samples in Clusters",
        "num_samples_in_mlst": "Samples in MLST Types",
    }

    for metrics, ax in zip(
        [
            ["mi", "nmi", "ari"],
            ["num_clusters", "num_mlst"],
            ["num_samples", "num_samples_in_clusters", "num_samples_in_mlst"],
        ],
        axs,
    ):
        for metric in metrics:
            sns.lineplot(
                data=clustering_results,
                x="eps",
                y=metric,
                label=metric_labels[metric],
                marker="o",
                ax=ax,
            )
        ax.set_ylim(bottom=0)
        ax.axvline(best_eps, color="red", linestyle="--", label="Best Epsilon")
        ax.set_xscale("log")
        ax.legend(
            loc="upper left", fontsize="small", frameon=False, bbox_to_anchor=(1, 1)
        )
        ax.grid(True)

    plt.tight_layout()
    plt.show()


# %%

# species = "hpylori"
species = "paeruginosa"

mash_fname = f"../results/{species}/mash_triangle.tsv"
tree_fname = f"../results/{species}/attotree.nwk"
info_fname = f"../results/{species}/combined_metadata.csv"
tree, mash_df, info = load_data(mash_fname, tree_fname, info_fname)
best_eps_dict = {
    "bpertussis": 0.005,
    "ecoli": 0.003,
    "hpylori": 0.005,
    "kpneumoniae": 0.002,
    "mtuberculosis": 0.005,
    "paeruginosa": 0.004,
    "saureus": 0.005,
    "sflexneri": 0.005,
    "spneumoniae": 0.005,
    "spyogenes": 0.005,
}
best_eps = best_eps_dict[species]
# %%

# distance distribution plot
plt.figure(figsize=(10, 6))
sns.ecdfplot(mash_df.values.flatten(), color="k")
# overlay for MLST
top_mlst = get_top_MLST(info, top_n=20)
for mlst in top_mlst:
    sns.ecdfplot(
        mash_df.loc[
            info[info["MLST"] == mlst].index, info[info["MLST"] == mlst].index
        ].values.flatten(),
        label=mlst,
    )
plt.xscale("log")
plt.title("Mash Distance Distribution")
plt.xlabel("Mash Distance")
plt.ylabel("Density")
plt.legend(title="MLST", loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()

# %%

mlst_series = info["MLST"].copy()
mlst_counts = mlst_series.value_counts()
small_mlst = mlst_counts[mlst_counts < 10].index
mlst_conf_series = mlst_series.copy()
mlst_conf_series.loc[mlst_conf_series.isin(small_mlst)] = "small"
mlst_conf_series.fillna("NA", inplace=True)


clustering_results = []
for eps in np.logspace(-3.5, -1.5, 25):
    results, _ = evaluate_clustering(eps)
    clustering_results.append(results)
# Convert results to DataFrame
clustering_results_df = pd.DataFrame(clustering_results)
plot_clustering_results(clustering_results_df, best_eps=best_eps)

# %%


# pick best epsilon
cluster_res, cluster_df = evaluate_clustering(best_eps)

# visualize confusion matrix
confusion_matrix = pd.crosstab(
    cluster_df["MLST"],
    cluster_df["Cluster"],
    # margins=True,
    normalize="index",  # Normalize by row to get proportions
    # normalize="columns",  # Normalize by column to get proportions
)

ordered_mlst = cluster_df["MLST"].value_counts().index.tolist()
ordered_clusters = cluster_df["Cluster"].value_counts().index.tolist()
confusion_matrix = confusion_matrix.reindex(ordered_mlst, axis=0)
confusion_matrix = confusion_matrix.reindex(ordered_clusters, axis=1)

plt.figure(figsize=(12, 8))
sns.heatmap(
    confusion_matrix,
    # annot=True,
    cmap="Greys",
    cbar_kws={"label": "Proportion"},
    linewidths=0.5,
    linecolor="black",
)
plt.title(
    f"Confusion Matrix ({cluster_res['num_clusters']} Clusters, {cluster_res['num_mlst']} MLSTs, {cluster_res['num_samples']} Samples)"
)
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
draw_MLST(tree, mlst_series, ax)
ax.set_title("MLST")

ax = axs[1]
# plot tree
Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
# plot cluster labels
top_clusters = (
    cluster_df[cluster_df["Cluster"] != -1]["Cluster"]
    .value_counts()
    .head(20)
    .index.tolist()
)
draw_MLST(tree, cluster_df["Cluster"], ax)
ax.set_title("Cluster Labels")
plt.tight_layout()
plt.show()

# %%

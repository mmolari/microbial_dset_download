import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
import argparse
import json

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo
from sklearn.cluster import HDBSCAN


from utils_cluster import get_xy_positions, load_data
from sklearn.metrics import (
    mutual_info_score,
    normalized_mutual_info_score,
    adjusted_rand_score,
)


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
        "--eps",
        required=True,
        type=float,
        help="Epsilon parameter for HDBSCAN clustering",
    )
    parser.add_argument(
        "--assignment-threshold-freq",
        required=True,
        type=float,
        help="Threshold for ST assignment frequency in clusters",
    )
    return parser.parse_args()


def mash_distance_overview(
    mash_df: pd.DataFrame, info_df: pd.DataFrame, best_eps: float, svname: str
):
    # distance distribution plot
    plt.figure(figsize=(10, 6))
    sns.ecdfplot(mash_df.values.flatten(), color="k")
    # overlay for MLST
    mlst_count = info_df["MLST"].value_counts()
    top_mlst = mlst_count.head(20).index
    for mlst in top_mlst:
        sns.ecdfplot(
            mash_df.loc[
                info_df[info_df["MLST"] == mlst].index,
                info_df[info_df["MLST"] == mlst].index,
            ].values.flatten(),
            label=f"ST {mlst} (n={mlst_count[mlst]})",
        )
    plt.axvline(best_eps, color="red", linestyle="--", label="selected epsilon")
    plt.xscale("log")
    plt.title("Mash Distance Distribution")
    plt.xlabel("Mash Distance")
    plt.ylabel("Density")
    plt.legend(title="MLST", loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(svname, dpi=300)
    plt.close()


def refined_mlst_series(info_df: pd.DataFrame, thr_size: int) -> pd.Series:
    """Refine MLST series to compare with clustering"""

    mlst_series = info_df["MLST"].copy()
    mlst_counts = mlst_series.value_counts()
    small_mlst = mlst_counts[mlst_counts < thr_size].index
    mlst_series = mlst_series.where(~mlst_series.isin(small_mlst), "small")
    mlst_series = mlst_series.fillna("NA")

    return mlst_series


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


def evaluate_clustering(
    eps: float, mlst_series: pd.Series, mash_df: pd.DataFrame, min_size: int
) -> tuple:
    """
    Evaluate clustering performance for a given epsilon value.

    Args:
        eps: Epsilon parameter for HDBSCAN clustering

    Returns:
        tuple: (clustering_results_dict, labels_series)
    """
    labels_series = cluster_distance_matrix_hdbscan(
        mash_df, eps=eps, min_samples=min_size
    )
    cluster_df = pd.concat(
        {
            "MLST": mlst_series,
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
        "num_mlst": len(set(mlst_series) - {"small", "NA"}),
        "num_samples": len(labels_series),
        "num_samples_in_clusters": len(labels_series[labels_series != -1]),
        "num_samples_in_mlst": len(mlst_series[~mlst_series.isin(["small", "NA"])]),
        "avg_cluster_size": len(labels_series[labels_series != -1])
        / len(set(labels_series) - {-1}),
    }

    return clustering_results, cluster_df


def clustering_benchmark_df(
    eps_values, mlst_series, mash_df, min_size: int
) -> pd.DataFrame:
    """
    Create a DataFrame to benchmark clustering performance across multiple epsilon values.

    Args:
        eps_values: List of epsilon values to evaluate
        mlst_series: Series containing MLST information
        mash_df: DataFrame containing mash distances

    Returns:
        DataFrame with clustering results for each epsilon value
    """
    results = []
    for eps in eps_values:
        clustering_results, _ = evaluate_clustering(eps, mlst_series, mash_df, min_size)
        results.append(clustering_results)

    return pd.DataFrame(results)


def plot_clustering_benchmark(
    clustering_results: pd.DataFrame, best_eps: float, thr_size: int, svname: str
):
    """
    Plot clustering results for different epsilon values.
    """
    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    metric_labels = {
        "mi": "Mutual Information",
        "nmi": "Normalized Mutual Information",
        "ari": "Adjusted Rand Index",
        "num_clusters": "Number of Clusters",
        "num_mlst": f"Number of large (>{thr_size}) MLST Types",
        "num_samples": "Total isolates",
        "num_samples_in_clusters": "N. isolates in Clusters",
        "num_samples_in_mlst": "N. isolates in MLST Types",
        "avg_cluster_size": "Average Cluster Size",
    }

    for metrics, ax in zip(
        [
            ["mi", "nmi", "ari"],
            ["num_clusters", "num_mlst"],
            [
                "num_samples",
                "num_samples_in_clusters",
                "num_samples_in_mlst",
                "avg_cluster_size",
            ],
        ],
        axs,
    ):
        for metric in metrics:
            sns.lineplot(
                data=clustering_results,
                x="eps",
                y=metric,
                label=metric_labels[metric],
                marker=".",
                ax=ax,
            )
        ax.set_ylim(bottom=0)
        ax.axvline(best_eps, color="red", linestyle="--", label="Selected Epsilon")
        ax.set_xscale("log")
        ax.set_ylabel("")
        ax.legend(
            loc="upper left", fontsize="small", frameon=False, bbox_to_anchor=(1, 1)
        )
        ax.grid(True)

    plt.tight_layout()
    plt.savefig(svname, dpi=300)
    plt.close(fig)


def plot_confusion_matrix(
    cluster_df: pd.DataFrame, cluster_res: dict, eps: float, svname: str
):
    """
    Plot confusion matrix of clusters vs MLST types.
    """
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
        cbar_kws={"label": "proportion of MLST isolates"},
        linewidths=0.5,
        linecolor="black",
    )
    plt.title(
        f"Confusion Matrix ({cluster_res['num_clusters']} Clusters, {cluster_res['num_mlst']} (large) MLSTs, {cluster_res['num_samples']} Samples)"
    )
    plt.xlabel("Cluster Labels")
    plt.ylabel("MLST")
    plt.tight_layout()
    plt.savefig(svname, dpi=300)
    plt.close()


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


def plot_tree_clustering(
    tree: Phylo.BaseTree.Tree, cluster_df: pd.DataFrame, svname: str
):
    """
    Plot phylogenetic tree with MLST and cluster labels.
    """
    fig, axs = plt.subplots(1, 2, figsize=(12, 20))

    ax = axs[0]
    # plot tree
    Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
    draw_MLST(tree, cluster_df["MLST"], ax)
    ax.set_title("MLST")

    ax = axs[1]
    # plot tree
    Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
    # plot cluster labels
    draw_MLST(tree, cluster_df["Cluster"], ax)
    ax.set_title("Clusters")
    plt.tight_layout()
    plt.savefig(svname, dpi=300)
    plt.close(fig)


def summary_stats(
    info_df: pd.DataFrame,
    mash_df: pd.DataFrame,
    isolates: list,
    species: str,
    assignment_threshold: float,
) -> dict:
    """
    Generate summary statistics for a cluster.
    """

    # list of pairwise mash distances
    mash = mash_df.loc[isolates, isolates]
    mash = mash.values[np.triu_indices(len(isolates), k=1)]

    # list of STs in the cluster vs all STs in the species
    mlst_all = info_df["MLST"].value_counts()
    mlst_cluster = info_df.loc[isolates, "MLST"].value_counts()

    # dictionary of number of isolates with a given ST
    st_dict = {st: (int(ct), int(mlst_all[st])) for st, ct in mlst_cluster.items()}
    st_assignment = ""
    st_included = []
    for st, (ct, tot) in st_dict.items():
        if ct / tot > assignment_threshold:
            st_included.append(st)

    summary_stats = {
        "num_isolates": len(isolates),
        "species": species,
        "mash_dist_avg": mash.mean(),
        "mash_dist_std": mash.std(),
        "ST_n_different": len(mlst_cluster),  # number of different STs in the cluster
        "ST_included": st_included,  # list of STs included in the cluster
        "ST_all": st_dict,  # dictionary of STs with counts in the cluster and total
        "ST_n_isolates_with_type": int(
            mlst_cluster.sum()
        ),  # number of isolates with STs in the cluster
    }
    # if there is a ST_majority, add precision and recall
    if len(mlst_cluster) > 0:
        st_majority = mlst_cluster.idxmax()  # most common ST in the cluster
        # frequency of the most common ST in the cluster
        st_max_freq = mlst_cluster.max() / len(isolates)

        focal_st_counts = st_dict[st_majority]
        precision = focal_st_counts[0] / len(isolates)
        recall = focal_st_counts[0] / focal_st_counts[1]
        if precision > assignment_threshold and recall > assignment_threshold:
            st_assignment = st_majority
        else:
            st_assignment = None
        added_stats = {
            "ST_most_common": st_majority,  # most common ST in the cluster
            "ST_most_common_freq": st_max_freq,  # frequency of the most common ST in the cluster
            "ST_assignment": st_assignment,  # ST assignment based on frequency threshold
            "precision": precision,  # precision of the ST assignment
            "recall": recall,  # recall of the ST assignment
        }
    else:
        added_stats = {
            "ST_most_common": None,
            "ST_most_common_freq": 0.0,
            "precision": None,
            "recall": None,
            "ST_assignment": None,
        }
    summary_stats.update(added_stats)
    return summary_stats


def save_clustering_results(
    clust_df: pd.DataFrame,
    info_df: pd.DataFrame,
    mash_df: pd.DataFrame,
    clust_res: dict,
    species: str,
    assignment_threshold: float,
    fld: Path,
):
    """
    Save clustering results to a JSON file.
    """
    fname = fld / "clustering_results.json"
    with open(fname, "w") as f:
        json.dump(clust_res, f, indent=4)

    # for each cluster:
    # - save a txt file with the list of isolates
    # - save a json file with summary statistics
    # - save a csv file with metadata
    for n_clust in clust_df["Cluster"].unique():
        if n_clust == -1:
            continue
        isolates = clust_df[clust_df["Cluster"] == n_clust].index.tolist()
        # save isolates as txt file
        fname = fld / f"CL_{n_clust}_isolates.txt"
        np.savetxt(fname, isolates, fmt="%s")
        # save metadata as csv file
        fname = fld / f"CL_{n_clust}_metadata.csv"
        clust_info = info_df.loc[isolates]
        clust_info.to_csv(fname, index=True)
        # save summary statistics as json file
        stats = summary_stats(
            info_df=info_df,
            mash_df=mash_df,
            isolates=isolates,
            species=species,
            assignment_threshold=assignment_threshold,
        )
        fname = fld / f"CL_{n_clust}_info.json"
        with open(fname, "w") as f:
            json.dump(stats, f, indent=4)


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()

    # Load data
    print(f"Loading data for species: {args.species}")
    tree, mash_df, info = load_data(args.mash_file, args.tree_file, args.metadata_file)

    # create output directories if they do not exist
    figs_fld = Path(args.output_figs)
    figs_fld.mkdir(parents=True, exist_ok=True)
    clusters_fld = Path(args.output_clusters)
    clusters_fld.mkdir(parents=True, exist_ok=True)

    # overview of mash distances
    print("Plotting mash distance overview...")
    fname = figs_fld / "cl_mash_distance.png"
    mash_distance_overview(mash_df, info, args.eps, fname)

    # Refine MLST series
    mlst_series = refined_mlst_series(info, thr_size=args.thr_size)

    # benchmark clustering performance
    print("Benchmarking clustering performance...")
    eps_values = np.logspace(-5, -1.5, 100)
    clust_benchmark = clustering_benchmark_df(
        eps_values=eps_values,
        mlst_series=mlst_series,
        mash_df=mash_df,
        min_size=args.thr_size,
    )

    # Display clustering benchmark results
    print("Plotting clustering benchmark results...")
    fname = figs_fld / "cl_clustering_benchmark.png"
    plot_clustering_benchmark(clust_benchmark, args.eps, args.thr_size, fname)

    # calculate best clustering
    clust_res, clust_df = evaluate_clustering(
        eps=args.eps, mlst_series=mlst_series, mash_df=mash_df, min_size=args.thr_size
    )

    # plot confusion matrix
    print("Plotting confusion matrix...")
    fname = figs_fld / "cl_confusion_matrix.png"
    plot_confusion_matrix(
        cluster_df=clust_df,
        cluster_res=clust_res,
        eps=args.eps,
        svname=fname,
    )

    # plot phylogenetic tree with MLST and cluster labels
    print("Plotting phylogenetic tree with MLST and cluster labels...")
    fname = figs_fld / "cl_tree_clustering.png"
    plot_tree_clustering(tree, clust_df, fname)

    # save clustering results
    print("Saving clustering results...")
    save_clustering_results(
        clust_df=clust_df,
        info_df=info,
        mash_df=mash_df,
        clust_res=clust_res,
        species=args.species,
        assignment_threshold=args.assignment_threshold_freq,
        fld=clusters_fld,
    )

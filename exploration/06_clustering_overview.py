# %%
import json
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_clusters(base_dir, subfolder, glob_pattern="ST_*_info.json"):
    """
    Parse ST cluster JSON files for all species and sequence types.
    """
    results = {}

    # Get all species directories
    base_path = Path(base_dir)
    if not base_path.exists():
        print(f"Base directory {base_dir} does not exist")
        return results

    for species_dir in base_path.iterdir():
        if species_dir.is_dir():
            species = species_dir.name
            results[species] = {}

            # Look for ST_clusters directory
            clusters_dir = species_dir / subfolder
            if clusters_dir.exists():
                # Parse all ST JSON files
                for json_file in clusters_dir.glob(glob_pattern):
                    cl_name = json_file.stem.replace("_info", "")

                    try:
                        with open(json_file, "r") as f:
                            data = json.load(f)
                        results[species][cl_name] = data
                        print(f"Parsed {species}/{cl_name}")
                    except Exception as e:
                        print(f"Error parsing {json_file}: {e}")

    return results


# Parse all ST cluster files
st_data = parse_clusters("../clusters", "ST_clusters", "ST_*_info.json")

# Convert to DataFrame for easier manipulation
cl_records = []
for species, st_dict in st_data.items():
    for st, info in st_dict.items():
        record = {
            "species": species,
            "ST": st,
            "n_isolates": info.get("n_isolates", 0),
            "n_ST_isolates": info.get("n_ST_isolates", 0),
            "mean_mash_dist": info.get("mash_dist_avg", 0.0),
            "std_mash_dist": info.get("mash_dist_std", 0.0),
            "precision": info["precision"],
            "recall": info["recall"],
        }
        cl_records.append(record)
st_df = pd.DataFrame(cl_records)
st_df
# %%
fig, ax = plt.subplots(figsize=(8, 4))
sns.scatterplot(
    data=st_df,
    x="n_isolates",
    y="mean_mash_dist",
    hue="species",
    style="species",
)
plt.xlabel("Number of isolates")
plt.ylabel("Mean MASH distance")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figs/f06_clustering_overview_st.png", dpi=300)
plt.show()

# %%

# Calculate number of clusters per species
cluster_counts = (
    st_df.groupby("species")
    .size()
    .reset_index(name="n_clusters")
    .sort_values("n_clusters", ascending=False)
)
species_order = cluster_counts["species"].tolist()


def create_boxplot_with_stripplot(data, species_order, title, filename):
    """
    Create a boxplot with stripplot overlay for isolate counts by species.
    """
    sns.boxplot(
        data=data,
        y="species",
        x="n_isolates",
        width=0.5,
        fliersize=0,
        color="lightgray",
        order=species_order,
    )
    sns.stripplot(
        data=data,
        y="species",
        x="n_isolates",
        jitter=0.3,
        alpha=0.8,
        hue="species",
        marker="o",
        edgecolor="black",
        linewidth=0.5,
        order=species_order,
    )


def add_cluster_count_annotations(cluster_counts):
    """
    Add cluster count annotations to the right side of the plot.
    """
    for i, (species, count) in enumerate(cluster_counts.values):
        plt.text(
            plt.xlim()[1] * 0.8,
            i,
            f"n={count}",
            verticalalignment="center",
            fontsize=10,
            fontweight="bold",
        )


def format_boxplot(title, filename):
    """
    Apply common formatting to boxplot figures.
    """
    plt.grid(True, which="major", axis="x", alpha=0.6, linewidth=1.0)
    plt.grid(True, which="minor", axis="x", alpha=0.3, linewidth=0.5)
    plt.xlabel("Number of isolates")
    plt.ylabel("Species")
    plt.xscale("log")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()


create_boxplot_with_stripplot(
    st_df,
    species_order,
    "Number of isolates per species in ST clusters",
    "figs/f06_clustering_overview_st_isolates.png",
)
add_cluster_count_annotations(cluster_counts)
format_boxplot(
    "Number of isolates per species in ST clusters",
    "figs/f06_clustering_overview_st_isolates.png",
)
# %%
# Create precision-recall plot for ST clusters
fig, ax = plt.subplots(figsize=(8, 6))

sns.scatterplot(
    data=st_df,
    x="recall",
    y="precision",
    hue="species",
    style="species",
    s=60,
    alpha=0.7,
)

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.xlim(-0.05, 1.05)  # Set x-axis limits to include 0 and 1
plt.ylim(-0.05, 1.05)  # Set y-axis limits
plt.title("Precision-Recall for ST Clusters")
plt.grid(True, alpha=0.3)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("figs/f06_clustering_overview_st_precision_recall.png", dpi=300)
plt.show()

# %%

# Parse all ST cluster files
cl_data = parse_clusters("../clusters", "hdbscan_clusters", "CL_*_info.json")

# Convert to DataFrame for easier manipulation
cl_records = []
for species, cl_dict in cl_data.items():
    for cl, info in cl_dict.items():
        max_iso = info["ST_most_common"]
        record = {
            "species": info["species"],
            "CL": cl,
            "ST": info["ST_assignment"],
            "n_isolates": info["num_isolates"],
            "mash_dist_avg": info["mash_dist_avg"],
            "mash_dist_std": info["mash_dist_std"],
            "assigned": True if info["ST_assignment"] is not None else False,
            "ST_max_freq": info["ST_most_common_freq"],
            "precision": info["precision"],
            "recall": info["recall"],
        }
        cl_records.append(record)
cl_df = pd.DataFrame(cl_records)
cl_df

# %%
fig, ax = plt.subplots(figsize=(8, 4))
sns.scatterplot(
    data=cl_df,
    x="n_isolates",
    y="mash_dist_avg",
    hue="species",
    style="species",
    size="assigned",
    sizes=(50, 20),
)
plt.xlabel("Number of isolates")
plt.ylabel("Mean MASH distance")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figs/f06_clustering_overview_hdbscan.png", dpi=300)
plt.show()

# %%
# Calculate number of clusters per species
cluster_counts = (
    cl_df.groupby("species")
    .size()
    .reset_index(name="n_clusters")
    .sort_values("n_clusters", ascending=False)
)
species_order = cluster_counts["species"].tolist()

create_boxplot_with_stripplot(
    cl_df,
    species_order,
    "Number of isolates per species in HDBSCAN clusters",
    "figs/f06_clustering_overview_hdbscan_isolates.png",
)
add_cluster_count_annotations(cluster_counts)
format_boxplot(
    "Number of isolates per species in HDBSCAN clusters",
    "figs/f06_clustering_overview_hdbscan_isolates.png",
)

# %%
sns.scatterplot(
    data=cl_df,
    x="n_isolates",
    y="ST_max_freq",
    hue="species",
    style="species",
    size="assigned",
    sizes=(50, 20),
)
plt.xlabel("Number of isolates")
plt.ylabel("Max ST frequency in cluster")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xscale("log")
plt.grid(True, alpha=0.3)
plt.grid(True, which="major", axis="x", alpha=0.6, linewidth=1.0)
plt.grid(True, which="minor", axis="x", alpha=0.3, linewidth=0.5)
plt.tight_layout()
plt.savefig("figs/f06_clustering_overview_hdbscan_ST_max_freq.png", dpi=300)
plt.show()

# %%
# Create precision-recall plot
fig, ax = plt.subplots(figsize=(8, 6))

# Filter out rows where recall is None
plot_data = cl_df.dropna(subset=["recall"])

sns.scatterplot(
    data=plot_data,
    x="recall",
    y="precision",
    hue="species",
    style="species",
    s=60,
    alpha=0.7,
)

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.xlim(-0.05, 1.05)  # Set x-axis limits to include 0 and 1
plt.ylim(-0.05, 1.05)  # Set y-axis limits
plt.axvline(0.6, color="gray", linestyle="--")
plt.axhline(0.6, color="gray", linestyle="--")
plt.title("Precision-Recall for HDBSCAN Clusters")
plt.grid(True, alpha=0.3)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("figs/f06_clustering_overview_hdbscan_precision_recall.png", dpi=300)
plt.show()

# %%

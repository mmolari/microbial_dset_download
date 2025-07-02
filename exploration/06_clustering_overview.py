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
            "ST_completeness": info.get("ST_completeness", 0.0),
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
    size="ST_completeness",
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

# Parse all ST cluster files
cl_data = parse_clusters("../clusters", "hdbscan_clusters", "CL_*_info.json")

# Convert to DataFrame for easier manipulation
cl_records = []
for species, cl_dict in cl_data.items():
    for cl, info in cl_dict.items():
        record = {
            "species": info["species"],
            "CL": cl,
            "ST": info["ST_assignment"],
            "n_isolates": info["num_isolates"],
            "mash_dist_avg": info["mash_dist_avg"],
            "mash_dist_std": info["mash_dist_std"],
            "assigned": len(info["ST_assignment"]) > 0,
            "ST_max_freq": info["ST_max_freq"],
        }
        cl_records.append(record)
cl_df = pd.DataFrame(cl_records)
cl_df

# %%
fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(
    data=cl_df,
    x="n_isolates",
    y="mash_dist_avg",
    hue="species",
    size="ST_max_freq",
    style="species",
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

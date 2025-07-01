# %%
import json
import os
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_st_clusters(base_dir):
    """
    Parse ST cluster JSON files for all species and sequence types.

    Args:
        base_dir (str): Base directory containing species folders

    Returns:
        dict: Nested dictionary with structure {species: {st: data}}
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
            st_clusters_dir = species_dir / "ST_clusters"
            if st_clusters_dir.exists():
                # Parse all ST JSON files
                for json_file in st_clusters_dir.glob("ST_*_info.json"):
                    st_name = json_file.stem.replace("_info", "")

                    try:
                        with open(json_file, "r") as f:
                            data = json.load(f)
                        results[species][st_name] = data
                        print(f"Parsed {species}/{st_name}")
                    except Exception as e:
                        print(f"Error parsing {json_file}: {e}")

    return results


# Parse all ST cluster files
st_data = parse_st_clusters("../clusters")

# Convert to DataFrame for easier manipulation
st_records = []
for species, st_dict in st_data.items():
    for st, info in st_dict.items():
        record = {
            "species": species,
            "ST": st,
            "n_isolates": info.get("n_isolates", 0),
            "n_ST_isolates": info.get("n_ST_isolates", 0),
            "mean_mash_dist": info.get("mean_mash_dist", 0.0),
            "std_mash_dist": info.get("std_mash_dist", 0.0),
            "ST_completeness": info.get("ST_completeness", 0.0),
        }
        st_records.append(record)
st_df = pd.DataFrame(st_records)
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
plt.savefig("figs/f06_clustering_overview.png")
plt.show()

# %%

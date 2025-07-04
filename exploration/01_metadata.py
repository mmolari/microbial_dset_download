# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# %%
def load_metadata(fname):
    """Load metadata dataframe from TSV file."""
    df_meta = pd.read_csv(
        fname, sep="\t", index_col=0, parse_dates=["Assembly Release Date"]
    )
    return df_meta


def load_mlst(fname):
    """Load MLST dataframe from TSV file."""
    df_mlst = pd.read_csv(fname, sep="\t", header=None)
    # first column: take the string and keep only the part after '/' and before '.fa'
    df_mlst[0] = df_mlst[0].str.split("/").str[-1].str.split(".fa").str[0]
    # keep only first two columns
    df_mlst = df_mlst[[0, 2]]
    df_mlst.columns = ["Assembly Name", "MLST"]
    # set index to "Assembly Name"
    df_mlst.set_index("Assembly Name", inplace=True)
    return df_mlst


def load_chromosome_mapping(fname):
    """Load assembly <-> chromosome mapping dataframe from TSV file."""
    df_chrom = pd.read_csv(fname, sep="\t", index_col=0)
    return df_chrom


def decorate_metadata(df_meta, df_chrom, df_mlst=None):
    """Add chromosome accession and MLST to metadata dataframe."""
    # Add chromosome accession by merging
    df_meta = df_meta.merge(
        df_chrom,
        left_index=True,
        right_index=True,
        how="right",
        validate="1:1",
    )

    # Add MLST if provided
    if df_mlst is not None:
        # merge MLST on "Assembly Name" and df_meta on "chromosome_acc"
        # check for 1-1 mapping
        df_meta = df_meta.merge(
            df_mlst,
            left_on="chromosome_acc",
            right_index=True,
            how="right",
            validate="1:1",
        )

    # Map all "not applicable", "Missing", "unknown" and "NA" values to np.nan
    df_meta = df_meta.replace(
        to_replace=[
            "not applicable",
            "Missing",
            "unknown",
            "NA",
            "Not Applicable",
            "missing",
            "Unknown",
            "na",
            "-",
            "N/A",
            "n/a",
            "n.a.",
            "n.a",
            "N.A.",
            "not collected",
        ],
        value=np.nan,
    )

    return df_meta


species = "saureus"

# Load data
df_meta = load_metadata(f"../data/species/{species}/info.tsv")
df_mlst = load_mlst(f"../data/species/{species}/mlst.tsv")
df_chrom = load_chromosome_mapping(f"../data/species/{species}/assembly_to_chrom.tsv")

# Decorate metadata
df_meta = decorate_metadata(df_meta, df_chrom, df_mlst)
df_meta


# %%
# Combined 4-panel plot
fig, axes = plt.subplots(
    2, 2, figsize=(12, 12), gridspec_kw={"height_ratios": [1, 1.7]}
)

# Add general title with total number of assemblies
fig.suptitle(f"{species} - {len(df_meta)} assemblies", fontsize=14, y=0.98)

# 1. Assembly release date histogram
ax = axes[0, 0]
sns.histplot(df_meta["Assembly Release Date"].dt.year, discrete=True, ax=ax)
ax.set_title("Assembly Release Date Distribution")
ax.set_xlabel("Year")
ax.set_ylabel("Count")
ax.tick_params(axis="x")

# 2. Assembly N50 histogram
ax = axes[0, 1]
sns.ecdfplot(df_meta["Assembly Stats Contig N50"] / 1e6, ax=ax, stat="count")
ax.set_title("Assembly N50 Distribution")
ax.set_xlabel("N50 (Mbp)")
ax.set_ylabel("Count")

# 3. MLST distribution
ax = axes[1, 0]
top_30_mlst = df_meta["MLST"].value_counts().head(30).index
sns.countplot(data=df_meta, y="MLST", order=top_30_mlst, ax=ax)
ax.set_title("MLST Distribution (Top 30)")
ax.set_ylabel("MLST")
ax.set_xlabel("Count")

# 4. Countries distribution
ax = axes[1, 1]
df_meta_countries = df_meta[["Assembly BioSample Geographic location"]].copy()
df_meta_countries["Assembly BioSample Geographic location"] = (
    df_meta_countries["Assembly BioSample Geographic location"]
    .str.split(":", n=1)
    .str[0]
)
top_30_countries = (
    df_meta_countries["Assembly BioSample Geographic location"]
    .value_counts()
    .head(30)
    .index
)
sns.countplot(
    data=df_meta_countries,
    y="Assembly BioSample Geographic location",
    order=top_30_countries,
    ax=ax,
)
ax.set_title("Country Distribution (Top 30)")
ax.set_ylabel("Country")
ax.set_xlabel("Count")

plt.tight_layout()
plt.show()
# %%

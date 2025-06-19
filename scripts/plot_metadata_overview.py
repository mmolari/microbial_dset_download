import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot metadata overview from assembly data"
    )
    parser.add_argument(
        "--df", required=True, help="Input CSV file containing metadata"
    )
    parser.add_argument("--species", required=True, help="Species name")
    parser.add_argument("--fig", required=True, help="Output figure filename")
    return parser.parse_args()


def plot_figure(df, species, fig_filename):
    # Combined 4-panel plot
    fig, axes = plt.subplots(
        2, 2, figsize=(12, 12), gridspec_kw={"height_ratios": [1, 1.7]}
    )

    # Add general title with total number of assemblies
    fig.suptitle(f"{species} - {len(df)} assemblies", fontsize=14, y=0.98)

    # 1. Assembly release date histogram
    ax = axes[0, 0]
    sns.histplot(df["Assembly Release Date"].dt.year, discrete=True, ax=ax)
    ax.set_title("Assembly Release Date Distribution")
    ax.set_xlabel("Year")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x")

    # 2. Assembly N50 histogram
    ax = axes[0, 1]
    sns.ecdfplot(df["Assembly Stats Contig N50"] / 1e6, ax=ax, stat="count")
    ax.set_title("Assembly N50 Distribution")
    ax.set_xlabel("N50 (Mbp)")
    ax.set_ylabel("Count")

    # 3. MLST distribution
    ax = axes[1, 0]
    df_mlst = df[["MLST"]].copy()
    df_mlst["MLST"] = pd.to_numeric(df_mlst["MLST"], errors="coerce")
    df_mlst = df_mlst.dropna(subset=["MLST"])
    df_mlst["MLST"] = df_mlst["MLST"].astype(int)
    top_30_mlst = df_mlst["MLST"].value_counts().head(30).index
    sns.countplot(data=df_mlst, y="MLST", order=top_30_mlst, ax=ax)
    ax.set_title("MLST Distribution (Top 30)")
    ax.set_ylabel("MLST")
    ax.set_xlabel("Count")

    # 4. Countries distribution
    ax = axes[1, 1]
    df_meta_countries = df[["Assembly BioSample Geographic location"]].copy()
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

    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_filename, dpi=250)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_arguments()
    df = pd.read_csv(args.df, index_col=0, parse_dates=["Assembly Release Date"])
    plot_figure(df, args.species, args.fig)

import argparse
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Midpoint-root and ladderize a Newick tree"
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree file")
    parser.add_argument("--metadata", required=True, help="Metadata dataframe file")
    parser.add_argument("--fig", required=True, help="Output figure")
    return parser.parse_args()


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
    df_mlst["MLST"] = pd.to_numeric(df_mlst["MLST"], errors="coerce")
    df_mlst = df_mlst.dropna(subset=["MLST"])
    df_mlst["MLST"] = df_mlst["MLST"].astype(int)
    return df_mlst["MLST"].value_counts().head(top_n).index.tolist()


def draw_MLST(tree, df, top_mlst, ax, positions):
    """Draw MLST labels on the tree."""
    cmap = plt.get_cmap("tab20")
    mlst_colors = {
        mlst: cmap(i % 20) for i, mlst in enumerate(top_mlst)
    }  # Assign a color to each MLST
    for clade in tree.get_terminals():
        mlst = df.loc[clade.name, "MLST"]
        if pd.notna(mlst) and int(mlst) in top_mlst:
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


def draw_countries(ax, df, positions, n_top=5, x=0):
    """Draw country labels on the tree."""

    # for countries in the form of "country:region", we only take the country part
    df["Assembly BioSample Geographic location"] = (
        df["Assembly BioSample Geographic location"].str.split(":", n=1).str[0]
    )

    top_countries = (
        df["Assembly BioSample Geographic location"]
        .value_counts()
        .head(n_top)
        .index.tolist()
    )
    country_cmap = plt.get_cmap("tab10")
    country_colors = {
        country: country_cmap(i % 10) for i, country in enumerate(top_countries)
    }  # Assign a color to each country

    for clade in tree.get_terminals():
        country = df.loc[clade.name, "Assembly BioSample Geographic location"]
        if country in top_countries:
            _, y = positions[clade.name]
            ax.plot(x, y, "_", color=country_colors[country])

    # add legend
    handles = [
        plt.Line2D(
            [],
            [],
            marker="_",
            color=country_colors[country],
            label=country,
            markersize=10,
        )
        for country in top_countries
    ]
    legend = ax.legend(
        handles=handles,
        title="Country",
        bbox_to_anchor=(1.05, 0),
        loc="lower left",
        frameon=False,
    )
    ax.add_artist(legend)


def draw_year(ax, df, positions, x=1):
    cmap = plt.get_cmap("Blues")
    years = df["Assembly Release Date"].dt.year.unique()
    year_color_norm = plt.Normalize(years.min(), years.max())
    for k, (_, y) in positions.items():
        year = df.loc[k, "Assembly Release Date"].year
        color = cmap(year_color_norm(year))
        ax.plot(x, y, "_", color=color)

    # add a colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=year_color_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", shrink=0.2)
    cbar.set_label("Assembly Release Year")


def draw_host(ax, df, positions, n_top=5, x=2):
    """Draw host labels on the tree."""
    top_hosts = df["Assembly BioSample Host"].value_counts().head(n_top).index.tolist()
    host_cmap = plt.get_cmap("tab10")
    host_colors = {
        host: host_cmap(i % 10) for i, host in enumerate(top_hosts)
    }  # Assign a color to each host

    for clade in tree.get_terminals():
        host = df.loc[clade.name, "Assembly BioSample Host"]
        if host in top_hosts:
            _, y = positions[clade.name]
            ax.plot(x, y, "_", color=host_colors[host])

    # add legend
    handles = [
        plt.Line2D(
            [],
            [],
            marker="_",
            color=host_colors[host],
            label=host,
            markersize=10,
        )
        for host in top_hosts
    ]
    ax.legend(
        handles=handles,
        title="Host",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        frameon=False,
    )


def plot_tree(tree, df, svname):
    top_MLST = get_top_MLST(df, top_n=20)
    positions = get_xy_positions(tree)

    fig, axs = plt.subplots(
        1,
        2,
        figsize=(10, 15),
        sharey=True,
        gridspec_kw={"width_ratios": [4, 1]},
    )
    ax = axs[0]
    Phylo.draw(tree, do_show=False, label_func=lambda x: "", axes=ax)
    draw_MLST(tree, df, top_MLST, ax, positions)

    ax = axs[1]
    sns.despine(ax=ax, left=True, bottom=True)
    draw_countries(ax, df, positions, n_top=5, x=0)
    draw_year(ax, df, positions, x=1)
    draw_host(ax, df, positions, n_top=5, x=2)

    # remove horizontal spacing between subplots
    plt.subplots_adjust(wspace=0.05)
    # set x-axis ticks
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Country", "Year", "Host"])
    # set horizontal grid lines
    for ax in axs:
        ax.yaxis.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(svname, bbox_inches="tight", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_arguments()

    # load metadata (not used in tree operations but parsed as input)
    df = pd.read_csv(
        args.metadata,
        index_col=["chromosome_acc"],
        parse_dates=["Assembly Release Date"],
    )

    # read tree
    tree = Phylo.read(args.tree, "newick")

    # plot the tree with metadata
    plot_tree(tree, df, args.fig)

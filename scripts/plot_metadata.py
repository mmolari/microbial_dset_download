import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=str)
    parser.add_argument("--out_fig", type=str)
    return parser.parse_args()


def load_metadata(fname):
    df = pd.read_csv(fname, sep="\t", low_memory=False)

    # interpret as year
    lab = "Assembly BioSample Collection date "
    df[lab] = pd.to_datetime(df[lab], errors="coerce")

    # location
    lab = "Assembly BioSample Geographic location "
    nlab = "Geographic Location"
    df[nlab] = df[lab].str.split(":").str[0].str.strip()
    df[nlab] = df[nlab].astype("category")

    mp = {
        "unknown": None,
        "Unknown": None,
        "not collected": None,
        "not applicable": None,
        "Not Applicable": None,
        "not determined": None,
        "missing": None,
        "Missing": None,
        "-": None,
    }
    df[nlab] = df[nlab].replace(mp)

    return df


if __name__ == "__main__":

    args = parse_args()

    df = load_metadata(args.metadata)

    fig, axs = plt.subplots(1, 2, figsize=(8, 9))

    ax = axs[0]
    lab = "Assembly BioSample Collection date "
    sns.histplot(data=df, y=lab, bins=30, ax=ax)
    ax.set_ylabel("collection date")

    ax = axs[1]
    lab = "Geographic Location"
    cat_ord = df[lab].value_counts().sort_values(ascending=False).index
    sns.countplot(data=df, y=lab, ax=ax, order=cat_ord)
    ax.set_xscale("log")

    sns.despine()
    plt.tight_layout()
    plt.savefig(args.out_fig)
    plt.close()

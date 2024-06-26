import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlst", type=str)
    parser.add_argument("--out_fig", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    df = ut.load_mlst(args.mlst)

    mst = df["ST"].value_counts().index[:20]
    df["main ST"] = df["ST"].map(lambda x: str(x) if x in mst else "other")
    mst = list(mst) + ["other"]

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.countplot(df, y="main ST", order=mst)
    plt.xlabel("number of isolates")
    plt.title("Top 20 STs")
    # plt.ylabel("ST")
    sns.despine()
    plt.tight_layout()
    plt.savefig(args.out_fig)
    plt.close()

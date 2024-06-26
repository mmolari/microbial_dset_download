import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlst", type=str)
    parser.add_argument("--mash", type=str)
    parser.add_argument("--ST_num", type=str)
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--out_fig", type=str)
    return parser.parse_args()


def mash_dist_plot(df, threshold, st, out_fig):

    mc_ST = df["ST"].value_counts().index[:5]
    mc_ST = set(mc_ST) - {"-"}
    df["main ST"] = df["ST"].map(lambda x: str(x) if x in mc_ST else "other")

    mask_thr = df["dist"] <= threshold
    mask_st = df["ST"] == st

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.histplot(
        df,
        x="dist",
        bins=100,
        hue="main ST",
        multiple="stack",
        element="step",
        linewidth=0.5,
    )

    ax.axvline(threshold, color="red", linestyle="--")
    ax.text(0.05, 0.9, f"N dist = {mask_thr.sum()} / {len(df)}", transform=ax.transAxes)
    ax.text(0.05, 0.8, f"N ST = {mask_st.sum()} / {len(df)}", transform=ax.transAxes)
    ax.text(
        0.05,
        0.7,
        f"N ST & dist = {(mask_thr & mask_st).sum()} / {len(df)}",
        transform=ax.transAxes,
    )

    plt.xlabel("mash distance")
    plt.title(f"Mash distances to ST{st} ref genome")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out_fig)
    plt.close()


if __name__ == "__main__":

    args = parse_args()

    df = ut.load_mash_mlst_df(args.mlst, args.mash)

    mash_dist_plot(df, args.threshold, args.ST_num, args.out_fig)

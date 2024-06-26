# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(exist_ok=True, parents=True)


fname = "../data/mash_dist/saureus/GCF_000009665.1.tsv"
df = pd.read_csv(fname, sep="\t", header=None)
df.columns = ["ref", "qry", "dist", "pval", "hashes"]
df["qry"] = df["qry"].str.split("/").str[-1]
df["ref"] = df["ref"].str.split("/").str[-1]

tp = fname.split("/")[-1].split(".")[0]

fname = "../data/species/saureus_mlst.tsv"
mlst = pd.read_csv(fname, sep="\t", header=None)
mlst.columns = ["qry", "mlst", "ST"] + [f"g{i+1}" for i in range(7)]
mlst["qry"] = mlst["qry"].str.split("/").str[-1]

df = df.merge(mlst[["qry", "ST"]], on="qry", validate="1:1")

# %%

sns.countplot(df, y="ST", order=df["ST"].value_counts().index[:20])

# %%


mc_ST = df["ST"].value_counts().index[:5]
df["main ST"] = df["ST"].map(lambda x: str(x) if x in mc_ST else "other")
df["main ST"]

# %%

threshold = 0.005
N = (df["dist"] <= threshold).sum()

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
ax.text(0.05, 0.9, f"N = {N} / {len(df)}", transform=ax.transAxes)

plt.xlabel("mash distance")
plt.title(f"Mash distances to {tp} ref genome")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"{tp}_mash_dist_hist.png")
plt.show()


# %%
mask = (df["ST"] == "5") & (df["dist"] < threshold)
mask.sum()

# %%

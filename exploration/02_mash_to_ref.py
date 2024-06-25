# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(exist_ok=True, parents=True)


fname = "../results/mash_dist/ST5.tsv"
df = pd.read_csv(fname, sep="\t", header=None)
df.columns = ["query", "ref", "dist", "pval", "hashes"]
df["query"] = df["query"].str.split("/").str[-1]
df["ref"] = df["ref"].str.split("/").str[-1]

tp = fname.split("/")[-1].split(".")[0]
# %%

threshold = 0.005
N = (df["dist"] <= threshold).sum()

fig, ax = plt.subplots(figsize=(8, 4))
sns.histplot(df["dist"], bins=100)

ax.axvline(threshold, color="red", linestyle="--")
ax.text(0.05, 0.9, f"N = {N} / {len(df)}", transform=ax.transAxes)

plt.xlabel("mash distance")
plt.title(f"Mash distances to {tp} ref genome")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"{tp}_mash_dist_hist.png")
plt.show()


# %%

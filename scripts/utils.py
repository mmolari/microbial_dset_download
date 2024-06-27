import pandas as pd


def load_mash(df_fname):
    df = pd.read_csv(df_fname, sep="\t", header=None)
    df.columns = ["ref", "qry", "dist", "pval", "hashes"]
    df["qry"] = df["qry"].str.split("/").str[-1]
    df["ref"] = df["ref"].str.split("/").str[-1]
    return df


def load_mlst(df_fname):
    df = pd.read_csv(df_fname, sep="\t", header=None)
    C = len(df.columns)
    df.columns = ["qry", "mlst", "ST"] + [f"g{i+1}" for i in range(C - 3)]
    df["qry"] = df["qry"].str.split("/").str[-1]
    return df


def load_mash_mlst_df(mlst_fname, mash_fname):
    mlst_df = load_mlst(mlst_fname)
    mash_df = load_mash(mash_fname)
    mash_df = mash_df.merge(mlst_df[["qry", "ST"]], on="qry", validate="1:1")
    return mash_df

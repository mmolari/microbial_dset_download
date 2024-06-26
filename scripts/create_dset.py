import argparse
import pandas as pd
import utils as ut
import shutil
import pathlib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlst", type=str)
    parser.add_argument("--mash", type=str)
    parser.add_argument("--mtd", type=str)
    parser.add_argument("--acc_map", type=str)
    parser.add_argument("--fas", type=str)
    parser.add_argument("--ST", type=str)
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--out_dir", type=str)
    parser.add_argument("--out_mtd", type=str)
    parser.add_argument("--out_mlst", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    df = ut.load_mash_mlst_df(mlst_fname=args.mlst, mash_fname=args.mash)

    mask = (df["ST"] == args.ST) & (df["dist"] < args.threshold)
    df = df[mask]

    # export df for selected isolates
    df.to_csv(args.out_mlst, sep="\t", index=False)

    isolates = df["qry"].unique()

    # copy fasta files in destination folder
    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    src_dir = pathlib.Path(args.fas)

    for isolate in isolates:
        src = src_dir / isolate
        dst = out_dir / isolate
        shutil.copy(src, dst)

    # create metadata sub-index
    acc_df = pd.read_csv(args.acc_map, sep="\t")
    acc_df.set_index("chromosome_acc", inplace=True)
    assembly_acc = acc_df.loc[isolates, "assembly_acc"].values

    # load and select metadata
    mtd = pd.read_csv(args.mtd, sep="\t")
    mtd.set_index("Assembly Accession", inplace=True)
    mtd = mtd.loc[assembly_acc]

    mtd.to_csv(args.out_mtd, sep="\t")

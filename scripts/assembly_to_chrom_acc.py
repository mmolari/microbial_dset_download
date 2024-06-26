import argparse
from Bio import SeqIO
import pathlib
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fas", type=str)
    parser.add_argument("--out", type=str)
    return parser.parse_args()


def get_fa_acc(fa_fname):
    acc = pathlib.Path(fa_fname).stem
    acc = "_".join(acc.split("_")[:-2])
    return acc


def fa_first_entry_acc(fa_fname):
    with open(fa_fname) as f:
        entry = SeqIO.read(f, "fasta")
        acc = entry.id
    return acc


if __name__ == "__main__":
    args = parse_args()

    fa_fld = pathlib.Path(args.fas)

    df = []
    for fa_fname in fa_fld.glob("*.fna"):
        acc = get_fa_acc(fa_fname)
        first_acc = fa_first_entry_acc(fa_fname)
        df.append([acc, first_acc])
    df = pd.DataFrame(df, columns=["assembly_acc", "chromosome_acc"])
    df.to_csv(args.out, sep="\t", index=False)

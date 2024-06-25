import argparse
import pathlib
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Extract chromosome from a genome")
    parser.add_argument("--fa", type=str, help="input genome file")
    parser.add_argument("--out_dir", type=str, help="output chromosome file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    records = []
    with open(args.fa, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            records.append(record)

    # select record of maximal length
    r = max(records, key=lambda x: len(x.seq))

    # format id and description
    r.id = r.id.split()[0]
    r.description = ""
    fname = pathlib.Path(args.out_dir) / f"{r.id}.fa"
    with open(fname, "w") as f:
        SeqIO.write(r, f, "fasta")

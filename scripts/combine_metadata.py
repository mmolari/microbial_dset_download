import pandas as pd
import numpy as np
import argparse


def load_metadata(fname):
    """Load metadata dataframe from TSV file."""
    df_meta = pd.read_csv(
        fname, sep="\t", index_col=0, parse_dates=["Assembly Release Date"]
    )
    return df_meta


def load_mlst(fname):
    """Load MLST dataframe from TSV file."""
    df_mlst = pd.read_csv(fname, sep="\t", header=None)
    # first column: take the string and keep only the part after '/' and before '.fa'
    df_mlst[0] = df_mlst[0].str.split("/").str[-1].str.split(".fa").str[0]
    # keep only first two columns
    df_mlst = df_mlst[[0, 2]]
    df_mlst.columns = ["Assembly Name", "MLST"]
    # set index to "Assembly Name"
    df_mlst.set_index("Assembly Name", inplace=True)
    return df_mlst


def load_chromosome_mapping(fname):
    """Load assembly <-> chromosome mapping dataframe from TSV file."""
    df_chrom = pd.read_csv(fname, sep="\t", index_col=0)
    return df_chrom


def decorate_metadata(df_meta, df_chrom, df_mlst=None):
    """Add chromosome accession and MLST to metadata dataframe."""
    # Add chromosome accession by merging
    df_meta = df_meta.merge(
        df_chrom,
        left_index=True,
        right_index=True,
        how="right",
        validate="1:1",
    )

    # Add MLST if provided
    if df_mlst is not None:
        # merge MLST on "Assembly Name" and df_meta on "chromosome_acc"
        # check for 1-1 mapping
        df_meta = df_meta.merge(
            df_mlst,
            left_on="chromosome_acc",
            right_index=True,
            how="right",
            validate="1:1",
        )

    # Map all "not applicable", "Missing", "unknown" and "NA" values to np.nan
    df_meta = df_meta.replace(
        to_replace=[
            "not applicable",
            "Missing",
            "unknown",
            "NA",
            "Not Applicable",
            "missing",
            "Unknown",
            "na",
            "-",
            "N/A",
            "n/a",
            "n.a.",
            "n.a",
            "N.A.",
            "not collected",
        ],
        value=np.nan,
    )

    return df_meta


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine metadata with MLST and chromosome mapping data"
    )
    parser.add_argument(
        "--metadata_file", required=True, help="Path to metadata TSV file"
    )
    parser.add_argument(
        "--chromosome_mapping_file",
        required=True,
        help="Path to chromosome mapping TSV file",
    )
    parser.add_argument("--mlst_file", help="Path to MLST TSV file (optional)")
    parser.add_argument(
        "--output", required=True, help="Output filename for combined metadata"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Load data
    df_meta = load_metadata(args.metadata_file)
    df_chrom = load_chromosome_mapping(args.chromosome_mapping_file)

    # Load MLST if provided
    df_mlst = None
    if args.mlst_file:
        df_mlst = load_mlst(args.mlst_file)

    # Combine metadata
    df_combined = decorate_metadata(df_meta, df_chrom, df_mlst)

    # Save combined metadata to output file
    df_combined.to_csv(args.output)

import numpy as np
import pandas as pd
from pathlib import Path


def parse_mash_triangle_output(file_path: str) -> pd.DataFrame:
    """
    Parse a lower-triangular matrix file produced by `mash triangle` into a pandas DataFrame.
    """
    with open(file_path, "r") as f:
        # First line: number of sequences
        n = int(f.readline().split()[-1])

        # Preallocate an n×n float array filled with zeros
        mat = np.zeros((n, n), dtype=float)
        ids = []

        for i, line in enumerate(f):
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            fp, *vals = parts
            ids.append(Path(fp).stem)

            if vals:
                # parse into 1d array of length i
                row_vals = np.fromiter(
                    (float(x) for x in vals), dtype=float, count=len(vals)
                )
                # fill lower triangle [i, 0:i] and mirror into [0:i, i]
                mat[i, : len(row_vals)] = row_vals
                mat[: len(row_vals), i] = row_vals

    # wrap in DataFrame
    return pd.DataFrame(mat, index=ids, columns=ids)


def get_xy_positions(tree):
    """Get x, y positions for each clade in the tree."""
    positions = {}
    depths = tree.depths()
    for i, clade in enumerate(tree.get_terminals()):
        positions[clade.name] = (depths[clade], i + 1)
    return positions

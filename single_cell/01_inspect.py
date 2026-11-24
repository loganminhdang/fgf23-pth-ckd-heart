#!/usr/bin/env python3
"""Report the structure of the atlas file without reading the expression matrix.

    python3 01_inspect.py data/HeartMap_V1.0.h5ad

Prints the encoding, shape and number of stored values of the expression
matrix, the available layers, the positions of the analysed genes, and every
obs column with its distinct values. Use the output to confirm that the
column names and labels set in atlas_io.py match the file.
"""

import sys

import h5py
import numpy as np

import atlas_io as io


def main(path):
    with h5py.File(path, "r") as h5:
        print("Top-level groups:", list(h5.keys()))

        X = h5["X"]
        encoding = io.attr(X, "encoding-type") or ("csr_matrix" if "indptr" in X else "array")
        print(f"\nX encoding: {encoding}")
        if isinstance(X, h5py.Group):
            shape = X.attrs.get("shape", None)
            print(f"X shape: {tuple(shape) if shape is not None else 'not recorded'}")
            print(f"Stored values: {X['data'].shape[0]:,}")
        else:
            print(f"X is dense, shape {X.shape}; 02_extract.py requires a sparse matrix")

        if "layers" in h5:
            print("Layers:", list(h5["layers"].keys()))

        genes = np.char.upper(io.gene_names(h5))
        print(f"\nGenes: {len(genes):,}")
        for g in io.GENES:
            hit = np.nonzero(genes == g)[0]
            print(f"  {g:6s} {'index ' + str(hit[0]) if len(hit) else 'NOT FOUND'}")

        print("\nObs columns:")
        n_nuclei = None
        for column in h5["obs"].keys():
            if column.startswith("_"):
                continue
            try:
                categories, codes = io.read_categorical(h5, column)
            except Exception as error:
                print(f"\n  {column}: not readable ({error})")
                continue
            n_nuclei = n_nuclei or len(codes)
            missing = int((codes < 0).sum())
            shown = ", ".join(categories[:12].tolist())
            more = f" (+{len(categories) - 12} more)" if len(categories) > 12 else ""
            print(f"\n  {column}: {len(categories)} values, {missing:,} missing")
            print(f"    {shown}{more}")
        print(f"\nNuclei: {n_nuclei:,}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])

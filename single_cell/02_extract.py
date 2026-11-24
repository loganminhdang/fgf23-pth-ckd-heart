#!/usr/bin/env python3
"""Extract the analysed genes and the nucleus annotations from the atlas into a
small cache.

The atlas (about 44 GB) is read once, in parallel blocks, and only the values
of the analysed genes are kept. The cache (a few megabytes) is the input to
03_analyse.py, so the atlas does not need to be read again.

    python3 02_extract.py data/HeartMap_V1.0.h5ad \\
        --matrices X layers/raw_count_cellbender --out data/heartmap_cache

Optional: --benchmark times 1 to 32 workers on 2% of the matrix and reports
the fastest setting for --workers.

Writes metadata.npz (obs columns as categories and codes) and one
expr_<matrix>.npz per matrix (for each gene, the indices and values of the
nuclei with a non-zero value).
"""

import argparse
import multiprocessing as mp
import os
import sys
import time

import h5py
import numpy as np

import atlas_io as io

_WORKER = {}


def _init(path, matrix_path, targets, n_genes):
    """Open the file once per worker and build the gene lookup table."""
    f = h5py.File(path, "r")
    X = f[matrix_path]
    lookup = np.zeros(n_genes, dtype=bool)
    lookup[np.asarray(targets, np.int64)] = True
    _WORKER.update(X=X, indptr=X["indptr"][:].astype(np.int64), lookup=lookup,
                   targets=np.asarray(targets, np.int64))


def _read_block(job):
    """Collect (row, value) pairs for the target genes in one range of stored values."""
    start, stop, chunk = job
    X, indptr = _WORKER["X"], _WORKER["indptr"]
    lookup, targets = _WORKER["lookup"], _WORKER["targets"]
    found = {int(g): ([], []) for g in targets}
    s = start
    while s < stop:
        e = min(s + chunk, stop)
        columns = X["indices"][s:e]
        hit = lookup[columns]
        if hit.any():
            position = np.nonzero(hit)[0]
            values = X["data"][s:e][position]
            gene = columns[position].astype(np.int64)
            rows = np.searchsorted(indptr, s + position, side="right") - 1
            for g in targets:
                m = gene == g
                if m.any():
                    found[int(g)][0].append(rows[m].astype(np.int64))
                    found[int(g)][1].append(values[m].astype(np.float32))
        s = e
    return {g: (np.concatenate(r) if r else np.zeros(0, np.int64),
                np.concatenate(v) if v else np.zeros(0, np.float32))
            for g, (r, v) in found.items()}


def read_csr(path, matrix_path, gene_index, n_genes, n_values, workers, chunk, quiet=False):
    """Read a CSR matrix in parallel and return {gene index: (rows, values)}."""
    n_blocks = max(workers * 4, 1)
    edges = np.linspace(0, n_values, n_blocks + 1).astype(np.int64)
    jobs = [(int(edges[i]), int(edges[i + 1]), chunk) for i in range(n_blocks)]
    merged = {g: ([], []) for g in gene_index}
    t0, done = time.time(), 0

    def collect(result):
        nonlocal done
        for g, (r, v) in result.items():
            merged[g][0].append(r)
            merged[g][1].append(v)
        done += 1
        if not quiet:
            elapsed = time.time() - t0
            print(f"\r  {done}/{n_blocks} blocks, {elapsed / 60:.1f} min elapsed",
                  end="", file=sys.stderr, flush=True)

    if workers == 1:
        _init(path, matrix_path, list(gene_index), n_genes)
        for job in jobs:
            collect(_read_block(job))
    else:
        context = mp.get_context("spawn" if sys.platform == "win32" else "fork")
        with context.Pool(workers, initializer=_init,
                          initargs=(path, matrix_path, list(gene_index), n_genes)) as pool:
            for result in pool.imap_unordered(_read_block, jobs):
                collect(result)
    if not quiet:
        print(file=sys.stderr)
    return {g: (np.concatenate(r) if r else np.zeros(0, np.int64),
                np.concatenate(v) if v else np.zeros(0, np.float32))
            for g, (r, v) in merged.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("h5ad")
    parser.add_argument("--matrices", nargs="+", default=["X"])
    parser.add_argument("--genes", nargs="+", default=io.GENES)
    parser.add_argument("--groupby", nargs="+", default=io.GROUPBY)
    parser.add_argument("--out", default="cache")
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count()))
    parser.add_argument("--chunk-nnz", type=int, default=100_000_000,
                        help="stored values read per step within a block")
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    with h5py.File(args.h5ad, "r") as h5:
        names = io.gene_names(h5)
        upper = np.char.upper(names)
        gene_index, genes = [], []
        for g in args.genes:
            hit = np.nonzero(upper == g.upper())[0]
            if len(hit):
                gene_index.append(int(hit[0]))
                genes.append(g)
            else:
                print(f"{g} not found in the atlas", file=sys.stderr)
        n_genes = len(names)

        if args.benchmark:
            X = h5[args.matrices[0]]
            subset = int(X["data"].shape[0] * 0.02)
            print(f"Timing on 2% of {args.matrices[0]} ({subset:,} values)", file=sys.stderr)
            best = (None, float("inf"))
            for w in [1, 2, 4, 8, 16, 32]:
                if w > mp.cpu_count():
                    break
                t = time.time()
                read_csr(args.h5ad, args.matrices[0], gene_index, n_genes, subset,
                         w, args.chunk_nnz, quiet=True)
                seconds = time.time() - t
                print(f"  {w:2d} workers: {seconds:6.1f} s; full matrix about "
                      f"{seconds * 50 / 60:.0f} min", file=sys.stderr)
                if seconds < best[1]:
                    best = (w, seconds)
            print(f"Fastest: --workers {best[0]}", file=sys.stderr)
            return

        # Metadata.
        meta = {}
        for column in args.groupby:
            categories, codes = io.read_categorical(h5, column)
            meta[f"{column}__cat"] = categories
            meta[f"{column}__code"] = codes
        lengths = {len(meta[f"{c}__code"]) for c in args.groupby}
        if len(lengths) != 1:
            raise SystemExit("obs columns differ in length")
        n_nuclei = lengths.pop()
        np.savez_compressed(os.path.join(args.out, "metadata.npz"),
                            columns=np.array(args.groupby), n_cells=n_nuclei, **meta)
        print(f"Nuclei: {n_nuclei:,}; obs columns: {', '.join(args.groupby)}", file=sys.stderr)

        # Expression, one pass per matrix.
        for matrix_path in args.matrices:
            if matrix_path not in h5:
                print(f"{matrix_path} not in the file; skipped", file=sys.stderr)
                continue
            X = h5[matrix_path]
            encoding = io.attr(X, "encoding-type") or ("csr_matrix" if "indptr" in X else "")
            n_values = int(X["data"].shape[0])
            print(f"\n{matrix_path}: {encoding}, {n_values:,} stored values, "
                  f"{args.workers} workers", file=sys.stderr)
            t0 = time.time()
            if "csc" in encoding:
                indptr = X["indptr"][:]
                hits = {g: (X["indices"][int(indptr[g]):int(indptr[g + 1])].astype(np.int64),
                            X["data"][int(indptr[g]):int(indptr[g + 1])].astype(np.float32))
                        for g in gene_index}
            else:
                hits = read_csr(args.h5ad, matrix_path, gene_index, n_genes, n_values,
                                args.workers, args.chunk_nnz)

            save = {"genes": np.array(genes), "n_cells": n_nuclei, "matrix": matrix_path}
            for g_index, g in zip(gene_index, genes):
                rows, values = hits[g_index]
                order = np.argsort(rows)
                save[f"{g}__cells"] = rows[order].astype(np.int32)
                save[f"{g}__vals"] = values[order]
            path = os.path.join(args.out, f"expr_{matrix_path.replace('/', '_')}.npz")
            np.savez_compressed(path, **save)
            print(f"  {(time.time() - t0) / 60:.1f} min; written to {path}", file=sys.stderr)
            for g in genes:
                n = len(save[f"{g}__cells"])
                print(f"    {g:6s} detected in {n:>9,} nuclei ({100 * n / n_nuclei:.3f}%)",
                      file=sys.stderr)


if __name__ == "__main__":
    main()

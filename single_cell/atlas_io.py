"""Reading the HeartMap atlas (.h5ad) and the extraction cache.

The atlas is read with h5py rather than anndata so that the expression
matrix can be read in parts without loading it into memory.
"""

import glob
import os

import h5py
import numpy as np
import pandas as pd

# Genes analysed, in the order reported in Table S5.
GENES = ["FGFR1", "FGFR2", "FGFR3", "KLB", "KL", "FGFR4", "FGF23"]

# Observation columns stored in the cache.
GROUPBY = ["cell_type__ontology_label", "donor_id", "study",
           "disease__ontology_label", "organ__ontology_label",
           "library_preparation_protocol__ontology_label", "sex"]

CELL_TYPE = "cell_type__ontology_label"
DONOR = "donor_id"
STUDY = "study"
DISEASE = "disease__ontology_label"
REGION = "organ__ontology_label"
CHEMISTRY = "library_preparation_protocol__ontology_label"
VENTRICULAR_CM = "regular ventricular cardiac myocyte"

# Matrix used for all reported tables. Detection frequencies are identical in
# X and in layers/raw_count_cellbender; mean expression is reported from X.
PRIMARY_MATRIX = "X"


# ---------------------------------------------------------------- h5ad access

def decode(values):
    return np.array([v.decode() if isinstance(v, bytes) else v for v in values])


def attr(node, key, default=""):
    value = node.attrs.get(key, default)
    return value.decode() if isinstance(value, bytes) else value


def index_key(group):
    """Name of the dataset holding the index of an obs or var group."""
    key = attr(group, "_index", "_index")
    if key in group:
        return key
    return "_index" if "_index" in group else list(group.keys())[0]


def read_strings(node):
    """Return (values as strings, mask of missing values) for an array node.

    Handles a dataset of fixed- or variable-length strings (or numbers) and
    the 'nullable-string-array' group (datasets 'values' and 'mask') that
    recent AnnData versions write for string columns and indexes.
    """
    if isinstance(node, h5py.Group):
        values = node["values"][:]
        missing = node["mask"][:].astype(bool) if "mask" in node else np.zeros(len(values), bool)
    else:
        values = node[:]
        missing = np.zeros(len(values), bool)
    if values.dtype.kind in "SO":
        values = decode(values)
    return values.astype(str), missing


def read_categorical(h5, column):
    """Return (categories, codes) for one obs column.

    Handles the current AnnData encoding (a group with 'categories' and
    'codes'), the legacy encoding (codes with categories under
    obs/__categories), plain string or numeric columns, and nullable string
    columns. A code of -1 denotes a missing value and is kept as -1.
    """
    node = h5["obs"][column]
    if isinstance(node, h5py.Group) and "categories" in node:
        return decode(node["categories"][:]).astype(str), node["codes"][:].astype(np.int32)
    legacy = f"obs/__categories/{column}"
    if isinstance(node, h5py.Dataset) and legacy in h5:
        return decode(h5[legacy][:]).astype(str), node[:].astype(np.int32)
    values, missing = read_strings(node)
    categories, inverse = np.unique(values[~missing], return_inverse=True)
    codes = np.full(len(values), -1, dtype=np.int32)
    codes[~missing] = inverse
    return categories.astype(str), codes


def gene_names(h5):
    var = h5["var"]
    return read_strings(var[index_key(var)])[0]


# ---------------------------------------------------------------- cache access

def load_metadata(cache):
    """Return the cached obs columns as a DataFrame, and the number of nuclei.

    Missing values (code -1) become NaN and are therefore excluded from
    every grouped summary.
    """
    z = np.load(os.path.join(cache, "metadata.npz"), allow_pickle=True)
    columns = [str(c) for c in z["columns"]]
    n_nuclei = int(z["n_cells"])
    frame = {c: pd.Categorical.from_codes(z[f"{c}__code"],
                                          categories=[str(x) for x in z[f"{c}__cat"]])
             for c in columns}
    return pd.DataFrame(frame), n_nuclei


def load_matrices(cache, n_nuclei):
    """Return {matrix name: (genes, {gene: dense vector over all nuclei})}."""
    matrices = {}
    for path in sorted(glob.glob(os.path.join(cache, "expr_*.npz"))):
        z = np.load(path, allow_pickle=True)
        genes = [str(g) for g in z["genes"]]
        values = {}
        for g in genes:
            v = np.zeros(n_nuclei, dtype=np.float32)
            v[z[f"{g}__cells"]] = z[f"{g}__vals"]
            values[g] = v
        matrices[str(z["matrix"])] = (genes, values)
    if not matrices:
        raise SystemExit(f"no expr_*.npz files in {cache}; run 02_extract.py first")
    return matrices


def apply_subsets(meta, matrices, subsets):
    """Restrict to nuclei matching every 'column=value' condition."""
    keep = np.ones(len(meta), dtype=bool)
    for condition in subsets:
        column, value = condition.split("=", 1)
        keep &= (meta[column].astype(str) == value.strip("\"'")).to_numpy()
    meta = meta[keep].reset_index(drop=True)
    matrices = {name: (genes, {g: v[keep] for g, v in values.items()})
                for name, (genes, values) in matrices.items()}
    return meta, matrices


def summarise(meta, values, genes, keys):
    """Detection and mean expression per group, one row per group and gene."""
    rows = []
    for name, index in meta.groupby(keys, observed=True).indices.items():
        name = name if isinstance(name, tuple) else (name,)
        for g in genes:
            v = values[g][index]
            detected = v > 0
            rows.append({**dict(zip(keys, map(str, name))), "gene": g,
                         "n_nuclei": len(v), "n_detected": int(detected.sum()),
                         "detection_pct": 100 * detected.mean(),
                         "mean_expr": float(v.mean()),
                         "mean_expr_detected": float(v[detected].mean()) if detected.any() else 0.0})
    return pd.DataFrame(rows)

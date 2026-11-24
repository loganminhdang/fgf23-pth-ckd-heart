#!/usr/bin/env python3
"""Detection tables from the cache (Supplementary Table S5, Figure 5,
Supplementary Figure S3, and the single-nucleus values in the text).

    python3 03_analyse.py data/heartmap_cache --outdir outputs/single_cell/tables

Detection frequency: percentage of nuclei in a group with any transcript of
the gene (non-zero value). Mean expression: mean value over all nuclei in the
group. All tables use the atlas's main expression matrix (X). When the
CellBender-corrected counts are also cached, detection in the two matrices is
compared (matrix_concordance.csv).

Tables written:
  table_s5.csv                     detection and mean expression by cell type (Table S5),
                                   four decimal places
  table_s5_printed.csv             Table S5 as printed in the Supplementary Information
                                   (nuclei labelled "cell" omitted; detection to three and
                                   mean expression to four decimal places, rounded from
                                   the unrounded values)
  by_celltype.csv                  the same in long format (Figure 5a)
  vcm_by_study_chemistry.csv       ventricular cardiomyocytes, by study and library chemistry (Figure 5b)
  vcm_chaffin_by_disease.csv       ventricular cardiomyocytes in Chaffin et al., by disease state (Figure 5c)
  by_study.csv, by_chemistry.csv,
  by_region.csv, by_disease.csv    all nuclei, by stratum (Figure S3a-d)
  donor_support.csv                percentage of donors with any detection, by cell type
  matrix_concordance.csv           detection in each cached matrix
  atlas_summary.csv                for each cached obs column, the number of distinct values
                                   and of nuclei with and without a value
  text_values.csv                  the single-nucleus values reported in the main text and
                                   the Figure 5 legend, unrounded (checked by 04_check_table_s5.py)

Other groupings: --by column [column ...], optionally with --subset column=value.
"""

import argparse
import os
import re
import sys

import pandas as pd

import atlas_io as io

UNRESOLVED = "cell"     # atlas label of nuclei without a resolved cell type


def wide(long, key, genes, decimals=4):
    """One row per group with detection and mean expression for each gene."""
    out = None
    for g in genes:
        part = (long[long["gene"] == g][[key, "detection_pct", "mean_expr"]]
                .rename(columns={"detection_pct": f"{g}_detect%", "mean_expr": f"{g}_mean"}))
        out = part if out is None else out.merge(part, on=key)
    return out if decimals is None else out.round(decimals)


def detection_pivot(long, index, genes):
    """Detection (%) with one row per group and one column per gene."""
    table = long.pivot_table(index=index, columns="gene", values="detection_pct", observed=True)
    return table[[g for g in genes if g in table.columns]]


def printed_table_s5(by_cell_type, genes):
    """Table S5 in the layout of the Supplementary Information, formatted from
    the unrounded values: detection (%) to three and mean expression to four
    decimal places; nuclei without a resolved cell type are omitted."""
    t = wide(by_cell_type, io.CELL_TYPE, genes, decimals=None)
    t = t[t[io.CELL_TYPE].astype(str) != UNRESOLVED].sort_values(io.CELL_TYPE)
    out = pd.DataFrame({"Cell type": [c[:1].upper() + c[1:] for c in t[io.CELL_TYPE].astype(str)]})
    for g in genes:
        out[f"{g} detected (%)"] = [f"{v:.3f}" for v in t[f"{g}_detect%"]]
    for g in genes:
        out[f"{g} mean"] = [f"{v:.4f}" for v in t[f"{g}_mean"]]
    return out


def atlas_summary(meta, n_nuclei):
    """Distinct values and nuclei with and without a value, per cached obs column."""
    rows = []
    for column in meta.columns:
        present = int(meta[column].notna().sum())
        rows.append({"column": column, "distinct_values": int(meta[column].nunique()),
                     "nuclei_with_value": present, "nuclei_missing": n_nuclei - present})
    return pd.DataFrame(rows)


def text_values(meta, n_nuclei, by_cell_type, chaffin_by_disease):
    """Values reported in the main text and the Figure 5 legend, unrounded."""
    rows = []

    def add(item, value, description, source):
        rows.append({"item": item, "value": value, "description": description, "source": source})

    add("atlas_nuclei", n_nuclei, "nuclei in the atlas", "atlas_summary.csv, nuclei_with_value + nuclei_missing")
    for item, column, what in [("atlas_donors", io.DONOR, "individuals (distinct donor_id)"),
                               ("atlas_studies", io.STUDY, "constituent studies"),
                               ("atlas_regions", io.REGION, "cardiac regions"),
                               ("atlas_disease_states", io.DISEASE, "disease states"),
                               ("atlas_chemistries", io.CHEMISTRY, "library chemistries")]:
        add(item, int(meta[column].nunique()), what, f"atlas_summary.csv, {column}")

    vcm = by_cell_type[by_cell_type[io.CELL_TYPE].astype(str) == io.VENTRICULAR_CM].set_index("gene")
    if len(vcm):
        add("vcm_nuclei", int(vcm["n_nuclei"].iloc[0]), "ventricular cardiomyocyte nuclei",
            "by_celltype.csv, n_nuclei")
        for g in vcm.index:
            add(f"vcm_detection_{g}", float(vcm.loc[g, "detection_pct"]),
                f"{g}, % of ventricular cardiomyocyte nuclei detected", "by_celltype.csv, detection_pct")
        for g in ("FGFR1", "FGF23"):
            if g in vcm.index and "FGFR4" in vcm.index and vcm.loc["FGFR4", "n_detected"] > 0:
                add(f"vcm_ratio_{g}_to_FGFR4",
                    float(vcm.loc[g, "n_detected"] / vcm.loc["FGFR4", "n_detected"]),
                    f"{g} detection relative to FGFR4 detection, ventricular cardiomyocytes",
                    "by_celltype.csv, n_detected")

    for label, short in [("normal", "normal"), ("hypertrophic cardiomyopathy", "hcm"),
                         ("dilated cardiomyopathy", "dcm")]:
        part = chaffin_by_disease[(chaffin_by_disease[io.DISEASE].astype(str).str.lower() == label)
                                  & (chaffin_by_disease["gene"] == "FGFR4")]
        if len(part):
            add(f"chaffin_vcm_detection_FGFR4_{short}", float(part["detection_pct"].iloc[0]),
                f"FGFR4, % of ventricular cardiomyocyte nuclei detected, Chaffin et al., {label}",
                "vcm_chaffin_by_disease.csv, detection_pct")
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cache")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--matrix", default=io.PRIMARY_MATRIX)
    parser.add_argument("--by", nargs="+", default=None)
    parser.add_argument("--subset", nargs="*", default=[])
    args = parser.parse_args()
    out = args.outdir or os.path.join(args.cache, "tables")
    os.makedirs(out, exist_ok=True)

    meta, n_nuclei = io.load_metadata(args.cache)
    matrices = io.load_matrices(args.cache, n_nuclei)
    if args.matrix not in matrices:
        raise SystemExit(f"matrix {args.matrix} not cached; available: {list(matrices)}")
    meta, matrices = io.apply_subsets(meta, matrices, args.subset)
    cached_genes, values = matrices[args.matrix]
    genes = [g for g in io.GENES if g in cached_genes]
    print(f"{len(meta):,} nuclei; matrix {args.matrix}", file=sys.stderr)

    # Custom grouping.
    if args.by:
        long = io.summarise(meta, values, genes, args.by)
        tag = "_".join(args.by + [re.sub(r"\W+", "-", s) for s in args.subset])
        long.to_csv(os.path.join(out, f"custom_{tag}.csv"), index=False)
        print(detection_pivot(long, args.by, genes).round(4).to_string())
        return

    # Structure of the cached metadata (of the subset, if --subset is given).
    summary = atlas_summary(meta, len(meta))
    summary.to_csv(os.path.join(out, "atlas_summary.csv"), index=False)
    print("Cached obs columns\n" + summary.to_string(index=False) + "\n")

    # Concordance of detection between matrices.
    if len(matrices) > 1:
        rows = [{"gene": g, **{name: 100 * (vals[g] > 0).mean()
                               for name, (_, vals) in matrices.items()}} for g in genes]
        concordance = pd.DataFrame(rows).round(4)
        concordance.to_csv(os.path.join(out, "matrix_concordance.csv"), index=False)
        print("Detection (%) by matrix\n" + concordance.to_string(index=False) + "\n")

    # Table S5 and Figure 5a.
    by_cell_type = io.summarise(meta, values, genes, [io.CELL_TYPE])
    by_cell_type.rename(columns={io.CELL_TYPE: "cell_type"}).to_csv(
        os.path.join(out, "by_celltype.csv"), index=False)
    wide(by_cell_type, io.CELL_TYPE, genes).rename(columns={io.CELL_TYPE: "cell_type"}) \
        .to_csv(os.path.join(out, "table_s5.csv"), index=False)
    printed_table_s5(by_cell_type, genes).to_csv(os.path.join(out, "table_s5_printed.csv"),
                                                 index=False)
    print("Detection (%) by cell type\n"
          + detection_pivot(by_cell_type, io.CELL_TYPE, genes).round(4).to_string() + "\n")

    # Figure 5b and 5c: ventricular cardiomyocytes.
    vcm = (meta[io.CELL_TYPE].astype(str) == io.VENTRICULAR_CM).to_numpy()
    vcm_meta = meta[vcm].reset_index(drop=True)
    vcm_values = {g: v[vcm] for g, v in values.items()}
    io.summarise(vcm_meta, vcm_values, genes, [io.STUDY, io.CHEMISTRY]) \
        .to_csv(os.path.join(out, "vcm_by_study_chemistry.csv"), index=False)
    chaffin = (vcm_meta[io.STUDY].astype(str) == "Chaffin").to_numpy()
    chaffin_by_disease = io.summarise(vcm_meta[chaffin].reset_index(drop=True),
                                      {g: v[chaffin] for g, v in vcm_values.items()},
                                      genes, [io.DISEASE])
    chaffin_by_disease.to_csv(os.path.join(out, "vcm_chaffin_by_disease.csv"), index=False)
    print(f"Ventricular cardiomyocytes, Chaffin et al. ({int(chaffin.sum()):,} nuclei), detection (%)\n"
          + detection_pivot(chaffin_by_disease, io.DISEASE, genes).round(4).to_string() + "\n")

    # Figure S3: all nuclei by stratum.
    for column, name in [(io.STUDY, "by_study"), (io.CHEMISTRY, "by_chemistry"),
                         (io.REGION, "by_region"), (io.DISEASE, "by_disease")]:
        io.summarise(meta, values, genes, [column]).to_csv(
            os.path.join(out, f"{name}.csv"), index=False)

    # Percentage of donors with any detection, by cell type.
    per_donor = io.summarise(meta, values, genes, [io.CELL_TYPE, io.DONOR])
    support = (per_donor.assign(any_detected=per_donor["n_detected"] > 0)
               .groupby([io.CELL_TYPE, "gene"], observed=True)
               .agg(donors_total=(io.DONOR, "nunique"),
                    donors_detecting=("any_detected", "sum"))
               .reset_index())
    support["donor_pct"] = (100 * support["donors_detecting"] / support["donors_total"]).round(1)
    support.to_csv(os.path.join(out, "donor_support.csv"), index=False)

    # Values reported in the text.
    text = text_values(meta, len(meta), by_cell_type, chaffin_by_disease)
    text.to_csv(os.path.join(out, "text_values.csv"), index=False)
    print("Values reported in the text (unrounded)\n"
          + text[["item", "value"]].to_string(index=False) + "\n")

    print(f"Tables written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()

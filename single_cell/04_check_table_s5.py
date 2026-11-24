#!/usr/bin/env python3
"""Check the regenerated Table S5, and the single-nucleus values reported in
the text, against the reported values.

    python3 04_check_table_s5.py outputs/single_cell/tables/table_s5.csv

Three comparisons, with the files written by 03_analyse.py in the folder of
the given table_s5.csv:
  table_s5.csv          every cell type and value in reference/table_s5_published.csv
                        (Table S5 to four decimal places, including nuclei with
                        the atlas label "cell", which the printed table omits)
  table_s5_printed.csv  Table S5 as printed in the Supplementary Information
                        (reference/table_s5_supplement.csv), character for character
  text_values.csv       the values reported in the main text, the Figure 5 legend and
                        the Supplementary Figure S3 legend
                        (reference/text_values_published.csv), rounded as printed
The script exits with an error if any value differs.
"""

import sys
from pathlib import Path

import pandas as pd

REFERENCE = Path(__file__).resolve().parent / "reference"
TOLERANCE = 0.00005 + 1e-9   # half a unit in the fourth decimal place
WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
         8: "eight", 9: "nine", 10: "ten"}


def check_four_decimals(path):
    published = pd.read_csv(REFERENCE / "table_s5_published.csv").set_index("cell_type")
    regenerated = pd.read_csv(path).set_index("cell_type")

    missing_rows = sorted(set(published.index) - set(regenerated.index))
    missing_cols = sorted(set(published.columns) - set(regenerated.columns))
    if missing_rows or missing_cols:
        return [f"table_s5.csv: missing cell types {missing_rows} or columns {missing_cols}"], 0

    regenerated = regenerated.loc[published.index, published.columns]
    differs = (regenerated - published).abs() > TOLERANCE
    failures = [f"table_s5.csv, {cell_type}, {column}: published {published.loc[cell_type, column]}, "
                f"regenerated {regenerated.loc[cell_type, column]}"
                for cell_type in published.index for column in published.columns
                if differs.loc[cell_type, column]]
    return failures, published.size


def check_printed(path):
    published = pd.read_csv(REFERENCE / "table_s5_supplement.csv", dtype=str,
                            keep_default_na=False).set_index("Cell type")
    if not path.exists():
        return [f"{path.name} not found; run 03_analyse.py"], 0
    regenerated = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("Cell type")
    if list(regenerated.columns) != list(published.columns):
        return [f"{path.name}: columns {list(regenerated.columns)} differ from the printed table"], 0
    if sorted(regenerated.index) != sorted(published.index):
        return [f"{path.name}: cell types {sorted(set(regenerated.index) ^ set(published.index))} "
                f"differ from the printed table"], 0
    failures = [f"{path.name}, {cell_type}, {column}: printed {published.loc[cell_type, column]}, "
                f"regenerated {regenerated.loc[cell_type, column]}"
                for cell_type in published.index for column in published.columns
                if regenerated.loc[cell_type, column] != published.loc[cell_type, column]]
    return failures, published.size


def formatted(value, style):
    if style == "grouped":
        return f"{round(value):,}".replace(",", " ")
    if style == "word":
        return WORDS.get(round(value), str(round(value)))
    if style.startswith("decimals:"):
        return f"{value:.{int(style.split(':')[1])}f}"
    raise ValueError(f"unknown format {style}")


def check_text(path):
    published = pd.read_csv(REFERENCE / "text_values_published.csv", dtype=str,
                            keep_default_na=False)
    if not path.exists():
        return [f"{path.name} not found; run 03_analyse.py"], 0
    regenerated = pd.read_csv(path).set_index("item")["value"]
    failures = []
    for row in published.itertuples(index=False):
        if row.item not in regenerated.index:
            failures.append(f"{path.name}, {row.item}: not computed")
            continue
        value = float(regenerated[row.item])
        text = formatted(value, row.format)
        if text != row.printed:
            failures.append(f"{path.name}, {row.item} ({row.where}): printed {row.printed}, "
                            f"regenerated {text} (unrounded {value:.6g})")
    return failures, len(published)


def main(path):
    path = Path(path)
    results = [("Table S5, four decimal places", check_four_decimals(path)),
               ("Table S5 as printed", check_printed(path.parent / "table_s5_printed.csv")),
               ("values in the text", check_text(path.parent / "text_values.csv"))]
    failed = False
    for name, (failures, checked) in results:
        if failures:
            failed = True
            print(f"{name}: {len(failures)} of {checked} values differ from the published values:")
            for line in failures:
                print(f"  {line}")
        else:
            print(f"{name}: all {checked} values reproduced.")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])

#!/usr/bin/env bash
# Runs all analyses and draws Figures 1-4 and Supplementary Figures S1-S2.
# Needs the study data files in data/ (not public; see data/README.md).
# The single-nucleus analysis (Figure 5, Supplementary Figure S3, Table S5)
# is run separately: bash single_cell/run_single_cell.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"

# 1. Analyses: each script writes its tables to outputs/ and stops if a result
#    differs from the value reported in the manuscript.
cd "${ROOT}/code"
for script in 01_primary_models.py 03_make_table_s1c.py \
              04_marker_comparison.py 05_spline_and_contrasts.py \
              06_diastolic_analyses.py 07_assay_batches_and_quantification.py \
              08_descriptives.py 09_pth_anatomy_and_robustness.py \
              10_lv_mass_decomposition.py; do
    echo "=== ${script} ==="
    python3 "${script}"
done

# 2. The tables as printed in the manuscript (outputs/tables/).
echo "=== 11_publication_tables.py ==="
python3 11_publication_tables.py

# 3. Figures (outputs/figures/submission/) and their source data.
cd "${ROOT}/figures"
for script in figure1.py figure2.py figureS1.py submission_figures.py export_source_data.py; do
    echo "=== ${script} ==="
    python3 "${script}"
done

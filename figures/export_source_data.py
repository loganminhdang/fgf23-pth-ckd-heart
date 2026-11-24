"""Copy the source data of Figures 1 to 4 and Supplementary Figures S1 and S2
from outputs/ to figure_source_data/ (file names listed in
figure_source_data/README.md).

Run after the analysis and figure scripts (run_all.sh does this). The files
hold summary values only (counts, histograms, estimates, fitted curves), no
participant-level data.
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS, DEST = ROOT / "outputs", ROOT / "figure_source_data"

SOURCES = {   # name in figure_source_data/ : file in outputs/
    "figure1a_flow.csv": "figure1a_flow.csv",
    "figure1b_histogram.csv": "figure1b_histogram.csv",
    "figure1b_stages.csv": "figure1b_stages.csv",
    "figure1c_histogram.csv": "figure1c_histogram.csv",
    "figure1d_geometry.csv": "figure1d_geometry.csv",
    "figure1e_ee_histogram.csv": "figure1e_ee_histogram.csv",
    "figure1e_lavi_histogram.csv": "figure1e_lavi_histogram.csv",
    "figure1f_trajectories.csv": "figure1f_trajectories.csv",
    "figure2a_fit.csv": "figure2a_fit.csv",
    "figure2b_phenotypes.csv": "figure2b_phenotypes.csv",
    "figure2c_repeated.csv": "figure2c_repeated.csv",
    "figure2d_assay_sensitivity.csv": "07_figure2d.csv",
    "figure3_lv_mass_decomposition.csv": "10_lv_mass_decomposition.csv",
    "figure3a_markers.csv": "04_table_s3a.csv",
    "figure3b_pth_vs_fgf23.csv": "04_figure3b.csv",
    "figure3c_dimensions.csv": "04_table_s3c.csv",
    "figure4ab_curves.csv": "05_figure4ab_curves.csv",
    "figure4ab_strips.csv": "05_figure4ab_strips.csv",
    "figure4c_all_observations.csv": "09_table_s4a.csv",
    "figure4c_normal_geometry.csv": "09_table_s4b.csv",
    "figure4d_diastolic.csv": "06_table_s4c.csv",
    "figure4e_pth_reference_limit.csv": "05_pth_uln_categories.csv",
    "figureS1a_batch_by_visit.csv": "07_batch_by_visit.csv",
    "figureS1b_batch_distribution.csv": "07_batch_distribution.csv",
    "figureS2a_specifications.csv": "09_table_s4d.csv",
    "figureS2b_contrasts.csv": "05_spline_contrasts.csv",
    "figureS2b_curve.csv": "05_spline_curve.csv",
    "figureS2b_pth_histogram.csv": "05_pth_histogram.csv",
}


def main():
    sources = SOURCES
    missing = [src for src in sources.values() if not (OUTPUTS / src).exists()]
    if missing:
        raise SystemExit(f"run the analysis and figure scripts first; missing {missing}")
    DEST.mkdir(exist_ok=True)
    for name, src in sources.items():
        shutil.copyfile(OUTPUTS / src, DEST / name)
    print(f"{len(sources)} source-data files written to {DEST}")


if __name__ == "__main__":
    main()

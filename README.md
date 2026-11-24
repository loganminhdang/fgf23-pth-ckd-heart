# Mineral metabolism and myocardial structure in non-dialysis CKD

Analysis code for the manuscript on PTH, FGF23 and left ventricular (LV)
structure in 200 adults with non-dialysis chronic kidney disease (CKD) and no
established heart disease, and on FGF receptor expression in a public human
heart single-nucleus atlas.

The code produces every table and figure of the manuscript and the
Supplementary Information from the study dataset, and checks each key result
against the value reported in the manuscript.

## What this repository contains

| Folder | Content |
|---|---|
| `code/` | Statistical analyses (one script per group of results) |
| `figures/` | Figures 1 to 4 and Supplementary Figures S1 and S2; Figure 5 and Supplementary Figure S3 from the single-nucleus tables |
| `single_cell/` | Analysis of the public heart atlas (Table S5, Figure 5, Supplementary Figure S3) |
| `outputs/` | Summary results written by the analyses (estimates, confidence intervals, counts); `outputs/tables/` holds each table as printed |
| `figure_source_data/` | Values plotted in Figures 1 to 4 and Supplementary Figures S1 and S2 |
| `data/` | Data dictionary only (`data/README.md`) |

**No participant data are included.** The pseudonymized study dataset is
available from the corresponding author on reasonable request, subject to
approval by the Institutional Review Board of Vinmec Times City International
Hospital and a data use agreement (see the Data Sharing Statement). Its files
are described in `data/README.md`; to run the analyses, place them in `data/`.
Result files that list individual participants (Supplementary Table S1b) are
written by the code but are not included here.

**Data preparation is not included.** The code that built the study dataset
from the clinical records (screening log, laboratory result forms,
pseudonymization) works on identifiable records. It is available, with the
dataset, to researchers who have a data use agreement. Two of its results are
used here: the participant-flow counts in `outputs/02_participant_flow.csv`
(Figure 1a) and Supplementary Table S1a (patients excluded at screening).

**The single-nucleus atlas is public.** HeartMap (Datar et al.) is available
from the Broad Institute Single Cell Portal, accession SCP3689. The code in
`single_cell/` reproduces Table S5, Figure 5 and Supplementary Figure S3 from
it; the values plotted can be obtained from the atlas directly.

## How to run

Python 3.12 or later.

```bash
pip install -r requirements.txt
bash run_all.sh                    # all analyses, tables, Figures 1-4, S1 and S2

pip install -r requirements-single-cell.txt
bash single_cell/run_single_cell.sh data/HeartMap_V1.0.h5ad   # Table S5, Figure 5, Figure S3
```

`run_all.sh` takes about a minute. The first single-cell run reads the atlas
(about 44 GB) and saves the analysed genes in a small cache
(`data/heartmap_cache/`); later runs use the cache only. After the single-cell
run, run `python3 code/11_publication_tables.py` again to write Table S5 as
printed.

## Built-in checks

Each analysis script stops with an error if a key result differs from the
value reported in the manuscript or Supplementary Tables (the `REPORTED_`
values at the top of each script). Script 01 also checks the analysis set:
200 participants, 865 visits, 524 echocardiograms, 523 passing quality
control. A successful run of `run_all.sh` therefore confirms the reported
results.

## Where each result comes from

| Result | Script | Output |
|---|---|---|
| Table 1 | `code/08_descriptives.py` | `outputs/tables/table_1.csv` |
| Figure 1 | `figures/figure1.py`, `figures/submission_figures.py` | `outputs/figures/submission/Figure 1.tif` |
| Figure 2 | `code/01_primary_models.py`, `code/07_assay_batches_and_quantification.py`, `figures/figure2.py`, `figures/submission_figures.py` | `outputs/figures/submission/Figure 2.tif` |
| Figure 3 | `code/04_marker_comparison.py`, `code/10_lv_mass_decomposition.py`, `figures/submission_figures.py` | `outputs/figures/submission/Figure 3.tif` |
| Figure 4 | `code/05_spline_and_contrasts.py`, `code/06_diastolic_analyses.py`, `code/09_pth_anatomy_and_robustness.py`, `figures/submission_figures.py` | `outputs/figures/submission/Figure 4.tif` |
| Figure 5 | `single_cell/run_single_cell.sh`, `figures/figure5.py` | `outputs/figures/submission/Figure 5.tif` |
| Supplementary Table S1a | data preparation (not public) | |
| Supplementary Table S1b | `code/08_descriptives.py` | `outputs/tables/table_s1b.csv` (not included) |
| Supplementary Table S1c | `code/03_make_table_s1c.py` | `outputs/tables/table_s1c.csv` |
| Supplementary Table S2 | `code/07_assay_batches_and_quantification.py` | `outputs/tables/table_s2.csv` |
| Supplementary Tables S3a–S3d | `code/04_marker_comparison.py` | `outputs/tables/table_s3a.csv` to `table_s3d.csv` |
| Supplementary Tables S4a, S4b, S4d | `code/09_pth_anatomy_and_robustness.py` | `outputs/tables/table_s4a.csv`, `table_s4b.csv`, `table_s4d.csv` |
| Supplementary Table S4c | `code/06_diastolic_analyses.py` | `outputs/tables/table_s4c.csv` |
| Supplementary Table S4e | `code/05_spline_and_contrasts.py` | `outputs/tables/table_s4e.csv` |
| Supplementary Table S5 | `single_cell/03_analyse.py` | `outputs/tables/table_s5.csv` |
| Supplementary Figure S1 | `code/07_assay_batches_and_quantification.py`, `figures/figureS1.py`, `figures/submission_figures.py` | `outputs/figures/submission/Supplementary Figure S1.tif` |
| Supplementary Figure S2 | `code/05_spline_and_contrasts.py`, `code/09_pth_anatomy_and_robustness.py`, `figures/submission_figures.py` | `outputs/figures/submission/Supplementary Figure S2.tif` |
| Supplementary Figure S3 | `single_cell/run_single_cell.sh`, `figures/figureS3.py` | `outputs/figures/submission/Supplementary Figure S3.tif` |

Numbers quoted in the text but not shown in a table are written to
`outputs/` by the same scripts, for example `01_primary_models.csv`,
`07_linkage_summary.csv`, `08_lvmi_by_visit.csv` and `09_kidney_function.csv`.
Figures are also saved as PNG, PDF and SVG (editable text).

## Analysis scripts (`code/`)

| Script | What it does |
|---|---|
| `prep.py` | Shared steps: reads the dataset, applies echocardiographic quality control, derives LVH and LV geometry |
| `01_primary_models.py` | Primary PTH and FGF23 models; PTH at each visit; within- and between-patient PTH associations |
| `03_make_table_s1c.py` | Characteristics at each study visit (Table S1c) |
| `04_marker_comparison.py` | Each mineral marker against LV mass and LV dimensions; PTH versus FGF23 (Table S3, Figure 3) |
| `05_spline_and_contrasts.py` | Shape of the PTH–LV mass association; PTH relative to the upper reference limit (Figure 4a, b, e; Table S4e; Figure S2b) |
| `06_diastolic_analyses.py` | Diastolic function (Table S4c, Figure 4d) |
| `07_assay_batches_and_quantification.py` | FGF23 assay batches and limits of quantification (Table S2, Figure 2d, Figure S1) |
| `08_descriptives.py` | Baseline characteristics (Table 1), LVH at baseline (Table S1b), LV mass index over follow-up |
| `09_pth_anatomy_and_robustness.py` | PTH and LV dimensions; robustness of the PTH association; kidney function (Tables S4a, S4b, S4d; Figure 4c; Figure S2a) |
| `10_lv_mass_decomposition.py` | Contributions of wall thickness and cavity size to the LV mass associations (Figure 3 legend) |
| `11_publication_tables.py` | Every table as printed, rounded from the unrounded results in `outputs/unrounded/` |

Methods in brief: each echocardiogram is one observation. LV measures are
related to each marker by linear regression adjusted for age, sex, BMI,
systolic and diastolic blood pressure, heart rate and eGFR, with standard
errors that allow for repeated echocardiograms of the same patient. PTH and
FGF23 are analysed on the log10 scale, so estimates are per 10-fold higher
concentration unless given per standard deviation (SD).

## Figure scripts (`figures/`)

| Script | What it does |
|---|---|
| `figure1.py` | Panels of Figure 1 and their source data |
| `figure2.py` | Panels of Figure 2 and their source data |
| `figureS1.py` | Panels of Supplementary Figure S1 |
| `submission_figures.py` | Figures 1 to 4 and Supplementary Figures S1 and S2 as submitted |
| `figure5.py`, `figureS3.py` | Figure 5 and Supplementary Figure S3 from the single-nucleus tables |
| `export_source_data.py` | Copies the plotted values to `figure_source_data/` |
| `kipanel.py`, `style.py`, `vector_svg.py` | Shared figure settings; fully editable SVG copies of Figure 1 and Supplementary Figure S1 |

## Single-cell scripts (`single_cell/`)

| Script | What it does |
|---|---|
| `run_single_cell.sh` | Runs the steps below in order |
| `01_inspect.py` | Lists the structure of the atlas file |
| `02_extract.py` | Saves the analysed genes and the nucleus annotations in a small cache |
| `03_analyse.py` | Percentage of nuclei with each transcript, and mean expression, by cell type, study, chemistry, region and disease state (Table S5, Figure 5, Figure S3, values in the text) |
| `04_check_table_s5.py` | Checks Table S5 and the values in the text against the reported values in `reference/` |
| `atlas_io.py` | Reading the atlas and the cache |

## Software

Package versions are fixed in `requirements.txt` (results were produced with
Python 3.13). The bootstrap in `code/04_marker_comparison.py` uses a fixed
random seed (20261021), so results are identical on every run. Figures are
sized for Helvetica; where it is not installed, Arial or Liberation Sans is
used and the figures differ slightly in appearance.

## License

MIT (see `LICENSE`).

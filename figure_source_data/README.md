# Figure source data

The values plotted in Figures 1 to 4 and Supplementary Figures S1 and S2,
written by the scripts in `code/` and `figures/` and copied here by
`figures/export_source_data.py`. The files hold summary values only (counts,
histogram bins, estimates with 95% confidence intervals, fitted curves); no
participant-level data. The points in Figure 2a are individual
echocardiograms and are therefore not included.

| Files | Figure |
|---|---|
| `figure1a_flow.csv`, `figure1b_histogram.csv`, `figure1b_stages.csv`, `figure1c_histogram.csv`, `figure1d_geometry.csv`, `figure1e_ee_histogram.csv`, `figure1e_lavi_histogram.csv`, `figure1f_trajectories.csv` | Figure 1 |
| `figure2a_fit.csv`, `figure2b_phenotypes.csv`, `figure2c_repeated.csv`, `figure2d_assay_sensitivity.csv` | Figure 2 |
| `figure3a_markers.csv`, `figure3b_pth_vs_fgf23.csv`, `figure3c_dimensions.csv` | Figure 3 (the eGFR row of panel a and the relative wall thickness rows of panel c are not drawn) |
| `figure3_lv_mass_decomposition.csv` | Figure 3 legend (wall and cavity contributions) |
| `figure4ab_curves.csv`, `figure4ab_strips.csv` | Figure 4a, b |
| `figure4c_all_observations.csv`, `figure4c_normal_geometry.csv` | Figure 4c |
| `figure4d_diastolic.csv` | Figure 4d (rows with Exposure = PTH and Covariates = primary) |
| `figure4e_pth_reference_limit.csv` | Figure 4e |
| `figureS1a_batch_by_visit.csv`, `figureS1b_batch_distribution.csv` | Supplementary Figure S1 |
| `figureS2a_specifications.csv` | Supplementary Figure S2a (the height^2.7 row is not drawn) |
| `figureS2b_curve.csv`, `figureS2b_contrasts.csv`, `figureS2b_pth_histogram.csv` | Supplementary Figure S2b |

Figure 5 and Supplementary Figure S3 are drawn from the public HeartMap atlas
(Broad Institute Single Cell Portal, accession SCP3689); `single_cell/`
reproduces them from it.

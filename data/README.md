# Data dictionary

The study data are not public. The pseudonymized dataset is available from the
corresponding author on reasonable request (see the Data Sharing Statement of
the manuscript). To run the analyses, place the three files below in this
folder. Participants are identified only by a study pseudonym.

## visits.csv

One row per participant visit: 200 participants, visits 1 to 5, 865 rows.
Echocardiography was performed at visits 1, 3 and 5.

| Column | Content |
|---|---|
| `Patient_ID` | participant pseudonym |
| `Visit` | study visit, 1 to 5 |
| `Age` | years, at the visit |
| `Sex` | Male or Female |
| `Height`, `Weight` | cm, kg |
| `BMI` | kg/m² |
| `SBP`, `DBP` | systolic and diastolic blood pressure, mm Hg |
| `Heart_rate` | beats per minute |
| `Creatinine` | µmol/l |
| `Calcium`, `Phosphate` | mmol/l |
| `Albumin` | g/l |
| `eGFR` | ml/min/1.73 m², CKD-EPI 2021 creatinine equation (without race) |
| `PTH`, `FGF23` | pg/ml. FGF23 values below the limit of quantification hold the concentration extrapolated from the standard curve, corrected for dilution |
| `LVIDd`, `LVIDs` | LV internal dimension in diastole and systole, mm |
| `IVSd`, `IVSs` | septal thickness in diastole and systole, mm |
| `PWd_used` | posterior wall thickness in diastole, mm |
| `LVMI` | LV mass index, g/m²: Devereux LV mass divided by Du Bois body surface area. Blank where the echocardiogram fails quality control |
| `RWT` | relative wall thickness, 2 × PWd_used / LVIDd |
| `E`, `A` | mitral E and A velocities, cm/s |
| `EA_reported` | E/A as reported (the analyses recalculate E/A from E and A) |
| `e_septal`, `e_lateral` | septal and lateral e′, cm/s |
| `E_e_avg` | average E/e′ |
| `TRV` | tricuspid regurgitation velocity, m/s |
| `LAVI` | left atrial volume index, ml/m² |
| `GLS` | LV global longitudinal strain, % |

## lab_linked_v2.csv

The FGF23 assay batch of each measurement: one row per FGF23 measurement (865).

| Column | Content |
|---|---|
| `pid` | participant pseudonym, as `Patient_ID` in `visits.csv` |
| `tp` | study visit, 1 to 5 |
| `batch` | assay batch, B01 to B15 in the order assayed |
| `plate_year` | year in which the sample was assayed |
| `curve_lo`, `curve_hi` | range of the standard curve, pg/ml |
| `factor` | dilution factor, where one was recorded |
| `limit` | "below" for values reported below the limit of quantification, otherwise blank |
| `bound` | for those values, the limit of quantification that applied, pg/ml (78, 39 or 20) |

## study_workbook.xlsx

Primary kidney disease and year of CKD diagnosis of the 200 participants (one
sheet; column labels in row 2, in Vietnamese, with English names in row 1).
Read by `code/08_descriptives.py` for Table 1.

## HeartMap_V1.0.h5ad (public)

The human heart single-nucleus atlas of Datar et al. (about 44 GB), from the
Broad Institute Single Cell Portal, accession SCP3689. Used by the scripts in
`single_cell/`, which save a small cache in `heartmap_cache/` in this folder.

# Interpretable Decision Tools for Multi-Timepoint Mortality Risk Stratification After Acute Myocardial Infarction

R code corresponding to the analyses reported in the manuscript *Interpretable Decision Tools for Multi-Timepoint Mortality Risk Stratification After Acute Myocardial Infarction Using Routinely Collected Admission Data*.

The repository covers data preparation, AIC-based logistic and Cox regression models, repeated five-fold cross-validation, nomogram construction, held-out test-set evaluation, ROSE sensitivity analyses, and k-prototypes clustering.

## Code

Run the scripts in numerical order from the repository root.

- `code/01_data_preparation.R` — complete-case preparation and 70:30 training/held-out test partitioning
- `code/02_model_development.R` — ROSE training data, AIC selection, repeated five-fold cross-validation (20 repeats), logistic/Cox models, and nomograms
- `code/03_model_evaluation.R` — held-out test-set evaluation reported in the manuscript
- `code/04_sensitivity_analysis.R` — ROSE sensitivity analysis corresponding to Table S2 and Fig. S1
- `code/05_clustering.R` — k-prototypes clustering corresponding to the reported cluster analyses

## Data

Patient-level registry data are not publicly distributed because of institutional and patient-confidentiality restrictions. The scripts expect authorised local copies at:

- `data/ami_30d.csv`
- `data/ami_longterm.csv`

Standardised analysis variables are:

`age`, `sex`, `residential_area`, `treatment_within_12h`, `onset_season`, `hypertension`, `diabetes`, `dyslipidemia`, `bmi`, `stemi`, `multivessel_disease`, `mortality_30d`, `followup_1y_days`, `event_1y`, `followup_5y_days`, and `event_5y`.

Coding follows the manuscript and supplementary information: sex (1 male, 2 female); residential area (0 rural, 1 urban); treatment within 12 hours (0 ≥12 hours, 1 <12 hours); onset season (1 spring, 2 summer, 3 autumn, 4 winter); hypertension, diabetes, dyslipidemia, STEMI, and multivessel disease (0 no/reference, 1 yes/exposed as defined in the manuscript).

## Software

R 4.4.2. Packages used in the analysis scripts: `caret`, `MASS`, `ROSE`, `rms`, `survival`, `riskRegression`, `pROC`, `dcurves`, `clustMixType`, `dplyr`, `readr`, `tidyr`, and `pheatmap`.

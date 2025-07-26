# AMI-mortality-risk-stratification-AI

This repository contains code and documentation for multi-timepoint mortality risk stratification in patients with acute myocardial infarction (AMI), developed using supervised regression and unsupervised clustering approaches.

## Objectives

* To construct an interpretable, modular pipeline for stratifying AMI mortality risk at distinct time points (30-day, 1-year, 5-year)
* To apply supervised modeling strategies (logistic regression for short-term; Cox regression for long-term) to support clinical decision-making
* To identify risk-based patient subtypes using unsupervised clustering derived from regression-based profiles

---

## Project Structure

```
AMI-mortality-risk-stratification-AI/
│
├── data/
│   └── dummy_input_30days.csv          # Simulated input (no real patient data)
│
├── code/
│   ├── short_term/                     # 30-day mortality (logistic regression)
│   │   ├── 1_data_preparation_short.R
│   │   ├── 2_ROSE_balancing_short.R
│   │   ├── 3_model_logistic_nomogram.R
│   │   └── 4_evaluation_logistic.R
│   │
│   └── long_term/                      # 1-year and 5-year mortality (Cox regression)
│       ├── 1_data_preparation_long.R
│       ├── 2_ROSE_balancing_long.R
│       ├── 3_model_cox_nomogram.R
│       └── 4_evaluation_cox.R
│
├── clustering/                        # Unsupervised subtyping
│   └── 5_clustering_kproto.R
│
├── README.md
└── LICENSE
```

---

## Software Requirements

This project requires **R (>= 4.2.0)** and the following CRAN packages:

* `haven`, `ROSE`, `rms`, `caret`, `pROC`, `survival`, `Hmisc`
* `clustMixType`, `cluster`, `dplyr`, `ggplot2`, `factoextra`

To install dependencies:

```r
install.packages(c("haven", "ROSE", "rms", "caret", "pROC", "survival", "Hmisc", 
                   "clustMixType", "cluster", "dplyr", "ggplot2", "factoextra"))
```

---

## Input Format

The pipeline requires a `.csv` file with the following structure:

```csv
Sex,Age_C,Time,Residential_area,HTN,DM,DYS,OB_C_bmi,STEMI,Bloodvessel_C,FU_30days,event_status_30
1,67,1,0,1,1,0,23.4,1,1,30,0
...
```

A dummy template (`data/dummy_input_30days.csv`) is included to illustrate the required input schema. This file contains column names only, with no patient data.

---

## Methodological Overview

* **Short-term model**: Logistic regression with ROSE balancing and nomogram visualization for 30-day mortality
* **Long-term model**: Cox proportional hazards regression applied to 1-year and 5-year survival estimation, sharing a unified pipeline
* **Unsupervised clustering**: k-prototypes clustering on variables selected via AIC to reveal patient subgroups

---

## Reproducibility Practices

This repository applies key reproducibility principles:

* Fixed random seeds (`set.seed(123)`)
* Stepwise and annotated scripts
* Open-source R packages only
* Dummy input schema provided for structural reference

These measures support consistent, transparent, and ethically compliant analysis workflows.

---

## Data

The dataset used in this study is not publicly available due to patient privacy protection and institutional regulations. The repository includes a dummy input file (`data/dummy_input_30days.csv`) that mimics the required structure without exposing any individual-level data.

---

## License

This project is licensed under the MIT License. See `LICENSE` for terms.

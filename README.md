# AMI-mortality-risk-stratification-AI

This repository contains code and documentation for a two-stage mortality risk stratification approach in patients with acute myocardial infarction (AMI), integrating supervised regression and unsupervised clustering approaches.

## Objectives

* To construct an interpretable，machine learning-based for stratifying AMI mortality risk at distinct time points (30-day, 1-year, 5-year)
* To apply supervised modeling strategies (logistic regression for short-term; Cox regression for long-term) to support clinical decision-making
* To identify risk-based patient subtypes using unsupervised clustering derived from regression-based profiles

---

## Project Structure

```
AMI-mortality-risk-stratification-AI/
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

The pipeline requires a `.csv` file with the following columns:

```csv
Sex,Age_C,Time,Residential_area,HTN,DM,DYS,OB_C_bmi,STEMI,Bloodvessel_C,FU_30days,event_status_30
```

Each row corresponds to a single patient record. Below is a sample row illustrating the expected format:

```csv
1,67,1,0,1,1,0,23.4,1,1,30,0
```

### Column Definitions:

| Column             | Description                                            |
| ------------------ | ------------------------------------------------------ |
| `Sex`              | 1 = Male, 2 = Female                                   |
| `Age_C`            | Age in years (continuous)                              |
| `Time`             | 0 = ≥12 hours, 1 = <12 hours from symptom to admission |
| `Residential_area` | 0 = Rural, 1 = Urban                                   |
| `HTN`              | 0 = No hypertension, 1 = Yes                           |
| `DM`               | 0 = No diabetes, 1 = Yes                               |
| `DYS`              | 0 = No dyslipidemia, 1 = Yes                           |
| `OB_C_bmi`         | Body Mass Index (BMI, numeric)                         |
| `STEMI`            | 0 = NSTEMI, 1 = STEMI                                  |
| `Bloodvessel_C`    | 0 = Single/none, 1 = Multi-vessel disease              |
| `FU_30days`        | Follow-up time in days (e.g., 30)                      |
| `event_status_30`  | 0 = Alive, 1 = Deceased within 30 days                 |

## Methodological Overview

* **Short-term model**: Logistic regression with ROSE balancing and nomogram visualization for 30-day mortality
* **Long-term model**: Cox proportional hazards regression applied to 1-year and 5-year survival estimation
* **Unsupervised clustering**: k-prototypes clustering on variables selected via AIC to reveal patient subgroups

---

## Data

The dataset used in this study is not publicly available due to patient privacy protection and institutional regulations. 

---

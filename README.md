# AMI-mortality-risk-stratification-AI
Code and documentation for AMI mortality risk stratification via nomogram and clustering.
## 📌 Objectives

* To develop an interpretable, reproducible pipeline for multi-timepoint mortality risk stratification following acute myocardial infarction (AMI)
* To apply supervised regression-based modeling (logistic for 30-day; Cox for 1-year) for clinically meaningful prediction
* To identify subgroups of AMI patients through unsupervised clustering based on nomogram-derived risk profiles

---

## 📁 Project Structure

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
│   └── long_term/                      # 1-year mortality (Cox regression)
│       ├── 1_data_preparation_long.R
│       ├── 2_ROSE_balancing_long.R
│       ├── 3_model_cox_nomogram.R
│       └── 4_evaluation_cox.R
│
├── clustering/                        # Unsupervised subtyping
│   └── 5_clustering_kproto.R
│
├── results/                           # Example visual outputs
│   ├── example_nomogram_30d.png
│   ├── example_nomogram_1y.png
│   └── example_heatmap_clusters.png
│
├── README.md
└── LICENSE
```

---

## 🔧 Dependencies

This project uses **R (>= 4.2.0)** and the following packages:

* `haven`, `ROSE`, `rms`, `caret`, `pROC`, `survival`, `Hmisc`
* `clustMixType`, `cluster`, `dplyr`, `ggplot2`, `factoextra`

To install all required packages:

```r
install.packages(c("haven", "ROSE", "rms", "caret", "pROC", "survival", "Hmisc", 
                   "clustMixType", "cluster", "dplyr", "ggplot2", "factoextra"))
```

---

## 📊 Input Format

The pipeline expects a `.csv` file with the following structure:

```csv
Sex,Age_C,Time,Residential_area,HTN,DM,DYS,OB_C_bmi,STEMI,Bloodvessel_C,FU_30days,event_status_30
1,67,1,0,1,1,0,23.4,1,1,30,0
...
```

> 📁 `data/dummy_input_30days.csv` is provided as a template with column names only.

---

## 🧪 Pipeline Summary

* **Short-term (30-day)**: Logistic regression + nomogram (internal validation, AUC)
* **Long-term (1-year)**: Cox regression + nomogram (C-index, calibration, ROC)
* **Clustering**: Unsupervised k-prototypes analysis based on selected features

---

## 🔁 Reproducibility

This repository follows best practices promoted by international reproducibility initiatives (e.g., **FAIR**, **TOP Guidelines**, **NIH Rigor and Reproducibility Framework**):

* All random operations use fixed seeds (`set.seed(123)`)
* Modular, stepwise code with filenames reflecting their function and execution order
* Transparent preprocessing, balancing, model training, validation, and evaluation steps
* Separation of internal (CV, C-index) and external (ROC, AUC, optimal threshold) validation
* Dummy input template provided to emulate workflow structure without compromising privacy
* No proprietary or institution-specific functions are used
* Code is compatible with open-source R ecosystem and can be executed on any academic machine

These practices ensure reproducibility, transparency, and extensibility for peer reviewers, collaborators, and the broader research community.

---

## 📄 Citation

If you use this repository or adapt its code, please cite the corresponding manuscript:

> Hsieh TC, et al. Two-stage machine learning for mortality risk stratification in acute myocardial infarction via complementary nomogram and clustering. *Submitted, 2025.*

---

## 📜 License

This project is released under the MIT License. See [LICENSE](./LICENSE) for details.

---

## 🙋 Contact

For questions about the code or requests for collaboration:

📧 [tchsieh@gms.tcu.edu.tw](mailto:tchsieh@gms.tcu.edu.tw)
👤 Corresponding Author: Tsung-Cheng Hsieh

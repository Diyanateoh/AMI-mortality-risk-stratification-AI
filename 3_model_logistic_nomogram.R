## Short-term Mortality: Logistic Regression & Nomogram (30-day)

```r
# ------------------------------------------------------------------------------
# Module: 3_model_logistic_nomogram.R
# Purpose: Construct logistic regression model and generate nomogram
# Input: ROSE-balanced dataset
# Output: Fitted model and visual nomogram object
# ------------------------------------------------------------------------------
# Description:
# This script fits a logistic regression model to predict 30-day mortality following
# acute myocardial infarction (AMI) using the balanced dataset from Module 2.
# The model includes variables selected via stepwise AIC criteria and is visualized
# as a nomogram using the rms package. This model supports interpretability and risk communication.
#
# Performance metrics such as AUC and threshold-based classification indices
# are computed in Module 4.
# ------------------------------------------------------------------------------
```

```r
library(rms)
library(Hmisc)

# Ensure data is available
if (!exists("rose_data")) stop("Balanced dataset 'rose_data' not found")

# Specify data distribution object for rms
dd <- datadist(rose_data)
options(datadist = "dd")

# Define model formula with AIC-selected variables
logit_model <- lrm(event_status_30 ~ Age_C + Time + Residential_area + DM + HTN + STEMI,
                   data = rose_data, x = TRUE, y = TRUE)

# Generate nomogram
nomo_short <- nomogram(logit_model,
                       fun = plogis,
                       funlabel = "30-day Mortality Probability",
                       lp = FALSE,
                       Age_C = seq(20, 95, by = 5),
                       fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))

# Plot nomogram
plot(nomo_short, xfrac = 0.3, cex.axis = 0.8, cex.var = 1, col.grid = gray(0.8),
     main = "Nomogram for 30-day Mortality Prediction")

# Output model summary
summary(logit_model)
```

# ------------------------------------------------------------------------------
# Module: 3_model_cox_nomogram.R
# Purpose: Fit Cox regression models for 1-year and 5-year mortality and create nomograms
# Input: Balanced datasets from Module 2 (ROSE output)
# Output: Cox models and nomogram objects for 1-year and 5-year mortality prediction
# ------------------------------------------------------------------------------

# Description:
# This script fits two Cox proportional hazards models to predict long-term mortality
# at 1-year and 5-year intervals following acute myocardial infarction (AMI).
# Predictors are selected using AIC-based variable selection (as per manuscript).
# Nomograms are generated using the rms package.

library(survival)
library(rms)

# Set datadist for nomogram plotting
options(datadist = 'dd')

# --- 1-Year Mortality Cox Model ---

# Assume balanced_1y is available (from Module 2)
dd <- datadist(balanced_1y)
options(datadist = 'dd')

cox_1y <- cph(Surv(FU_1y, event_status_1y) ~ Age_C + Sex + Time + Residential_area +
                HTN + DM + DYS + OB_C_bmi + STEMI + Bloodvessel_C,
              data = balanced_1y, x = TRUE, y = TRUE, surv = TRUE)

# Generate nomogram for 1-year mortality
surv_1y <- Survival(cox_1y)
nomogram_1y <- nomogram(cox_1y, fun = list(function(x) surv_1y(365, x)),
                        funlabel = "1-Year Survival Probability")

plot(nomogram_1y)

# --- 5-Year Mortality Cox Model ---

# Assume balanced_5y is available (from Module 2)
dd <- datadist(balanced_5y)
options(datadist = 'dd')

cox_5y <- cph(Surv(FU_5y, event_status_5y) ~ Age_C + Sex + Time + Residential_area +
                HTN + DM + DYS + OB_C_bmi + STEMI + Bloodvessel_C,
              data = balanced_5y, x = TRUE, y = TRUE, surv = TRUE)

# Generate nomogram for 5-year mortality
surv_5y <- Survival(cox_5y)
nomogram_5y <- nomogram(cox_5y, fun = list(function(x) surv_5y(1825, x)),
                        funlabel = "5-Year Survival Probability")

plot(nomogram_5y)

# The fitted models and nomograms are now ready for evaluation (Module 4).

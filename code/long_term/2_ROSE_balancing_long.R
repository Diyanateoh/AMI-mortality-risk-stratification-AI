# ------------------------------------------------------------------------------
# Module: 2_ROSE_balancing_long.R
# Purpose: Address class imbalance for long-term (1-year and 5-year) mortality prediction
# Input: Preprocessed dataset from Module 1
# Output: ROSE-balanced datasets for both 1-year and 5-year survival models
# ------------------------------------------------------------------------------

# Description:
# This script applies the ROSE (Random OverSampling Examples) technique to mitigate 
# outcome imbalance in long-term mortality datasets. Observed event rates were 7.7% 
# for 1-year and 19.0% for 5-year mortality. Balancing is performed separately for each timepoint.

library(ROSE)
library(dplyr)

# Assume preprocessed data from Module 1 is available as `data`

# --- 1-Year Mortality Balancing ------------------------------------------------
# Create a balanced version of the 1-year survival dataset

data_1y <- data %>%
  select(Sex, Age_C, Time, Residential_area, HTN, DM, DYS, OB_C_bmi, STEMI, Bloodvessel_C,
         FU_1y, event_status_1y)

set.seed(123)
balanced_1y <- ROSE(event_status_1y ~ ., data = data_1y, seed = 123)$data

# --- 5-Year Mortality Balancing ------------------------------------------------
# Create a balanced version of the 5-year survival dataset

data_5y <- data %>%
  select(Sex, Age_C, Time, Residential_area, HTN, DM, DYS, OB_C_bmi, STEMI, Bloodvessel_C,
         FU_5y, event_status_5y)

set.seed(123)
balanced_5y <- ROSE(event_status_5y ~ ., data = data_5y, seed = 123)$data

# --- Summary Outputs -----------------------------------------------------------

cat("\nOriginal 1-year event distribution:\n")
print(table(data$event_status_1y))

cat("\nBalanced 1-year event distribution:\n")
print(table(balanced_1y$event_status_1y))

cat("\nOriginal 5-year event distribution:\n")
print(table(data$event_status_5y))

cat("\nBalanced 5-year event distribution:\n")
print(table(balanced_5y$event_status_5y))

# The balanced datasets are now ready for Cox modeling (Module 3).

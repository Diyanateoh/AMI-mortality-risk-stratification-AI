# Short-term AMI: Data Preparation (30-day)

```r
# ------------------------------------------------------------------------------
# Module: 1_data_preparation_short.R
# Purpose: Preprocessing for 30-day mortality prediction following AMI
# Output: Cleaned dataset with harmonized variable types and structure
# ------------------------------------------------------------------------------
# Description:
# This script performs initial preprocessing of acute myocardial infarction (AMI)
# patient data for short-term (30-day) mortality analysis.
#
# Key transformations:
# - Exclusion of ID and non-modeled behavioral or seasonal variables
# - Labeling of categorical variables as factors with manuscript-consistent levels
# - Conversion of survival time (FU_30days) and binary event indicator (event_status_30)
# - Verification of missing values and variable structure post-transformation
#
# The processed dataset is ready for:
# - Class balancing via ROSE (2_ROSE_balancing_short.R)
# - Logistic regression modeling and nomogram construction (3_model_logistic_nomogram.R)

```r
library(haven)
library(dplyr)

# Load SPSS dataset (.sav format)
# Note: Replace with actual file path

data <- read_sav("your_data_path.sav")

# Remove ID column if present
if ("ID" %in% names(data)) {
  data$ID <- NULL
}

# Recode categorical variables to match manuscript definitions

data$Sex <- factor(data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
data$Time <- factor(data$Time, levels = c(0, 1), labels = c("Greater than 12 hours", "Within 12 hours"))
data$Residential_area <- factor(data$Residential_area, levels = c(1, 0), labels = c("Urban", "Rural"))
data$HTN <- factor(data$HTN, levels = c(0, 1), labels = c("No", "Yes"))
data$DM <- factor(data$DM, levels = c(0, 1), labels = c("No", "Yes"))
data$DYS <- factor(data$DYS, levels = c(0, 1), labels = c("No", "Yes"))
data$STEMI <- factor(data$STEMI, levels = c(0, 1), labels = c("NSTEMI", "STEMI"))
data$Bloodvessel_C <- factor(data$Bloodvessel_C, levels = c(0, 1), labels = c("Non-multi", "Multi diseased vessels"))

# Confirm presence and numeric conversion of continuous variables
if (!"Age_C" %in% names(data)) stop("Missing variable: Age_C")
if (!"OB_C_bmi" %in% names(data)) stop("Missing variable: OB_C_bmi")

data$Age_C <- as.numeric(data$Age_C)
data$OB_C_bmi <- as.numeric(data$OB_C_bmi)

# Ensure survival and outcome indicators are present
if (!all(c("FU_30days", "event_status_30") %in% names(data))) {
  stop("Missing variables: FU_30days and/or event_status_30")
}

data$FU_30days <- as.numeric(data$FU_30days)
data$event_status_30 <- as.numeric(data$event_status_30)

# Summary and missing value diagnostics
cat("\nVariable summary after preprocessing:\n")
print(summary(data))

cat("\nMissing values by variable:\n")
print(colSums(is.na(data)))

# The resulting dataset is ready for downstream modeling
```

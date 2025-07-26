# ------------------------------------------------------------------------------
# Module: 1_data_preparation_long.R
# Purpose: Data preprocessing for long-term (1-year and 5-year) mortality prediction
# Input: Raw dataset in SPSS (.sav) format
# Output: Cleaned and transformed dataset for modeling
# ------------------------------------------------------------------------------

# Description:
# This script performs data preprocessing for long-term survival models
# used in predicting 1-year and 5-year mortality following acute myocardial infarction (AMI).
# The script includes SPSS file loading, variable recoding, factor conversion, and missing value checks.

# Required libraries
library(haven)
library(dplyr)

# Load data (replace with actual file path)
data <- read_sav("your_data_path.sav")

# Remove patient ID if present
if ("ID" %in% names(data)) {
  data$ID <- NULL
}

# Factor conversion for categorical variables
# Based on variables used in AIC-based selection for final manuscript model

data$Sex <- factor(data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
data$Time <- factor(data$Time, levels = c(0, 1), labels = c(">12h", "<=12h"))
data$Residential_area <- factor(data$Residential_area, levels = c(1, 0), labels = c("Urban", "Rural"))
data$HTN <- factor(data$HTN, levels = c(0, 1), labels = c("No", "Yes"))
data$DM <- factor(data$DM, levels = c(0, 1), labels = c("No", "Yes"))
data$DYS <- factor(data$DYS, levels = c(0, 1), labels = c("No", "Yes"))
data$STEMI <- factor(data$STEMI, levels = c(0, 1), labels = c("NSTEMI", "STEMI"))
data$Bloodvessel_C <- factor(data$Bloodvessel_C, levels = c(0, 1), labels = c("Non-multi", "Multi diseased vessels"))

# Ensure numeric format for continuous variables
if (!"Age_C" %in% names(data)) stop("Missing variable: Age_C")
if (!"OB_C_bmi" %in% names(data)) stop("Missing variable: OB_C_bmi")

data$Age_C <- as.numeric(data$Age_C)
data$OB_C_bmi <- as.numeric(data$OB_C_bmi)

# Verify event and time-to-event variables are present
if (!all(c("FU_1y", "event_status_1y", "FU_5y", "event_status_5y") %in% names(data))) {
  stop("Missing survival time or status variables for 1-year or 5-year")
}

# Ensure numeric conversion
data$FU_1y <- as.numeric(data$FU_1y)
data$event_status_1y <- as.numeric(data$event_status_1y)
data$FU_5y <- as.numeric(data$FU_5y)
data$event_status_5y <- as.numeric(data$event_status_5y)

# Summary of processed data
cat("\nSummary statistics of processed long-term dataset:\n")
print(summary(data))

cat("\nMissing values by variable:\n")
print(colSums(is.na(data)))

# The preprocessed dataset is ready for class balancing and survival model development

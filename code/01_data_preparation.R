library(dplyr)
library(readr)
library(tidyr)

analysis_seed <- 123L

short <- read_csv("data/ami_30d.csv", show_col_types = FALSE)
long <- read_csv("data/ami_longterm.csv", show_col_types = FALSE)

short_required <- c(
  "age", "sex", "residential_area", "treatment_within_12h", "onset_season",
  "hypertension", "diabetes", "dyslipidemia", "bmi", "stemi",
  "multivessel_disease", "mortality_30d"
)

long_required <- c(
  "age", "sex", "residential_area", "treatment_within_12h",
  "hypertension", "diabetes", "dyslipidemia", "bmi", "stemi",
  "multivessel_disease", "followup_1y_days", "event_1y",
  "followup_5y_days", "event_5y"
)

stopifnot(all(short_required %in% names(short)))
stopifnot(all(long_required %in% names(long)))

encode_common <- function(x) {
  x %>%
    mutate(
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      sex = factor(as.integer(sex), levels = c(1, 2), labels = c("Male", "Female")),
      residential_area = factor(
        as.integer(residential_area),
        levels = c(0, 1),
        labels = c("Rural", "Urban")
      ),
      treatment_within_12h = factor(
        as.integer(treatment_within_12h),
        levels = c(0, 1),
        labels = c(">=12 hours", "<12 hours")
      ),
      hypertension = factor(as.integer(hypertension), levels = c(0, 1), labels = c("No", "Yes")),
      diabetes = factor(as.integer(diabetes), levels = c(0, 1), labels = c("No", "Yes")),
      dyslipidemia = factor(as.integer(dyslipidemia), levels = c(0, 1), labels = c("No", "Yes")),
      stemi = factor(as.integer(stemi), levels = c(0, 1), labels = c("NSTEMI", "STEMI")),
      multivessel_disease = factor(
        as.integer(multivessel_disease),
        levels = c(0, 1),
        labels = c("Non-multi", "Multi")
      )
    )
}

short <- short %>%
  select(all_of(short_required)) %>%
  drop_na() %>%
  encode_common() %>%
  mutate(
    onset_season = factor(
      as.integer(onset_season),
      levels = 1:4,
      labels = c("Spring", "Summer", "Autumn", "Winter")
    ),
    mortality_30d = as.integer(mortality_30d)
  )

long <- long %>%
  select(all_of(long_required)) %>%
  drop_na() %>%
  encode_common() %>%
  mutate(
    followup_1y_days = as.numeric(followup_1y_days),
    event_1y = as.integer(event_1y),
    followup_5y_days = as.numeric(followup_5y_days),
    event_5y = as.integer(event_5y)
  )

set.seed(analysis_seed)

idx_30d <- sample.int(nrow(short), size = floor(0.70 * nrow(short)))
idx_long <- sample.int(nrow(long), size = floor(0.70 * nrow(long)))

analysis_data <- list(
  short = short,
  long = long,
  train_30d = short[idx_30d, , drop = FALSE],
  test_30d = short[-idx_30d, , drop = FALSE],
  train_long = long[idx_long, , drop = FALSE],
  test_long = long[-idx_long, , drop = FALSE]
)

dir.create("derived", showWarnings = FALSE)
saveRDS(analysis_data, "derived/analysis_data.rds")

vapply(analysis_data[c("train_30d", "test_30d", "train_long", "test_long")], nrow, integer(1))

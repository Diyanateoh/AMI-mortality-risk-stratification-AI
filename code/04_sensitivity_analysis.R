library(dplyr)
library(pROC)
library(riskRegression)
library(ROSE)
library(survival)

analysis_seed <- 123L
d <- readRDS("derived/analysis_data.rds")
m <- readRDS("derived/models.rds")

rose_sample <- function(data, outcome, p, seed) {
  f <- reformulate(setdiff(names(data), outcome), response = outcome)
  ROSE(f, data = data, N = 2L * nrow(data), p = p, seed = seed)$data
}

status_at_horizon <- function(time, event, horizon) {
  ifelse(
    event == 1 & time <= horizon,
    1L,
    ifelse(time >= horizon, 0L, NA_integer_)
  )
}

metrics <- function(outcome, risk) {
  keep <- !is.na(outcome) & !is.na(risk)
  outcome <- as.integer(outcome[keep])
  risk <- as.numeric(risk[keep])

  roc_obj <- roc(outcome, risk, quiet = TRUE, direction = "<")
  threshold <- as.numeric(
    coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )[1, 1]
  )

  predicted <- as.integer(risk >= threshold)

  tp <- sum(predicted == 1L & outcome == 1L)
  tn <- sum(predicted == 0L & outcome == 0L)
  fp <- sum(predicted == 1L & outcome == 0L)
  fn <- sum(predicted == 0L & outcome == 1L)

  data.frame(
    AUC = as.numeric(auc(roc_obj)),
    AUC_lower = as.numeric(ci.auc(roc_obj)[1]),
    AUC_upper = as.numeric(ci.auc(roc_obj)[3]),
    threshold = threshold,
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    accuracy = (tp + tn) / (tp + tn + fp + fn)
  )
}

candidate_30d <- c(
  "age", "sex", "residential_area", "treatment_within_12h", "onset_season",
  "hypertension", "diabetes", "dyslipidemia", "bmi", "stemi",
  "multivessel_disease"
)

candidate_long <- c(
  "age", "sex", "residential_area", "treatment_within_12h",
  "hypertension", "diabetes", "dyslipidemia", "bmi", "stemi",
  "multivessel_disease"
)

fit_30d <- function(train) {
  glm(m$formulas$final_30d, data = train, family = binomial())
}

fit_1y <- function(train) {
  coxph(m$formulas$final_1y, data = train, x = TRUE, y = TRUE, model = TRUE)
}

fit_5y <- function(train) {
  coxph(m$formulas$final_5y, data = train, x = TRUE, y = TRUE, model = TRUE)
}

train_30d_original <- d$train_30d[, c(candidate_30d, "mortality_30d")]
train_1y_original <- d$train_long[, c(candidate_long, "followup_1y_days", "event_1y")]
train_5y_original <- d$train_long[, c(candidate_long, "followup_5y_days", "event_5y")]

train_30d_observed <- rose_sample(train_30d_original, "mortality_30d", 0.074, analysis_seed)
train_1y_observed <- rose_sample(train_1y_original, "event_1y", 0.077, analysis_seed)
train_5y_observed <- rose_sample(train_5y_original, "event_5y", 0.190, analysis_seed)

train_30d_equal <- rose_sample(train_30d_original, "mortality_30d", 0.5, analysis_seed)
train_1y_equal <- rose_sample(train_1y_original, "event_1y", 0.5, analysis_seed)
train_5y_equal <- rose_sample(train_5y_original, "event_5y", 0.5, analysis_seed)

models_30d <- list(
  original = fit_30d(train_30d_original),
  observed_event_rate = fit_30d(train_30d_observed),
  p_0_5 = fit_30d(train_30d_equal)
)

models_1y <- list(
  original = fit_1y(train_1y_original),
  observed_event_rate = fit_1y(train_1y_observed),
  p_0_5 = fit_1y(train_1y_equal)
)

models_5y <- list(
  original = fit_5y(train_5y_original),
  observed_event_rate = fit_5y(train_5y_observed),
  p_0_5 = fit_5y(train_5y_equal)
)

outcome_1y <- status_at_horizon(d$test_long$followup_1y_days, d$test_long$event_1y, 365)
outcome_5y <- status_at_horizon(d$test_long$followup_5y_days, d$test_long$event_5y, 1825)

evaluate_30d <- function(model, label) {
  risk <- predict(model, newdata = d$test_30d, type = "response")
  cbind(timepoint = "30-day", sampling = label, metrics(d$test_30d$mortality_30d, risk))
}

evaluate_1y <- function(model, label) {
  risk <- drop(predictRisk(model, newdata = d$test_long, times = 365))
  cbind(timepoint = "1-year", sampling = label, metrics(outcome_1y, risk))
}

evaluate_5y <- function(model, label) {
  risk <- drop(predictRisk(model, newdata = d$test_long, times = 1825))
  cbind(timepoint = "5-year", sampling = label, metrics(outcome_5y, risk))
}

sensitivity_results <- bind_rows(
  Map(evaluate_30d, models_30d, names(models_30d)),
  Map(evaluate_1y, models_1y, names(models_1y)),
  Map(evaluate_5y, models_5y, names(models_5y))
)

dir.create("outputs", showWarnings = FALSE)
write.csv(sensitivity_results, "outputs/sensitivity_results.csv", row.names = FALSE)

print(sensitivity_results)

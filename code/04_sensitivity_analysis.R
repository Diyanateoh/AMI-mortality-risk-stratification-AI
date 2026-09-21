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

classification_metrics <- function(outcome, risk) {
  keep <- !is.na(outcome) & !is.na(risk)
  outcome <- as.integer(outcome[keep])
  risk <- as.numeric(risk[keep])

  roc_obj <- roc(outcome, risk, quiet = TRUE, direction = "<")
  auc_ci <- ci.auc(roc_obj)

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

  list(
    roc = roc_obj,
    summary = data.frame(
      AUC = as.numeric(auc(roc_obj)),
      AUC_lower = as.numeric(auc_ci[1]),
      AUC_upper = as.numeric(auc_ci[3]),
      threshold = threshold,
      sensitivity = tp / (tp + fn),
      specificity = tn / (tn + fp),
      accuracy = (tp + tn) / (tp + tn + fp + fn)
    )
  )
}

fit_30d <- function(train) {
  glm(m$formulas$final_30d, data = train, family = binomial())
}

fit_1y <- function(train) {
  coxph(m$formulas$final_1y, data = train, x = TRUE, y = TRUE, model = TRUE)
}

fit_5y <- function(train) {
  coxph(m$formulas$final_5y, data = train, x = TRUE, y = TRUE, model = TRUE)
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

train_30d <- d$train_30d[, c(candidate_30d, "mortality_30d")]
train_1y <- d$train_long[, c(candidate_long, "followup_1y_days", "event_1y")]
train_5y <- d$train_long[, c(candidate_long, "followup_5y_days", "event_5y")]

training_sets_30d <- list(
  original = train_30d,
  observed_event_rate = rose_sample(train_30d, "mortality_30d", 0.074, analysis_seed),
  p_0_5 = rose_sample(train_30d, "mortality_30d", 0.5, analysis_seed)
)

training_sets_1y <- list(
  original = train_1y,
  observed_event_rate = rose_sample(train_1y, "event_1y", 0.077, analysis_seed),
  p_0_5 = rose_sample(train_1y, "event_1y", 0.5, analysis_seed)
)

training_sets_5y <- list(
  original = train_5y,
  observed_event_rate = rose_sample(train_5y, "event_5y", 0.190, analysis_seed),
  p_0_5 = rose_sample(train_5y, "event_5y", 0.5, analysis_seed)
)

models_30d <- lapply(training_sets_30d, fit_30d)
models_1y <- lapply(training_sets_1y, fit_1y)
models_5y <- lapply(training_sets_5y, fit_5y)

outcome_1y <- status_at_horizon(
  d$test_long$followup_1y_days,
  d$test_long$event_1y,
  365
)

outcome_5y <- status_at_horizon(
  d$test_long$followup_5y_days,
  d$test_long$event_5y,
  1825
)

evaluate_30d <- function(model, label) {
  risk <- predict(model, newdata = d$test_30d, type = "response")
  result <- classification_metrics(d$test_30d$mortality_30d, risk)

  list(
    row = data.frame(
      timepoint = "30-day",
      sampling = label,
      result$summary,
      check.names = FALSE
    ),
    roc = result$roc
  )
}

evaluate_1y <- function(model, label) {
  risk <- drop(predictRisk(model, newdata = d$test_long, times = 365))
  result <- classification_metrics(outcome_1y, risk)

  list(
    row = data.frame(
      timepoint = "1-year",
      sampling = label,
      result$summary,
      check.names = FALSE
    ),
    roc = result$roc
  )
}

evaluate_5y <- function(model, label) {
  risk <- drop(predictRisk(model, newdata = d$test_long, times = 1825))
  result <- classification_metrics(outcome_5y, risk)

  list(
    row = data.frame(
      timepoint = "5-year",
      sampling = label,
      result$summary,
      check.names = FALSE
    ),
    roc = result$roc
  )
}

res_30d <- Map(evaluate_30d, models_30d, names(models_30d))
res_1y <- Map(evaluate_1y, models_1y, names(models_1y))
res_5y <- Map(evaluate_5y, models_5y, names(models_5y))

all_results <- bind_rows(
  lapply(res_30d, `[[`, "row"),
  lapply(res_1y, `[[`, "row"),
  lapply(res_5y, `[[`, "row")
)

table_s2 <- all_results %>%
  filter(sampling %in% c("observed_event_rate", "p_0_5")) %>%
  mutate(
    sampling_p = case_when(
      sampling == "p_0_5" ~ 0.5,
      timepoint == "30-day" ~ 0.074,
      timepoint == "1-year" ~ 0.077,
      timepoint == "5-year" ~ 0.190
    )
  ) %>%
  select(
    timepoint, sampling_p, AUC, AUC_lower, AUC_upper,
    sensitivity, specificity, accuracy, threshold
  )

dir.create("outputs", showWarnings = FALSE)

write.csv(
  all_results,
  "outputs/sensitivity_all_models.csv",
  row.names = FALSE
)

write.csv(
  table_s2,
  "outputs/sensitivity_table_s2.csv",
  row.names = FALSE
)

pdf("outputs/sensitivity_roc_fig_s1.pdf", width = 10, height = 3.5)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

plot(
  res_30d$original$roc,
  legacy.axes = TRUE,
  main = "30-day",
  col = "red",
  lwd = 2
)
plot(
  res_30d$observed_event_rate$roc,
  add = TRUE,
  col = "blue",
  lwd = 2
)
legend(
  "bottomright",
  legend = c(
    sprintf("Original AUC = %.3f", res_30d$original$row$AUC),
    sprintf("ROSE AUC = %.3f", res_30d$observed_event_rate$row$AUC)
  ),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n",
  cex = 0.8
)

plot(
  res_1y$original$roc,
  legacy.axes = TRUE,
  main = "1-year",
  col = "red",
  lwd = 2
)
plot(
  res_1y$observed_event_rate$roc,
  add = TRUE,
  col = "blue",
  lwd = 2
)
legend(
  "bottomright",
  legend = c(
    sprintf("Original AUC = %.3f", res_1y$original$row$AUC),
    sprintf("ROSE AUC = %.3f", res_1y$observed_event_rate$row$AUC)
  ),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n",
  cex = 0.8
)

plot(
  res_5y$original$roc,
  legacy.axes = TRUE,
  main = "5-year",
  col = "red",
  lwd = 2
)
plot(
  res_5y$observed_event_rate$roc,
  add = TRUE,
  col = "blue",
  lwd = 2
)
legend(
  "bottomright",
  legend = c(
    sprintf("Original AUC = %.3f", res_5y$original$row$AUC),
    sprintf("ROSE AUC = %.3f", res_5y$observed_event_rate$row$AUC)
  ),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n",
  cex = 0.8
)

dev.off()

print(table_s2)

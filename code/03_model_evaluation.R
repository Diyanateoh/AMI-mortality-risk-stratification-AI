library(dcurves)
library(pROC)
library(riskRegression)
library(rms)
library(survival)

analysis_seed <- 123L
d <- readRDS("derived/analysis_data.rds")
m <- readRDS("derived/models.rds")

clip_probability <- function(x, eps = 1e-6) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

class_metrics <- function(outcome, risk) {
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
    auc = as.numeric(auc(roc_obj)),
    auc_ci = as.numeric(auc_ci[c(1, 3)]),
    threshold = threshold,
    confusion = matrix(
      c(tn, fn, fp, tp),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        Predicted = c("Survival", "Mortality"),
        Observed = c("Survival", "Mortality")
      )
    ),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    accuracy = (tp + tn) / (tp + tn + fp + fn)
  )
}

status_at_horizon <- function(time, event, horizon) {
  ifelse(
    event == 1 & time <= horizon,
    1L,
    ifelse(time >= horizon, 0L, NA_integer_)
  )
}

calibrate_binary <- function(outcome, risk, B = 1000L, seed = 123L) {
  risk <- clip_probability(risk)
  lp <- qlogis(risk)
  x <- data.frame(outcome = outcome, lp = lp)

  fit <- lrm(outcome ~ rcs(lp, 3), data = x, x = TRUE, y = TRUE)

  grid <- seq(
    max(0.001, min(risk)),
    min(0.999, max(risk)),
    length.out = 100
  )

  fitted <- plogis(
    predict(
      fit,
      newdata = data.frame(lp = qlogis(grid)),
      type = "lp"
    )
  )

  set.seed(seed)
  boot <- matrix(NA_real_, nrow = length(grid), ncol = B)

  for (b in seq_len(B)) {
    idx <- sample.int(nrow(x), replace = TRUE)

    fit_b <- try(
      lrm(
        outcome ~ rcs(lp, 3),
        data = x[idx, , drop = FALSE],
        x = TRUE,
        y = TRUE
      ),
      silent = TRUE
    )

    if (!inherits(fit_b, "try-error")) {
      boot[, b] <- plogis(
        predict(
          fit_b,
          newdata = data.frame(lp = qlogis(grid)),
          type = "lp"
        )
      )
    }
  }

  calibrated <- plogis(
    predict(fit, newdata = data.frame(lp = lp), type = "lp")
  )

  list(
    curve = data.frame(
      predicted = grid,
      observed = fitted,
      lower = apply(boot, 1, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(boot, 1, quantile, probs = 0.975, na.rm = TRUE)
    ),
    Eavg = mean(abs(calibrated - risk))
  )
}

calibrate_survival <- function(time, event, risk, horizon) {
  risk <- clip_probability(risk)
  lp <- qlogis(risk)
  x <- data.frame(time = time, event = event, lp = lp)

  fit <- cph(
    Surv(time, event) ~ rcs(lp, 3),
    data = x,
    x = TRUE,
    y = TRUE,
    surv = TRUE
  )

  grid <- seq(
    max(0.001, min(risk)),
    min(0.999, max(risk)),
    length.out = 100
  )

  pred_grid <- survest(
    fit,
    newdata = data.frame(lp = qlogis(grid)),
    times = horizon,
    conf.int = 0.95
  )

  pred_subject <- survest(
    fit,
    newdata = data.frame(lp = lp),
    times = horizon,
    conf.int = FALSE
  )

  calibrated <- 1 - as.numeric(pred_subject$surv)

  list(
    curve = data.frame(
      predicted = grid,
      observed = 1 - as.numeric(pred_grid$surv),
      lower = 1 - as.numeric(pred_grid$upper),
      upper = 1 - as.numeric(pred_grid$lower)
    ),
    Eavg = mean(abs(calibrated - risk))
  )
}

risk_30d <- predict(
  m$evaluation$model_30d,
  newdata = d$test_30d,
  type = "response"
)

risk_1y <- drop(
  predictRisk(
    m$evaluation$model_1y,
    newdata = d$test_long,
    times = 365
  )
)

risk_5y <- drop(
  predictRisk(
    m$evaluation$model_5y,
    newdata = d$test_long,
    times = 1825
  )
)

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

classification_30d <- class_metrics(
  d$test_30d$mortality_30d,
  risk_30d
)

classification_1y <- class_metrics(
  outcome_1y,
  risk_1y
)

classification_5y <- class_metrics(
  outcome_5y,
  risk_5y
)

brier_30d <- Score(
  list(Model = m$evaluation$model_30d),
  formula = mortality_30d ~ 1,
  data = d$test_30d,
  metrics = "brier",
  conf.int = 0.95,
  null.model = FALSE
)

brier_1y <- Score(
  list(Model = m$evaluation$model_1y),
  formula = Surv(followup_1y_days, event_1y) ~ 1,
  data = d$test_long,
  times = 365,
  metrics = "brier",
  conf.int = 0.95,
  null.model = FALSE,
  cens.method = "ipcw"
)

brier_5y <- Score(
  list(Model = m$evaluation$model_5y),
  formula = Surv(followup_5y_days, event_5y) ~ 1,
  data = d$test_long,
  times = 1825,
  metrics = "brier",
  conf.int = 0.95,
  null.model = FALSE,
  cens.method = "ipcw"
)

calibration_30d <- calibrate_binary(
  d$test_30d$mortality_30d,
  risk_30d,
  B = 1000L,
  seed = analysis_seed
)

calibration_1y <- calibrate_survival(
  d$test_long$followup_1y_days,
  d$test_long$event_1y,
  risk_1y,
  horizon = 365
)

calibration_5y <- calibrate_survival(
  d$test_long$followup_5y_days,
  d$test_long$event_5y,
  risk_5y,
  horizon = 1825
)

dca_30d_data <- transform(d$test_30d, model_risk = risk_30d)
dca_1y_data <- transform(d$test_long, model_risk = risk_1y)
dca_5y_data <- transform(d$test_long, model_risk = risk_5y)

dca_30d <- dca(
  mortality_30d ~ model_risk,
  data = dca_30d_data,
  thresholds = seq(0.01, 0.30, by = 0.01),
  label = list(model_risk = "Model")
)

dca_1y <- dca(
  Surv(followup_1y_days, event_1y) ~ model_risk,
  data = dca_1y_data,
  time = 365,
  thresholds = seq(0.01, 0.30, by = 0.01),
  label = list(model_risk = "Model")
)

dca_5y <- dca(
  Surv(followup_5y_days, event_5y) ~ model_risk,
  data = dca_5y_data,
  time = 1825,
  thresholds = seq(0.01, 0.30, by = 0.01),
  label = list(model_risk = "Model")
)

evaluation <- list(
  risk = list(
    risk_30d = risk_30d,
    risk_1y = risk_1y,
    risk_5y = risk_5y
  ),
  discrimination = list(
    result_30d = classification_30d[c("roc", "auc", "auc_ci")],
    result_1y = classification_1y[c("roc", "auc", "auc_ci")],
    result_5y = classification_5y[c("roc", "auc", "auc_ci")]
  ),
  brier = list(
    result_30d = brier_30d,
    result_1y = brier_1y,
    result_5y = brier_5y
  ),
  classification = list(
    result_30d = classification_30d,
    result_1y = classification_1y,
    result_5y = classification_5y
  ),
  calibration = list(
    result_30d = calibration_30d,
    result_1y = calibration_1y,
    result_5y = calibration_5y
  ),
  dca = list(
    result_30d = dca_30d,
    result_1y = dca_1y,
    result_5y = dca_5y
  )
)

saveRDS(evaluation, "derived/evaluation.rds")

print(classification_30d[c(
  "auc", "auc_ci", "threshold", "sensitivity", "specificity", "accuracy"
)])

print(classification_1y[c(
  "auc", "auc_ci", "threshold", "sensitivity", "specificity", "accuracy"
)])

print(classification_5y[c(
  "auc", "auc_ci", "threshold", "sensitivity", "specificity", "accuracy"
)])

print(brier_30d)
print(brier_1y)
print(brier_5y)

c(
  Eavg_30d = calibration_30d$Eavg,
  Eavg_1y = calibration_1y$Eavg,
  Eavg_5y = calibration_5y$Eavg
)

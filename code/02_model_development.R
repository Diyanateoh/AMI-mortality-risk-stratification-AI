library(caret)
library(MASS)
library(pROC)
library(riskRegression)
library(ROSE)
library(rms)
library(survival)

analysis_seed <- 123L
d <- readRDS("derived/analysis_data.rds")

rose_sample <- function(data, outcome, p, seed) {
  f <- reformulate(setdiff(names(data), outcome), response = outcome)
  ROSE(f, data = data, N = 2L * nrow(data), p = p, seed = seed)$data
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

rose_30d <- rose_sample(
  d$train_30d[, c(candidate_30d, "mortality_30d")],
  outcome = "mortality_30d",
  p = 0.074,
  seed = analysis_seed
)

rose_1y <- rose_sample(
  d$train_long[, c(candidate_long, "followup_1y_days", "event_1y")],
  outcome = "event_1y",
  p = 0.077,
  seed = analysis_seed
)

rose_5y <- rose_sample(
  d$train_long[, c(candidate_long, "followup_5y_days", "event_5y")],
  outcome = "event_5y",
  p = 0.190,
  seed = analysis_seed
)

full_30d <- glm(
  reformulate(candidate_30d, response = "mortality_30d"),
  data = rose_30d,
  family = binomial()
)

full_1y <- coxph(
  as.formula(
    paste(
      "Surv(followup_1y_days, event_1y) ~",
      paste(candidate_long, collapse = " + ")
    )
  ),
  data = rose_1y,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

full_5y <- coxph(
  as.formula(
    paste(
      "Surv(followup_5y_days, event_5y) ~",
      paste(candidate_long, collapse = " + ")
    )
  ),
  data = rose_5y,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

aic_30d <- stepAIC(full_30d, direction = "both", trace = FALSE)
aic_1y <- stepAIC(full_1y, direction = "both", trace = FALSE)
aic_5y <- stepAIC(full_5y, direction = "both", trace = FALSE)

final_30d <- mortality_30d ~ age + sex + stemi + bmi +
  multivessel_disease + residential_area + onset_season

final_1y <- Surv(followup_1y_days, event_1y) ~ age +
  treatment_within_12h + residential_area + hypertension + diabetes

final_5y <- Surv(followup_5y_days, event_5y) ~ age +
  treatment_within_12h + diabetes + sex + dyslipidemia

selected_terms <- function(model) sort(attr(terms(model), "term.labels"))

stopifnot(identical(
  selected_terms(aic_30d),
  sort(attr(terms(final_30d), "term.labels"))
))

stopifnot(identical(
  selected_terms(aic_1y),
  sort(attr(terms(final_1y), "term.labels"))
))

stopifnot(identical(
  selected_terms(aic_5y),
  sort(attr(terms(final_5y), "term.labels"))
))

glm_30d <- glm(final_30d, data = rose_30d, family = binomial())
cox_1y <- coxph(final_1y, data = rose_1y, x = TRUE, y = TRUE, model = TRUE)
cox_5y <- coxph(final_5y, data = rose_5y, x = TRUE, y = TRUE, model = TRUE)

cv_auc <- function(
    data, formula, event_var,
    model = c("logistic", "cox"),
    horizon = NULL, time_var = NULL,
    k = 5L, repeats = 20L, seed = 123L) {

  model <- match.arg(model)
  set.seed(seed)

  auc_values <- numeric(k * repeats)
  z <- 1L

  for (r in seq_len(repeats)) {
    folds <- caret::createFolds(
      data[[event_var]],
      k = k,
      list = TRUE
    )

    for (f in seq_len(k)) {
      train <- data[-folds[[f]], , drop = FALSE]
      valid <- data[folds[[f]], , drop = FALSE]

      if (length(unique(valid[[event_var]])) < 2L) {
        auc_values[z] <- NA_real_
        z <- z + 1L
        next
      }

      if (model == "logistic") {
        fit <- glm(formula, data = train, family = binomial())
        risk <- predict(fit, newdata = valid, type = "response")
        outcome <- valid[[event_var]]
      } else {
        fit <- coxph(formula, data = train)
        sf <- survfit(fit, newdata = valid)
        risk <- 1 - as.vector(summary(sf, times = horizon)$surv)
        outcome <- ifelse(
          valid[[time_var]] <= horizon & valid[[event_var]] == 1,
          1L,
          0L
        )
      }

      auc_values[z] <- as.numeric(
        auc(roc(outcome, risk, quiet = TRUE, direction = "<"))
      )
      z <- z + 1L
    }
  }

  auc_values <- auc_values[!is.na(auc_values)]
  mean_auc <- mean(auc_values)
  se_auc <- sd(auc_values) / sqrt(length(auc_values))

  c(
    mean_auc = mean_auc,
    lower_95 = mean_auc - qnorm(0.975) * se_auc,
    upper_95 = mean_auc + qnorm(0.975) * se_auc
  )
}

cv_stability <- rbind(
  "30-day" = cv_auc(
    rose_30d,
    final_30d,
    event_var = "mortality_30d",
    model = "logistic",
    k = 5L,
    repeats = 20L,
    seed = analysis_seed
  ),
  "1-year" = cv_auc(
    rose_1y,
    final_1y,
    event_var = "event_1y",
    model = "cox",
    horizon = 365,
    time_var = "followup_1y_days",
    k = 5L,
    repeats = 20L,
    seed = analysis_seed
  ),
  "5-year" = cv_auc(
    rose_5y,
    final_5y,
    event_var = "event_5y",
    model = "cox",
    horizon = 1825,
    time_var = "followup_5y_days",
    k = 5L,
    repeats = 20L,
    seed = analysis_seed
  )
)

dd_30d <- datadist(rose_30d)
options(datadist = "dd_30d")
lrm_30d <- lrm(final_30d, data = rose_30d, x = TRUE, y = TRUE)
nomogram_30d <- nomogram(
  lrm_30d,
  fun = plogis,
  funlabel = "30-day Mortality Risk",
  lp = FALSE
)

dd_1y <- datadist(rose_1y)
options(datadist = "dd_1y")
cph_1y <- cph(final_1y, data = rose_1y, x = TRUE, y = TRUE, surv = TRUE)
s_1y <- Survival(cph_1y)
nomogram_1y <- nomogram(
  cph_1y,
  fun = function(lp) s_1y(365, lp),
  funlabel = "1-year Survival Probability",
  lp = FALSE
)

dd_5y <- datadist(rose_5y)
options(datadist = "dd_5y")
cph_5y <- cph(final_5y, data = rose_5y, x = TRUE, y = TRUE, surv = TRUE)
s_5y <- Survival(cph_5y)
nomogram_5y <- nomogram(
  cph_5y,
  fun = function(lp) s_5y(1825, lp),
  funlabel = "5-year Survival Probability",
  lp = FALSE
)

options(datadist = NULL)

models <- list(
  formulas = list(
    final_30d = final_30d,
    final_1y = final_1y,
    final_5y = final_5y
  ),
  aic = list(
    model_30d = aic_30d,
    model_1y = aic_1y,
    model_5y = aic_5y
  ),
  evaluation = list(
    model_30d = glm_30d,
    model_1y = cox_1y,
    model_5y = cox_5y
  ),
  cross_validation = cv_stability,
  nomogram = list(
    model_30d = lrm_30d,
    model_1y = cph_1y,
    model_5y = cph_5y,
    nomogram_30d = nomogram_30d,
    nomogram_1y = nomogram_1y,
    nomogram_5y = nomogram_5y
  ),
  rose = list(
    data_30d = rose_30d,
    data_1y = rose_1y,
    data_5y = rose_5y
  )
)

saveRDS(models, "derived/models.rds")
print(cv_stability)

# ------------------------------------------------------------------------------
# Module: 4_evaluation_cox.R
# Purpose: Evaluate Cox model performance based on manuscript-specified metrics
# Models: Cox regression models for 1-year and 5-year mortality
# ------------------------------------------------------------------------------

# Description:
# This module evaluates prediction performance of Cox models using test set metrics
# explicitly reported in the manuscript:
# - Area under the ROC curve (AUC)
# - Optimal threshold (Youden's Index)
# - Sensitivity, Specificity
# - Positive Predictive Value (PPV), Negative Predictive Value (NPV), F1-score

library(survival)
library(caret)
library(pROC)

# External validation using test set
evaluate_external <- function(model, test_data, time_horizon, status_var, time_var) {
  # Predict survival probability at specified time point
  surv_fit <- survfit(model, newdata = test_data)
  surv_probs <- summary(surv_fit, times = time_horizon)$surv
  event_probs <- 1 - surv_probs

  # Define binary event outcome based on true survival status
  true_outcome <- ifelse(test_data[[time_var]] >= time_horizon & test_data[[status_var]] == 0, 0, 1)

  # ROC curve and optimal threshold
  roc_obj <- roc(true_outcome, event_probs, quiet = TRUE)
  best_thresh <- coords(roc_obj, "best", ret = "threshold")$threshold
  pred_binary <- ifelse(event_probs > best_thresh, 1, 0)

  # Confusion matrix and derived performance metrics
  cm <- confusionMatrix(factor(pred_binary, levels = c(0, 1)), factor(true_outcome, levels = c(0, 1)))

  list(
    AUC = auc(roc_obj),
    CI = ci.auc(roc_obj, method = "delong"),
    Threshold = best_thresh,
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Accuracy = cm$overall["Accuracy"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    F1 = cm$byClass["F1"]
  )
}

# Example usage:
# model_1y <- coxph(Surv(FU_1y, event_status_1y) ~ Age_C + Time + DM + HTN, data = train_data)
# model_5y <- coxph(Surv(FU_5y, event_status_5y) ~ Age_C + Time + DM + HTN, data = train_data)
#
# metrics_1y <- evaluate_external(model_1y, test_data, time_horizon = 365, status_var = "event_status_1y", time_var = "FU_1y")
# metrics_5y <- evaluate_external(model_5y, test_data, time_horizon = 1825, status_var = "event_status_5y", time_var = "FU_5y")

# Output can be printed or stored as needed for reporting.

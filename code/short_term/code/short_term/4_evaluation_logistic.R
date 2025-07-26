# ------------------------------------------------------------------------------
# Module: 4_evaluation_logistic.R
# Purpose: Evaluate model performance of 30-day mortality prediction
# Input: Logistic model from Module 3, external test set
# Output: AUC, optimal threshold (Youden Index), accuracy, sensitivity, specificity,
#         positive predictive value (PPV), negative predictive value (NPV), F1 score
# ------------------------------------------------------------------------------

# Description:
# This script evaluates the predictive performance of the logistic regression model
# constructed for 30-day mortality following acute myocardial infarction (AMI).
# Using an external test set, the model’s discrimination is assessed via AUC and ROC.
# Optimal threshold classification (Youden's index) is used to calculate accuracy,
# sensitivity, specificity, and additional metrics including PPV, NPV, and F1 score.
# These metrics are consistent with those reported in the manuscript.

# Note: Observed 30-day mortality rate in the test cohort was 7.4%
# ------------------------------------------------------------------------------

library(pROC)
library(caret)

# Load external test set (preprocessed in same format as training)
test_data <- read.csv("your_external_test_data.csv")

# Predict probabilities using logistic model
predicted_probs <- predict(logit_model, newdata = test_data, type = "fitted")

# ROC curve and AUC
roc_obj <- roc(test_data$event_status_30, predicted_probs)
auc_value <- auc(roc_obj)
cat("\nAUC:", round(auc_value, 3), "\n")

# Plot ROC curve
plot(roc_obj, print.auc = TRUE, col = "steelblue", main = "ROC Curve - 30-day Mortality")

# Optimal threshold (Youden’s Index)
opt_thresh <- coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)
threshold_value <- round(opt_thresh["threshold"], 3)
cat("\nOptimal threshold (Youden):", threshold_value, "\n")

# Classify using optimal threshold
predicted_class <- ifelse(predicted_probs >= threshold_value, 1, 0)

# Confusion matrix and metrics
conf_mat <- confusionMatrix(factor(predicted_class),
                            factor(test_data$event_status_30),
                            positive = "1")

cat("\nConfusion Matrix:\n")
print(conf_mat$table)

cat("\nAccuracy:", round(conf_mat$overall["Accuracy"], 3))
cat("\nSensitivity:", round(conf_mat$byClass["Sensitivity"], 3))
cat("\nSpecificity:", round(conf_mat$byClass["Specificity"], 3))
cat("\nPositive Predictive Value (PPV):", round(conf_mat$byClass["Pos Pred Value"], 3))
cat("\nNegative Predictive Value (NPV):", round(conf_mat$byClass["Neg Pred Value"], 3))
cat("\nF1 Score:", round(conf_mat$byClass["F1"], 3))

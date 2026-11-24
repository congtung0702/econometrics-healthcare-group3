# Load Data ----
data = read.csv("Fake_CRC_riskdata.csv")

# Libraries used
library(dplyr)
library(pROC)

# Functions ----

# Calibration plot, intercept, slope and Wald tests
calibration_plot = function(actual, predicted, name = "Model", bins = 10) {
  
  # Divide observations into equally sized bins
  calibration = data.frame(
    actual = actual,
    predicted = predicted,
    bin = ntile(predicted, bins)
  ) %>%
    group_by(bin) %>%
    summarise(
      predicted = mean(predicted),
      observed = mean(actual),
      n = n(),
      .groups = "drop"
    )
  
  print(calibration)
  
  # Calibration plot
  max_value = max(calibration$predicted, calibration$observed)
  
  plot(
    calibration$predicted,
    calibration$observed,
    xlab = "Predicted CRC probability",
    ylab = "Observed CRC proportion",
    main = paste("Calibration Plot -", name),
    pch = 19,
    xlim = c(0, max_value),
    ylim = c(0, max_value)
  )
  
  # Perfect calibration line
  abline(0, 1, col = "red", lty = 2)
  
  # Calibration intercept and slope
  logit_pred = qlogis(
    pmin(pmax(predicted, 1e-8), 1 - 1e-8)
  )
  
  cal_model = glm(
    actual ~ logit_pred,
    family = binomial(link = "logit")
  )
  
  # Estimated coefficients and standard errors
  estimates = coef(cal_model)
  se = sqrt(diag(vcov(cal_model)))
  
  # Null hypotheses: intercept = 0, slope = 1
  null_values = c(0, 1)
  
  # Wald tests
  z_values = (estimates - null_values) / se
  p_values = 2 * pnorm(-abs(z_values))
  
  # Results
  cal_results = data.frame(
    Estimate = estimates,
    Std_Error = se,
    H0 = null_values,
    z_value = z_values,
    p_value = p_values
  )
  
  print(round(cal_results, 4))
  
  return(invisible(list(
    calibration = calibration,
    model = cal_model,
    results = cal_results
  )))
}

# ROC curve and AUROC
roc_plot = function(actual, predicted, name = "Model") {
  
  # ROC curve
  roc_model = roc(
    response = actual,
    predictor = predicted,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  # Calculate AUROC
  auc_value = as.numeric(auc(roc_model))
  
  # Plot ROC curve
  plot(
    roc_model,
    legacy.axes = TRUE,
    xlim = c(1, 0),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    asp = NA,
    col = "blue",
    lwd = 3,
    main = paste("ROC Curve -", name),
    xlab = "False Positive Rate (1 - Specificity)",
    ylab = "True Positive Rate (Sensitivity)",
    identity = TRUE,
    identity.col = "red",
    identity.lty = 2
  )
  
  # Add AUROC to plot
  legend(
    "bottomright",
    legend = paste("AUROC =", round(auc_value, 4)),
    bty = "n"
  )
  
  print(paste(name, "- AUROC:", round(auc_value, 4)))
  
  return(auc_value)
}

# Question A. ----

# No missing values
colSums(is.na(data))
# General descriptive statistics
summary(data)
# Standard deviations
round(sqrt(diag(var(data))), 3)
# Correlation check
round(cor(data), 3)

# Question B & C. ----
set.seed(2026)

# Randomly select 70% of observations
train_index = sample(1:nrow(data), size = 0.7 * nrow(data))

# Split data
train = data[train_index, ]
test = data[-train_index, ]

# Baseline Model

model = glm(
  CRC ~ male + age + FIT_1 + FIT_2,
  family = binomial(link = "logit"),
  data = train
)

summary(model)

# Predictions
test$predicted = predict(model, newdata = test, type = "response")

# Question D: Calibration ----
cal1 = calibration_plot(test$CRC, test$predicted, "Baseline Model")

# Question E: ROC + AUROC ----
auc1 = roc_plot(test$CRC, test$predicted, "Baseline Model")

# Question H: Interaction Model ----

model2 = glm(
  CRC ~ (male + age + FIT_1 + FIT_2)^2,
  family = binomial(link = "logit"),
  data = train
)

summary(model2)

# Predictions
test$predicted2 = predict(model2, newdata = test, type = "response")

# Calibration
cal2 = calibration_plot(test$CRC, test$predicted2, "Interaction Model")

# ROC + AUROC
auc2 = roc_plot(test$CRC, test$predicted2, "Interaction Model")


# Model Comparison 

comparison = data.frame(
  Model = c("Baseline", "Interactions"),
  AUROC = c(auc1, auc2),
  Calibration_Intercept = c(
    coef(cal1$model)[1],
    coef(cal2$model)[1]
  ),
  Calibration_Slope = c(
    coef(cal1$model)[2],
    coef(cal2$model)[2]
  )
)

print(comparison)

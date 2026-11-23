# Load Data
data = read.csv("Fake_CRC_riskdata.csv")

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
train = data[train_index,]
test = data[-train_index,]

# Train logigistic regression and show coeffients
model = glm(CRC ~ male + age + FIT_1 + FIT_2, family = binomial(link = "logit"), data = train)
summary(model)

# Question D ----
# Predict probabilities
test$predicted = predict(
  model,
  newdata = test,
  type = "response"
)

# Divide test observations into 10 equally sized bins
library(dplyr)

test$bin = ntile(test$predicted, 10)

# Calculate observed and predicted probabilities per bin
calibration = test %>%
  group_by(bin) %>%
  summarise(
    predicted = mean(predicted),
    observed = mean(CRC),
    n = n()
  )

print(calibration)

# Calibration plot
plot(
  calibration$predicted,
  calibration$observed,
  xlab = "Predicted CRC probability",
  ylab = "Observed CRC proportion",
  main = "Calibration plot (10 bins)",
  pch = 19,
  xlim = c(0, max(calibration$predicted,
                  calibration$observed)),
  ylim = c(0, max(calibration$predicted,
                  calibration$observed))
)

# Perfect calibration line
abline(0, 1, col = "red", lty = 2)

# Calibration intercept and slope
test$logit_pred = qlogis(
  pmin(pmax(test$predicted, 1e-8), 1 - 1e-8)
)

cal_model = glm(
  CRC ~ logit_pred,
  family = binomial(link = "logit"),
  data = test
)

summary(cal_model)

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


# Question E. ----
library(pROC)

# ROC curve using test dataset
roc_model = roc(
  response = test$CRC,
  predictor = test$predicted,
  levels = c(0, 1),
  direction = "<"
)

# Area under the ROC curve
auc(roc_model)

# ROC-Curve plot 
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
  main = "ROC Curve",
  xlab = "False Positive Rate (1 - Specificity)",
  ylab = "True Positive Rate (Sensitivity)",
  identity = TRUE,
  identity.col = "red",
  identity.lty = 2
)

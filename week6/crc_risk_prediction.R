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


# Question B. ----

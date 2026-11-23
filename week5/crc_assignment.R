# Load Data
data = read.csv("MISCANPosterior.csv")

# Question A. ----
# Define selected pairs (based on earlier visual inspection)
selected_pairs = list(
  c("dwell_pcl1", "dwell_pcl2"),
  c("trans_clin_loc_12", "trans_clin_loc_14"),
  c("trans_clin_loc_22", "trans_clin_loc_24"),
  c("trans_clin_loc_32", "trans_clin_loc_34"),
  c("trans_clin_age_Pcl2100", "trans_clin_age_Pcl20"),
  c("trans_clin_age_Pcl3100", "trans_clin_age_Pcl30")
)

# Keep only selected variables
selected = data[, unique(unlist(selected_pairs))]

# Correlations for the six selected pairs
correlations = data.frame(
  Variable1 = sapply(selected_pairs, `[`, 1),
  Variable2 = sapply(selected_pairs, `[`, 2),
  Correlation = sapply(selected_pairs, function(p) {
    cor(data[[p[1]]], data[[p[2]]], use = "complete.obs")
  })
)

print(correlations)

# Six scatterplots in one figure
par(mfrow = c(2, 3))

for (i in 1:length(selected_pairs)) {
  p = selected_pairs[[i]]
  
  plot(
    data[[p[1]]], data[[p[2]]],
    xlab = p[1],
    ylab = p[2],
    main = paste("r =", round(correlations$Correlation[i], 3)),
    pch = 16,
    cex = 0.3
  )
}

par(mfrow = c(1, 1))

# Questions B, C and D ----

# Set up 2x2 figure
par(mfrow = c(2, 2),
    mar = c(4.5, 4.5, 3, 1),
    las = 1)

# Function for scatterplots
cor_plot = function(x, y, xlab, main) {
  
  # Calculate correlation
  r = cor(x, y, use = "complete.obs")
  
  # Scatterplot
  plot(x, y,
       main = main,
       xlab = xlab,
       ylab = "FIT sensitivity",
       pch = 16,
       cex = 0.5,
       col = rgb(0.4, 0.4, 0.4, 0.3))
  
  # Black regression line (optional)
  abline(lm(y ~ x), col = "black", lwd = 2)
  
  # Correlation coefficient
  legend("topright",
         legend = paste0("r = ", round(r, 3)),
         bty = "n",
         cex = 1.1)
}

# Question B: Stage 1 and Stage 2
cor_plot(data$dwell_pcl1, data$sensitivity_47FIT_b,
         "Dwell time stage 1", "B: Stage 1")

cor_plot(data$dwell_pcl2, data$sensitivity_47FIT_b,
         "Dwell time stage 2", "B: Stage 2")

# Question C: Stage 4
cor_plot(data$dwell_pcl4, data$sensitivity_47FIT_b,
         "Dwell time stage 4", "C: Stage 4")

# Question D: Hazard variance
cor_plot(data$hazard_variance, data$sensitivity_47FIT_b,
         "Hazard variance", "D: Hazard variance")
## Admixture
# K 1-10, 10 replicates each K for ADMIXTURE was done by hpc.
# ADMIXTURE is based on a maximum likelihood estimation framework and shares its foundational model with STRUCTURE, but it uses a faster numerical optimization algorithm.
# It includes robust cross-validation (CV) to help select the optimal number of ancestral populations (K) based on predictive accuracy.

# Choose the best K
#Use ADMIXTURE’s cross-validation procedure. A good value of K will exhibit a low cross-validation error compared to other K values. 
setwd("/scratch/fm2371/Winchell/SLF_spring23/SLF_spring23/raw.vcf/Final/Admixture")

# Read in the data
data <- list.files(pattern = "Final\\.log.*\\.out", full.names = TRUE)

# Use lapply to read in all log files at once
bigData <- lapply(data, readLines)

# Find the line that starts with "CV error"
library(stringr)
foundset_cv <- sapply(bigData, function(x) x[which(str_detect(x, "^CV error"))])

# Extract the numeric value from the "CV error" line
cv_errors <- as.numeric(sub("CV error \\(K=\\d+\\): (-?\\d+\\.\\d+).*", "\\1", foundset_cv))

# Create a data frame with K values and corresponding CV errors
cv_df <- data.frame(
  K = rep(1:10, each = 10),
  cv_error = cv_errors
)

# Calculate the mean CV error for each K
cv_summary <- aggregate(cv_error ~ K, data = cv_df, FUN = mean)

cv_summary_plot <- ggplot(cv_summary, aes(x = K, y = cv_error)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(color = "#E74C3C", size = 3) +
  #geom_text(aes(label = K), nudge_x = 0.2, nudge_y = 0.002, size = 4, color = "#2C3E50") +
  labs(x = "K", y = "Mean CV Error") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  ggtitle("Mean CV Error over 10 replicates") +
  scale_x_continuous(breaks = seq(1, max(cv_summary$K), by = 1)) +
  scale_y_continuous(limits = c(min(cv_summary$cv_error) - 0.005, max(cv_summary$cv_error) + 0.005))

# Save the plot. Supplementary Figure S3.C
ggsave(filename="cv_summary_plot.pdf", plot=cv_summary_plot, width=12, height=8, units="in")
ggsave(filename="cv_summary_plot.png", plot=cv_summary_plot,dpi = 300)
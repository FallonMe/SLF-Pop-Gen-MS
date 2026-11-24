# Stats for momi2.

library(dplyr)
library(readr)

# Read and clean AIC file
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Momi2")
aic_data_raw <- readLines("model_aic_and_likelihoods.txt")

# Filter out lines with bad values or repeated headers
aic_data_clean <- aic_data_raw[!grepl("^Model|^\\s*$", aic_data_raw)]

# Read it as a data frame (forcing proper parsing)
aic_data <- read.csv(text = aic_data_clean, header = FALSE, stringsAsFactors = FALSE)
colnames(aic_data) <- c("Replicate", "Model", "AIC", "LogLikelihood", "DeltaAIC", "AICWeight")

# Coerce numeric columns (bad values become NA)
aic_data <- aic_data %>%
  mutate(across(c(AIC, LogLikelihood, DeltaAIC, AICWeight), ~as.numeric(trimws(.))))

# Remove rows with any NAs
aic_data <- aic_data %>% filter(complete.cases(.))

# Optional: Add model labels
model_labels <- c(
  "Model1_Urban_Korea" = "Model1_Urban_Korea",
  "Model2_Urban_Direct" = "Model2_Urban_Direct",
  "Model3_Rural_Korea" = "Model3_Rural_Korea",
  "Model4_Rural_Direct" = "Model4_Rural_Direct",
  "Model5_Independent_Urban_Korea" = "Model5_Independent_Urban_Korea",
  "Model6_Independent_Urban_Direct" = "Model6_Independent_Urban_Direct",
  "Model7_Independent_Rural_Korea" = "Model7_Independent_Rural_Korea",
  "Model8_Independent_Rural_Direct" = "Model8_Independent_Rural_Direct"
)
aic_data$Model <- factor(aic_data$Model, levels = names(model_labels), labels = model_labels)

# Summarize
mean_values_df <- aic_data %>%
  group_by(Model) %>%
  summarise(
    Mean_AIC = mean(AIC),
    SD_AIC = sd(AIC),
    Mean_LogLikelihood = mean(LogLikelihood),
    Mean_DeltaAIC = mean(DeltaAIC),
    Mean_AICWeight = mean(AICWeight),
    N = n(),
    CI95_AIC_Lower = Mean_AIC - qt(0.975, df = N - 1) * SD_AIC / sqrt(N),
    CI95_AIC_Upper = Mean_AIC + qt(0.975, df = N - 1) * SD_AIC / sqrt(N),
    .groups = "drop"
  ) %>%
  arrange(Mean_AIC)

# Output
print(mean_values_df)
write.csv(mean_values_df, "mean_model_AIC_summary.csv", row.names = FALSE)
cat("Best-fitting model:", mean_values_df$Model[1], "\n")


# Calculate the average AIC serial stats and CI/Sd after over 700 bootstraps, find the best model of fit.

# The bootstrap analysis
# ```{r}
# Load necessary library
library(dplyr)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Momi2")

# Read as a comma-separated CSV
params <- read.csv("best_model_params.txt", header = TRUE, stringsAsFactors = FALSE)

# Check
head(params)

# After reading in params:
nrow(params)  # 7065

# Reassign Replicate numbers: 9 rows per replicate
params$Replicate <- rep(1:785, each = 9)

# Check
table(params$Replicate)  # Each replicate should have 9 rows

# 2. Fix Replicate numbers
params$Replicate <- rep(1:(nrow(params)/9), each = 9)

# Now convert to list of named lists (bootstrap_resultsA)
bootstrap_resultsA <- params %>%
  filter(BestModel == "Model5_Independent_Urban_Korea") %>%
  group_by(Replicate) %>%
  summarise(param_list = list(setNames(Value, Parameter))) %>%
  pull(param_list)

# Check the first replicate
bootstrap_resultsA[[1]]

# Define a helper to format one parameter set like Python dictionary
format_dict <- function(lst) {
  items <- paste0("'", names(lst), "': ", lst)
  paste0("({", paste(items, collapse = ", "), "})")
}

# Apply the formatter to each replicate's parameter set
formatted_list <- sapply(bootstrap_resultsA, format_dict)

# Combine into final Python-style list string
cat("bootstrap_resultsA=[", paste(formatted_list, collapse = ",\n "), "]\n")

# Save the Python-style bootstrap list to a text file
output_path <- "bootstrap_resultsA_list.txt"
writeLines(c("bootstrap_resultsA=[", paste(formatted_list, collapse = ",\n "), "]"), con = output_path)

# Convert list of named lists to data frame
param <- do.call(rbind, lapply(bootstrap_resultsA, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
param <- as.data.frame(lapply(param, as.numeric))  # Ensure numeric

# Function to calculate 95% CI
ci95 <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  n <- sum(!is.na(x))
  error <- qt(0.975, df = n - 1) * s / sqrt(n)
  return(c(mean = m, sd = s, min = min(x), max = max(x), lower_CI = m - error, upper_CI = m + error))
}

# Apply to each parameter
summary_stats <- as.data.frame(t(sapply(param, ci95)))

# Round for better presentation
summary_stats_rounded <- round(summary_stats, 4)

# Show results
summary_stats_rounded

write.csv(summary_stats_rounded, "model5_summary_stats_rounded.csv")
write.csv(summary_stats, "model5_summary_stats_full.csv")


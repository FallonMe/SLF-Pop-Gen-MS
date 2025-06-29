# ─────────────── 1. Load Libraries ───────────────
library(ggplot2)
library(dplyr)
library(boot)
library(readr)

# ─────────────── 2. Load Data ───────────────
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/ANGSD")

# Custom color schemes
country_colors <- c("USA" = "#3690c0", "CN" = "#c861a5")
habitat_colors <- c("urban" = "#c861a5", "rural" = "#009E73")

# Metadata and ANGSD output
metadata <- read_csv("../Final_SLF_MS_Info.csv")
het_data <- read_delim("ANGSD_heterozygosity.txt")

# ─────────────── 3. Merge and Arrange ───────────────
het_data <- merge(het_data, metadata, by.x = "SampleID", by.y = "SLF_ID")
het_data <- het_data %>%
  arrange(lat) %>%
  mutate(Location = factor(Location, levels = unique(Location)))

# ─────────────── 4. Summary Stats ───────────────
summary_stats <- het_data %>%
  group_by(Country) %>%
  summarize(
    mean_Ho = mean(Heterozygosity, na.rm = TRUE),
    sd_Ho   = sd(Heterozygosity, na.rm = TRUE),
    min_Ho  = min(Heterozygosity, na.rm = TRUE),
    max_Ho  = max(Heterozygosity, na.rm = TRUE),
    n       = n()
  )
write.csv(summary_stats, "ANGSD_summary_stats_Country.csv")
print(summary_stats)

# ─────────────── 5. Statistical Tests ───────────────
# Welch's t-test
t_test <- t.test(Heterozygosity ~ Country, data = het_data)
print(t_test)

# F-test for variance
print(var.test(Heterozygosity ~ Country, data = het_data))

# ─────────────── 6. Bootstrap CI of Mean Difference ───────────────
boot_diff_means <- function(data, indices) {
  d <- data[indices, ]
  mean(d$Heterozygosity[d$Country == "CN"]) - mean(d$Heterozygosity[d$Country == "USA"])
}
boot_res <- boot(het_data, boot_diff_means, R = 10000)
print(boot.ci(boot_res, type = "perc"))

# ─────────────── 7. Resampling Test for Subsampling Bias ───────────────
# Bootstrap mean heterozygosity of USA (sample size = 20)
set.seed(123)
usa_vals <- het_data %>% filter(Country == "USA") %>% pull(Heterozygosity)
cn_vals  <- het_data %>% filter(Country == "CN")  %>% pull(Heterozygosity)

subsample_means <- replicate(1000, mean(sample(usa_vals, size = 20, replace = FALSE)))
sub_mean <- mean(subsample_means)
sub_se   <- sd(subsample_means)

# Save for downstream
write.csv(data.frame(bootstrap_usa_means = subsample_means), "USA_Het_Subsample_Means.csv")

# ─────────────── 8. Plot: Violin with Subsample SE ───────────────
subsample_df <- data.frame(
  Country = "USA",
  Heterozygosity = sub_mean,
  ymin = sub_mean - sub_se,
  ymax = sub_mean + sub_se
)

p <- ggplot(het_data, aes(x = Country, y = Heterozygosity)) +
  geom_violin(aes(fill = Country), alpha = 0.7) +
  geom_point(aes(color = Country), alpha = 0.8, size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 95, alpha = 0.5) +
  geom_errorbar(data = subsample_df, aes(x = Country, ymin = ymin, ymax = ymax),
                color = "red", width = 0.2, size = 1.5) +
  geom_point(data = subsample_df, aes(x = Country, y = Heterozygosity),
             color = "red", size = 3, shape = 21, fill = "red") +
  scale_fill_manual(values = country_colors, guide = FALSE) +
  scale_color_manual(values = country_colors, guide = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Genome-wide Heterozygosity") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

ggsave("Country_ANGSD_Heterozygosity_violin_plot_errorbar_with_subsample.pdf", plot = p, width = 8, height = 6)
ggsave("Country_ANGSD_Heterozygosity_violin_plot_errorbar_with_subsample.png", plot = p, dpi = 300)

# ─────────────── 9. Density Visualization of Subsamples ───────────────
ggplot(data.frame(het = subsample_means), aes(het)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = sub_mean, color = "black") +
  geom_vline(xintercept = mean(usa_vals), color = "red") +
  geom_vline(xintercept = mean(cn_vals), color = "blue") +
  theme_bw() +
  labs(x = "USA Heterozygosity (subsampled, n=20)", y = "Density")

# ─────────────── 10. Export Merged Data ───────────────
write.csv(het_data, "ANGSD_merged_heterozygosity.csv", row.names = FALSE)

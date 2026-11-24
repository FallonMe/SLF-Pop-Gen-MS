# ========================================================================
# Comprehensive Runs of Homozygosity (ROH) Analysis for Spotted Lanternfly
# Without Filtering
# ========================================================================

# Set working directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/RoH_new/NO_filter")

# Load necessary libraries
library(ggplot2)       # For plotting
library(dplyr)         # For data manipulation
library(readr)         # For reading/writing files
library(tidyverse)     # For data manipulation and visualization
library(ggpubr)        # For statistical testing visualization
library(car)           # For ANOVA testing
library(lmtest)        # For robust statistical tests

# ========================================================================
# 1. SET UP COLOR SCHEMES AND LOAD DATA
# ========================================================================

# Define color schemes for different grouping variables
country_colors <- c(
  "USA" = "#3690c0",   # Blue for USA
  "CN" = "#c861a5"     # Pink for China
)

habitat_colors <- c(
  "urban" = "#c861a5", # Pink for urban areas
  "rural" = "#009E73"  # Green for rural areas
)

# More detailed location colors
location_colors <- c(
  "Bronx" = "#252525", 
  "Brooklyn" = "#e6ab02",
  "Camden" = "#a65628",
  "Connecticut" = "#e41a1c",
  "Governor Island" = "#08589e",
  "Manhattan" = "#c861a5",
  "Staten Island" = "#3690c0", 
  "West Orange" = "#009E73",
  "Jinqiao" = "#ff7f00",
  "Sheshan" = "#6a3d9a",
  "Tianma" = "#4daf4a"
)

# Define colors for ROH size classes
size_class_colors <- c(
  "0.5-1 Mb" = "#009E73",  # Green for small ROHs
  "1-5 Mb" = "#e6ab02",    # Yellow for medium ROHs
  ">5 Mb" = "#c861a5"      # Pink for large ROHs
)

# Load population metadata
popdata <- read.csv("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/Final_SLF_MS_Info.csv", header=TRUE)

# Read PLINK ROH data - using the NOFilter_roh_results.hom file
plink_roh_no_filter <- read.delim('ROH_noSexChr_results.hom', sep='')

# ========================================================================
# 2. PROCESS ROH DATA AND CATEGORIZE BY SIZE CLASS
# ========================================================================

# Categorize ROH segments into size classes
plink_roh_size_no_filter <- plink_roh_no_filter %>%
  mutate(size_class = case_when(
    KB >= 500 & KB < 1000 ~ "0.5-1 Mb",   # Small ROHs
    KB >= 1000 & KB < 5000 ~ "1-5 Mb",    # Medium ROHs
    KB >= 5000 ~ ">5 Mb"                  # Large ROHs
  ))

# Order size classes for plotting
plink_roh_size_no_filter$size_class <- factor(plink_roh_size_no_filter$size_class, 
                                              levels=c("0.5-1 Mb", "1-5 Mb", ">5 Mb"))

# Define genome length covered by SNPs for FROH calculation
LAUTO <- 1998219946 / 1000  # Convert to kb

# ========================================================================
# 3. CALCULATE ROH STATISTICS BY INDIVIDUAL
# ========================================================================

# Calculate count and total length of ROH segments for each individual
roh_individual_no_filter <- plink_roh_size_no_filter %>%
  group_by(FID, size_class) %>%
  summarise(
    roh_count = n(),               # Number of segments
    total_length = sum(KB),        # Total length in kb
    mean_length = mean(KB),        # Mean segment length
    min_length = min(KB),          # Minimum segment length
    max_length = max(KB),          # Maximum segment length
    .groups = 'drop'
  ) %>%
  left_join(popdata, by = c("FID" = "SLF_ID"))  # Add population info

# Calculate total number of ROH segments per individual (across all size classes)
total_roh_individual_no_filter <- plink_roh_size_no_filter %>%
  group_by(FID) %>%
  summarise(
    total_roh_count = n(),         # Total number of segments
    total_length = sum(KB),        # Total length of all segments
    mean_length = mean(KB),        # Mean segment length
    min_length = min(KB),          # Minimum segment length
    max_length = max(KB),          # Maximum segment length
    .groups = 'drop'
  ) %>%
  mutate(FROH = total_length / LAUTO) %>%  # Calculate FROH
  left_join(popdata, by = c("FID" = "SLF_ID"))  # Add population info

# ========================================================================
# 4. CALCULATE SUMMARY STATISTICS BY AREA AND COUNTRY
# ========================================================================

# ROH count statistics by area and size class
roh_area_size_stats_no_filter <- roh_individual_no_filter %>%
  group_by(area, size_class) %>%
  summarise(
    mean_count = mean(roh_count),                      # Mean count
    median_count = median(roh_count),                  # Median count
    sd_count = sd(roh_count),                          # Standard deviation
    se_count = sd(roh_count)/sqrt(n()),                # Standard error
    min_count = min(roh_count),                        # Minimum count
    max_count = max(roh_count),                        # Maximum count
    sample_size = n(),                                 # Sample size
    mean_length = mean(total_length),                  # Mean total length
    median_length = median(total_length),              # Median total length
    sd_length = sd(total_length),                      # Standard deviation of length
    se_length = sd(total_length)/sqrt(n()),            # Standard error of length
    min_length = min(total_length),                    # Minimum total length
    max_length = max(total_length),                    # Maximum total length
    ci_lower_count = mean_count - qt(0.975, sample_size-1) * se_count,  # 95% CI lower bound
    ci_upper_count = mean_count + qt(0.975, sample_size-1) * se_count,  # 95% CI upper bound
    .groups = 'drop'
  )

# ROH count statistics by country and size class
roh_country_size_stats_no_filter <- roh_individual_no_filter %>%
  group_by(Country, size_class) %>%
  summarise(
    mean_count = mean(roh_count),                      # Mean count
    median_count = median(roh_count),                  # Median count
    sd_count = sd(roh_count),                          # Standard deviation
    se_count = sd(roh_count)/sqrt(n()),                # Standard error
    min_count = min(roh_count),                        # Minimum count
    max_count = max(roh_count),                        # Maximum count
    sample_size = n(),                                 # Sample size
    mean_length = mean(total_length),                  # Mean total length
    median_length = median(total_length),              # Median total length
    sd_length = sd(total_length),                      # Standard deviation of length
    se_length = sd(total_length)/sqrt(n()),            # Standard error of length
    min_length = min(total_length),                    # Minimum total length
    max_length = max(total_length),                    # Maximum total length
    ci_lower_count = mean_count - qt(0.975, sample_size-1) * se_count,  # 95% CI lower bound
    ci_upper_count = mean_count + qt(0.975, sample_size-1) * se_count,  # 95% CI upper bound
    .groups = 'drop'
  )

# FROH statistics by area
froh_area_stats_no_filter <- total_roh_individual_no_filter %>%
  group_by(area) %>%
  summarise(
    mean_froh = mean(FROH),                           # Mean FROH
    median_froh = median(FROH),                       # Median FROH
    sd_froh = sd(FROH),                               # Standard deviation
    se_froh = sd(FROH)/sqrt(n()),                     # Standard error
    min_froh = min(FROH),                             # Minimum FROH
    max_froh = max(FROH),                             # Maximum FROH
    sample_size = n(),                                # Sample size
    ci_lower = mean_froh - qt(0.975, sample_size-1) * se_froh,  # 95% CI lower bound
    ci_upper = mean_froh + qt(0.975, sample_size-1) * se_froh,  # 95% CI upper bound
    .groups = 'drop'
  )

# FROH statistics by country
froh_country_stats_no_filter <- total_roh_individual_no_filter %>%
  group_by(Country) %>%
  summarise(
    mean_froh = mean(FROH),                           # Mean FROH
    median_froh = median(FROH),                       # Median FROH
    sd_froh = sd(FROH),                               # Standard deviation
    se_froh = sd(FROH)/sqrt(n()),                     # Standard error
    min_froh = min(FROH),                             # Minimum FROH
    max_froh = max(FROH),                             # Maximum FROH
    sample_size = n(),                                # Sample size
    ci_lower = mean_froh - qt(0.975, sample_size-1) * se_froh,  # 95% CI lower bound
    ci_upper = mean_froh + qt(0.975, sample_size-1) * se_froh,  # 95% CI upper bound
    .groups = 'drop'
  )

# ========================================================================
# 5. STATISTICAL TESTING FOR FROH DIFFERENCES
# ========================================================================

# Test for FROH differences between countries
country_test_no_filter <- aov(FROH ~ Country, data = total_roh_individual_no_filter)
country_summary_no_filter <- summary(country_test_no_filter)
country_tukey_no_filter <- TukeyHSD(country_test_no_filter)

# Test for FROH differences between areas
area_test_no_filter <- aov(FROH ~ area, data = total_roh_individual_no_filter)
area_summary_no_filter <- summary(area_test_no_filter)
area_tukey_no_filter <- TukeyHSD(area_test_no_filter)

# Non-parametric test (Kruskal-Wallis) for areas
kruskal_area_no_filter <- kruskal.test(FROH ~ area, data = total_roh_individual_no_filter)

# Pairwise Wilcoxon tests with Bonferroni correction
pairwise_area_no_filter <- pairwise.wilcox.test(total_roh_individual_no_filter$FROH, 
                                                total_roh_individual_no_filter$area, 
                                                p.adjust.method = "bonferroni")

# For ROH segment counts by size class
# Small ROH segments (0.5-1 Mb)
t_test_small <- t.test(roh_count ~ Country.x, 
                       data = subset(plink_roh_size_counts, size_class == "0.5-1 Mb"))

# Medium ROH segments (1-5 Mb)
t_test_medium <- t.test(roh_count ~ Country.x, 
                        data = subset(plink_roh_size_counts, size_class == "1-5 Mb"))

# For total ROH counts
total_roh_by_individual <- plink_roh_size_counts %>%
  group_by(FID, Country.x) %>%
  summarize(total_roh = sum(roh_count), .groups = 'drop')

# Then run the t-test
t_test_total_roh <- t.test(total_roh ~ Country.x, data = total_roh_by_individual)

# For FROH comparison
t_test_froh <- t.test(FROH ~ Country, data = total_roh_individual_no_filter)

# ========================================================================
# 6. PLOTTING
# ========================================================================

# Plot 1: ROH counts by country and size class
roh_count_country_no_filter <- ggplot(roh_country_size_stats_no_filter, 
                                      aes(x = Country, y = mean_count, fill = size_class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = size_class_colors) +
  theme_minimal(base_size = 15) + 
  labs(x = NULL,
       y = "Average number of ROH segments",
       fill = "Size Class",
       title = NULL) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
  
# Plot 2: ROH counts by area and size class
roh_count_area_no_filter <- ggplot(roh_area_size_stats_no_filter, 
                                   aes(x = area, y = mean_count, fill = size_class)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = size_class_colors) +
  theme_minimal(base_size = 15) + 
  labs(x = NULL,
       y = "Average number of ROH segments",
       fill = "Size Class",
       title = NULL) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

# Plot 3: FROH by country
froh_country_plot_no_filter <- ggplot(total_roh_individual_no_filter, 
                                      aes(x = Country, y = FROH, fill = Country)) +
  geom_violin(alpha = 0.7) +
  geom_point(aes(color = Country), alpha = 0.8, size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 95, color = "black", alpha = 0.5) +
  scale_fill_manual(values = country_colors) +
  scale_color_manual(values = country_colors) +
  labs(x = NULL, y = "F(ROH)", title = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

# Plot 4: FROH by area with statistical annotation
froh_area_plot_no_filter <- ggplot(total_roh_individual_no_filter, 
                                   aes(x = area, y = FROH, fill = area)) +
  geom_violin(alpha = 0.7) +
  geom_point(aes(color = area), alpha = 0.8, size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 95, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("SH_rural" = "#009E73", 
                               "SH_urban" = "#c861a5", 
                               "USA" = "#3690c0")) +
  scale_color_manual(values = c("SH_rural" = "#009E73", 
                                "SH_urban" = "#c861a5", 
                                "USA" = "#3690c0")) +
  labs(x = NULL, y = "F(ROH)", title = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

# Add statistical comparison brackets to FROH by area plot
froh_area_plot_stats_no_filter <- froh_area_plot_no_filter +
  stat_compare_means(method = "anova", label.y = max(total_roh_individual_no_filter$FROH) * 1.05) +
  stat_compare_means(comparisons = list(c("SH_rural", "SH_urban"), 
                                        c("SH_urban", "USA"),
                                        c("SH_rural", "USA")), 
                     method = "wilcox.test", 
                     label = "p.signif")

# ========================================================================
# 7. SAVE PLOTS AND RESULTS
# ========================================================================

# Save plots
ggsave("No_filter_ROH_count_by_country.pdf", roh_count_country_no_filter, width = 8, height = 6, dpi = 300)
ggsave("No_filter_ROH_count_by_country.png", roh_count_country_no_filter, width = 6, height = 6, dpi = 300)

ggsave("No_filter_ROH_count_by_area.pdf", roh_count_area_no_filter, width = 8, height = 6, dpi = 300)
ggsave("No_filter_ROH_count_by_area.png", roh_count_area_no_filter, width = 8, height = 6, dpi = 300)

ggsave("No_filter_FROH_by_country.pdf", froh_country_plot_no_filter, width = 8, height = 6, dpi = 300)
ggsave("No_filter_FROH_by_country.png", froh_country_plot_no_filter, width = 6, height = 6, dpi = 300)

ggsave("No_filter_FROH_by_area_with_stats.pdf", froh_area_plot_stats_no_filter, width = 8, height = 6, dpi = 300)
ggsave("No_filter_FROH_by_area_with_stats.png", froh_area_plot_stats_no_filter, width = 8, height = 6, dpi = 300)

# ========================================================================
# 8. SAVE SUMMARY STATISTICS TO CSV FILES
# ========================================================================

# Save ROH statistics
write_csv(roh_country_size_stats_no_filter, "No_filter_roh_stats_by_country_size.csv")
write_csv(roh_area_size_stats_no_filter, "No_filter_roh_stats_by_area_size.csv")
write_csv(froh_country_stats_no_filter, "No_filter_froh_stats_by_country.csv")
write_csv(froh_area_stats_no_filter, "No_filter_froh_stats_by_area.csv")

# Save all individual-level data
write_csv(total_roh_individual_no_filter, "No_filter_all_individual_roh_froh.csv")

# Save statistical test results
capture.output(
  cat("# FROH Differences Between Countries\n\n"),
  print(country_summary_no_filter),
  cat("\n\n# Tukey HSD for Countries\n\n"),
  print(country_tukey_no_filter),
  cat("\n\n# FROH Differences Between Areas\n\n"),
  print(area_summary_no_filter),
  cat("\n\n# Tukey HSD for Areas\n\n"),
  print(area_tukey_no_filter),
  cat("\n\n# Kruskal-Wallis Test for Areas\n\n"),
  print(kruskal_area_no_filter),
  cat("\n\n# Pairwise Wilcoxon Tests with Bonferroni Correction\n\n"),
  print(pairwise_area_no_filter),
  file = "No_filter_froh_statistical_tests.txt"
)

# ========================================================================
# 9. PRINT SUMMARY STATISTICS TO CONSOLE
# ========================================================================

# Print key summary statistics
cat("\n============ NO FILTER: FROH SUMMARY STATISTICS ============\n")
cat("\nBy Country:\n")
print(froh_country_stats_no_filter)

cat("\nBy Area:\n")
print(froh_area_stats_no_filter)

cat("\n============ NO FILTER: STATISTICAL TESTS ============\n")
cat("\nANOVA for Country Differences:\n")
print(country_summary_no_filter)

cat("\nANOVA for Area Differences:\n")
print(area_summary_no_filter)

cat("\nKruskal-Wallis for Area Differences:\n")
print(kruskal_area_no_filter)

cat("\nPairwise Wilcoxon Tests with Bonferroni Correction:\n")
print(pairwise_area_no_filter)

# ========================================================================
# 10. ADDITIONAL ANALYSES - ROH LENGTH DISTRIBUTION
# ========================================================================

# Calculate statistics for ROH lengths by size class
roh_length_stats_no_filter <- plink_roh_size_no_filter %>%
  group_by(size_class) %>%
  summarise(
    mean_length = mean(KB),
    median_length = median(KB),
    sd_length = sd(KB),
    se_length = sd(KB)/sqrt(n()),
    min_length = min(KB),
    max_length = max(KB),
    total_count = n(),
    .groups = 'drop'
  )

# Plot ROH length distribution by size class and country
roh_length_dist_no_filter <- ggplot(plink_roh_size_no_filter %>% 
                                      left_join(popdata, by = c("FID" = "SLF_ID")), 
                                    aes(x = KB, fill = Country)) +
  geom_histogram(position = "dodge", binwidth = 500, alpha = 0.7) +
  scale_fill_manual(values = country_colors) +
  facet_wrap(~size_class, scales = "free_y") +
  theme_minimal() +
  labs(x = "ROH Length (Kb)", y = "Count", 
       title = "Distribution of ROH Lengths by Size Class and Country") +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot and statistics
ggsave("No_filter_ROH_length_distribution.pdf", roh_length_dist_no_filter, width = 10, height = 8, dpi = 300)
write_csv(roh_length_stats_no_filter, "No_filter_roh_length_statistics.csv")

# ========================================================================
# 11. LOCATION-SPECIFIC ANALYSES (DETAILED LOCATIONS)
# ========================================================================

# Calculate ROH statistics by specific location and size class
roh_location_size_stats_no_filter <- roh_individual_no_filter %>%
  group_by(Location, size_class) %>%
  summarise(
    mean_count = mean(roh_count),                      # Mean count
    median_count = median(roh_count),                  # Median count
    sd_count = sd(roh_count),                          # Standard deviation
    se_count = sd(roh_count)/sqrt(n()),                # Standard error
    min_count = min(roh_count),                        # Minimum count
    max_count = max(roh_count),                        # Maximum count
    sample_size = n(),                                 # Sample size
    mean_length = mean(total_length),                  # Mean total length
    median_length = median(total_length),              # Median total length
    sd_length = sd(total_length),                      # Standard deviation of length
    se_length = sd(total_length)/sqrt(n()),            # Standard error of length
    min_length = min(total_length),                    # Minimum total length
    max_length = max(total_length),                    # Maximum total length
    ci_lower_count = mean_count - qt(0.975, sample_size-1) * se_count,  # 95% CI lower bound
    ci_upper_count = mean_count + qt(0.975, sample_size-1) * se_count,  # 95% CI upper bound
    .groups = 'drop'
  )

# FROH statistics by detailed location
froh_location_stats_no_filter <- total_roh_individual_no_filter %>%
  group_by(Location) %>%
  summarise(
    mean_froh = mean(FROH),                           # Mean FROH
    median_froh = median(FROH),                       # Median FROH
    sd_froh = sd(FROH),                               # Standard deviation
    se_froh = sd(FROH)/sqrt(n()),                     # Standard error
    min_froh = min(FROH),                             # Minimum FROH
    max_froh = max(FROH),                             # Maximum FROH
    sample_size = n(),                                # Sample size
    ci_lower = mean_froh - qt(0.975, sample_size-1) * se_froh,  # 95% CI lower bound
    ci_upper = mean_froh + qt(0.975, sample_size-1) * se_froh,  # 95% CI upper bound
    .groups = 'drop'
  )

# Statistical test for FROH differences between specific locations
# Only run if there are sufficient samples per location
if(min(froh_location_stats_no_filter$sample_size) >= 3) {
  location_test_no_filter <- aov(FROH ~ Location, data = total_roh_individual_no_filter)
  location_summary_no_filter <- summary(location_test_no_filter)
  
  # Only run Tukey if ANOVA is significant
  if(summary(location_test_no_filter)[[1]][["Pr(>F)"]][1] < 0.05) {
    location_tukey_no_filter <- TukeyHSD(location_test_no_filter)
  } else {
    location_tukey_no_filter <- "ANOVA not significant, Tukey test not performed"
  }
  
  # Non-parametric test (Kruskal-Wallis)
  kruskal_location_no_filter <- kruskal.test(FROH ~ Location, data = total_roh_individual_no_filter)
  
  # Pairwise Wilcoxon tests
  # Only perform if Kruskal-Wallis is significant
  if(kruskal_location_no_filter$p.value < 0.05) {
    pairwise_location_no_filter <- pairwise.wilcox.test(
      total_roh_individual_no_filter$FROH, 
      total_roh_individual_no_filter$Location, 
      p.adjust.method = "bonferroni"
    )
  } else {
    pairwise_location_no_filter <- "Kruskal-Wallis not significant, pairwise tests not performed"
  }
} else {
  location_test_no_filter <- "Insufficient samples for ANOVA by location"
  location_summary_no_filter <- "Insufficient samples for ANOVA by location"
  location_tukey_no_filter <- "Insufficient samples for Tukey test by location"
  kruskal_location_no_filter <- "Insufficient samples for Kruskal-Wallis test by location"
  pairwise_location_no_filter <- "Insufficient samples for pairwise tests by location"
}

# ========================================================================
# 12. LOCATION-SPECIFIC PLOTS
# ========================================================================

# Order locations for plotting
location_order <- total_roh_individual_no_filter %>%
  group_by(Location) %>%
  summarise(avg_lat = mean(lat, na.rm = TRUE)) %>%
  arrange(avg_lat) %>%
  pull(Location)

# Force factor order of locations
roh_location_size_stats_no_filter$Location <- factor(roh_location_size_stats_no_filter$Location, 
                                                     levels = location_order)
froh_location_stats_no_filter$Location <- factor(froh_location_stats_no_filter$Location, 
                                                 levels = location_order)
total_roh_individual_no_filter$Location <- factor(total_roh_individual_no_filter$Location, 
                                                  levels = location_order)

# Plot 5: ROH counts by detailed location and size class
# Check if too many locations for a standard bar plot (more than 7)
if(length(unique(total_roh_individual_no_filter$Location)) > 7) {
  # Create a faceted plot for many locations
  roh_count_location_no_filter <- ggplot(roh_location_size_stats_no_filter, 
                                         aes(x = size_class, y = mean_count, fill = size_class)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = ci_lower_count, ymax = ci_upper_count), 
                  width = 0.2, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = size_class_colors) +
    facet_wrap(~Location, ncol = 3) +
    theme_minimal(base_size = 15) + 
    labs(x = "Size Class",
         y = "Average number of ROH segments",
         fill = "Size Class",
         title = "No Filter") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 12))
} else {
  # Standard bar plot for fewer locations
  roh_count_location_no_filter <- ggplot(roh_location_size_stats_no_filter, 
                                         aes(x = Location, y = mean_count, fill = size_class)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = ci_lower_count, ymax = ci_upper_count), 
                  width = 0.2, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = size_class_colors) +
    theme_minimal(base_size = 15) + 
    labs(x = NULL,
         y = "Average number of ROH segments",
         fill = "Size Class",
         title = "No Filter") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 12))
}

# Plot 6: FROH by detailed location - Violin plot with SE only
froh_location_plot_no_filter <- ggplot(total_roh_individual_no_filter, 
                                       aes(x = Location, y = FROH, fill = Location)) +
  geom_violin(aes(fill = Location), alpha = 0.5) +
  geom_point(aes(color = Location), alpha = 0.8, size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, alpha = 0.5, color = "black") +  # Mean ± SE error bar
  stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 95, alpha = 0.5) +  # Mean line
  scale_fill_manual(values = location_colors) +
  scale_color_manual(values = location_colors) +
  labs(x = NULL, y = "F(ROH)", title = "No Filter") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

# ========================================================================
# 13. SAVE LOCATION-SPECIFIC RESULTS
# ========================================================================

# Save location-specific plots
ggsave("No_filter_ROH_count_by_location.pdf", roh_count_location_no_filter, width = 12, height = 10, dpi = 300)
ggsave("No_filter_ROH_count_by_location.png", roh_count_location_no_filter, width = 12, height = 10, dpi = 300)

ggsave("No_filter_FROH_by_location.pdf", froh_location_plot_no_filter, width = 10, height = 8, dpi = 300)
ggsave("No_filter_FROH_by_location.png", froh_location_plot_no_filter, width = 10, height = 8, dpi = 300)

# Save location-specific statistics
write_csv(roh_location_size_stats_no_filter, "No_filter_roh_stats_by_location_size.csv")
write_csv(froh_location_stats_no_filter, "No_filter_froh_stats_by_location.csv")

# Save statistical test results for locations
capture.output(
  cat("# FROH Differences Between Detailed Locations\n\n"),
  print(location_summary_no_filter),
  cat("\n\n# Tukey HSD for Locations\n\n"),
  print(location_tukey_no_filter),
  cat("\n\n# Kruskal-Wallis Test for Locations\n\n"),
  print(kruskal_location_no_filter),
  cat("\n\n# Pairwise Wilcoxon Tests for Locations\n\n"),
  print(pairwise_location_no_filter),
  file = "No_filter_froh_location_statistical_tests.txt"
)

# ========================================================================
# 14. PRINT LOCATION-SPECIFIC SUMMARY STATISTICS
# ========================================================================

cat("\n============ NO FILTER: FROH STATISTICS BY LOCATION ============\n")
print(froh_location_stats_no_filter)

cat("\n============ NO FILTER: STATISTICAL TESTS BY LOCATION ============\n")
if(!is.character(location_summary_no_filter)) {
  cat("\nANOVA for Location Differences:\n")
  print(location_summary_no_filter)
} else {
  cat("\nANOVA for Location Differences:", location_summary_no_filter, "\n")
}

if(!is.character(kruskal_location_no_filter)) {
  cat("\nKruskal-Wallis for Location Differences:\n")
  print(kruskal_location_no_filter)
} else {
  cat("\nKruskal-Wallis for Location Differences:", kruskal_location_no_filter, "\n")
}

# ========================================================================
# 15. COMPARISON BETWEEN FILTERED AND NON-FILTERED RESULTS
# ========================================================================

# Read the filtered ROH data for comparison
plink_roh <- read.delim('Combined_RoH.hom', sep='')

# Process filtered ROH data for comparison
plink_roh_size <- plink_roh %>%
  mutate(size_class = case_when(
    KB >= 500 & KB < 1000 ~ "0.5-1 Mb",
    KB >= 1000 & KB < 5000 ~ "1-5 Mb",
    KB >= 5000 ~ ">5 Mb"
  ))

# Calculate individual FROH for filtered data
total_roh_individual <- plink_roh_size %>%
  group_by(FID) %>%
  summarise(
    total_roh_count = n(),
    total_length = sum(KB),
    mean_length = mean(KB),
    min_length = min(KB),
    max_length = max(KB),
    .groups = 'drop'
  ) %>%
  mutate(FROH = total_length / LAUTO) %>%
  left_join(popdata, by = c("FID" = "SLF_ID"))

# Create a combined dataset for comparison
comparison_df <- bind_rows(
  mutate(total_roh_individual, Filter_Status = "Filtered"),
  mutate(total_roh_individual_no_filter, Filter_Status = "No Filter")
)

# Plot comparison of FROH by country with filter status
froh_comparison_country <- ggplot(comparison_df, 
                                  aes(x = Country, y = FROH, fill = Filter_Status)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Filtered" = "#56B4E9", "No Filter" = "#D55E00")) +
  labs(x = NULL, y = "F(ROH)", title = "Comparison of FROH by Country: Filtered vs. No Filter") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "top",
        legend.title = element_text(size = 12))

# Plot comparison of FROH by area with filter status
froh_comparison_area <- ggplot(comparison_df, 
                               aes(x = area, y = FROH, fill = Filter_Status)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Filtered" = "#56B4E9", "No Filter" = "#D55E00")) +
  labs(x = NULL, y = "F(ROH)", title = "Comparison of FROH by Area: Filtered vs. No Filter") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "top",
        legend.title = element_text(size = 12))

# Calculate summary statistics for both filtered and non-filtered data
comparison_stats <- comparison_df %>%
  group_by(Filter_Status, Country) %>%
  summarise(
    mean_froh = mean(FROH),
    median_froh = median(FROH),
    sd_froh = sd(FROH),
    se_froh = sd(FROH)/sqrt(n()),
    min_froh = min(FROH),
    max_froh = max(FROH),
    sample_size = n(),
    .groups = 'drop'
  )

# Calculate difference between filtered and non-filtered
comparison_diff <- comparison_stats %>%
  pivot_wider(
    id_cols = Country,
    names_from = Filter_Status,
    values_from = c(mean_froh, median_froh, sd_froh, se_froh, min_froh, max_froh)
  ) %>%
  mutate(
    diff_mean = mean_froh_No_Filter - mean_froh_Filtered,
    pct_diff_mean = (diff_mean / mean_froh_Filtered) * 100,
    diff_median = median_froh_No_Filter - median_froh_Filtered,
    pct_diff_median = (diff_median / median_froh_Filtered) * 100
  )

# Save comparison plots and statistics
ggsave("FROH_comparison_by_country.pdf", froh_comparison_country, width = 10, height = 8, dpi = 300)
ggsave("FROH_comparison_by_country.png", froh_comparison_country, width = 10, height = 8, dpi = 300)

ggsave("FROH_comparison_by_area.pdf", froh_comparison_area, width = 10, height = 8, dpi = 300)
ggsave("FROH_comparison_by_area.png", froh_comparison_area, width = 10, height = 8, dpi = 300)

write_csv(comparison_stats, "FROH_comparison_statistics.csv")
write_csv(comparison_diff, "FROH_comparison_differences.csv")

# Perform t-tests to compare FROH between filtered and non-filtered approaches
# For each country
for (country_val in unique(comparison_df$Country)) {
  subset_data <- comparison_df %>% filter(Country == country_val)
  t_result <- t.test(FROH ~ Filter_Status, data = subset_data)
  
  cat("\n=== T-test for", country_val, "===\n")
  print(t_result)
}

# For each area
for (area_val in unique(comparison_df$area)) {
  subset_data <- comparison_df %>% filter(area == area_val)
  t_result <- t.test(FROH ~ Filter_Status, data = subset_data)
  
  cat("\n=== T-test for", area_val, "===\n")
  print(t_result)
}

# Save t-test results
capture.output(
  cat("# T-tests comparing FROH between Filtered and Non-Filtered data\n\n"),
  cat("## By Country\n"),
  for (country_val in unique(comparison_df$Country)) {
    subset_data <- comparison_df %>% filter(Country == country_val)
    t_result <- t.test(FROH ~ Filter_Status, data = subset_data)
    
    cat("\n=== T-test for", country_val, "===\n")
    print(t_result)
  },
  cat("\n\n## By Area\n"),
  for (area_val in unique(comparison_df$area)) {
    subset_data <- comparison_df %>% filter(area == area_val)
    t_result <- t.test(FROH ~ Filter_Status, data = subset_data)
    
    cat("\n=== T-test for", area_val, "===\n")
    print(t_result)
  },
  file = "FROH_comparison_t_tests.txt"
)



# ========================================================================
# 16. INFER COALESCENCE TIME (TMRCA) BY ROH LENGTH CLASS
# ========================================================================
library(tibble)

# Add Country column to full ROH dataset (if not already joined)
plink_roh_size_country <- plink_roh_size_no_filter %>%
  left_join(popdata, by = c("FID" = "SLF_ID"))

# Convert KB to Mb
plink_roh_size_country <- plink_roh_size_country %>%
  mutate(Mb = KB / 1000)

# Compute median ROH size by country and size class
roh_medians_by_country <- plink_roh_size_country %>%
  group_by(Country, size_class) %>%
  summarise(
    median_Mb = median(Mb, na.rm = TRUE),
    .groups = 'drop'
  )


# ========================================================================
# 18. TMRCA Estimates with SLF-specific Recombination Rates (1, 2, 3 cM/Mb)
# ========================================================================
# Define updated recombination rates (SLF-specific range)
slf_recomb_rates <- c(1.0)

# Recalculate TMRCA using country-specific ROH median sizes
tmrca_slf <- roh_medians_by_country %>%
  crossing(recomb_rate_cM_per_Mb = slf_recomb_rates) %>%
  mutate(
    L_cM = median_Mb * recomb_rate_cM_per_Mb,
    TMRCA_generations = 100 / (2 * L_cM),
    TMRCA_years = TMRCA_generations * 1  # 1 generation per year
  ) %>%
  mutate(
    recomb_label = paste0(recomb_rate_cM_per_Mb, " cM/Mb"),
    TMRCA_generations = round(TMRCA_generations, 1),
    TMRCA_years = round(TMRCA_years, 1)
  ) %>%
  select(Country, size_class, median_Mb, recomb_label, TMRCA_generations, TMRCA_years)

# Save output
write_csv(tmrca_slf, "TMRCA_SLF_custom_recomb_range.csv")

# Print results
cat("\n============ TMRCA Estimates for SLF (1–3 cM/Mb) by ROH Class ============\n")
print(tmrca_slf)
















#A. Genome-wide Distribution of ROHs
# First, let's examine the genomic distribution of ROHs to identify potential ROH islands or hotspots:
# Add chromosome and position information to ROH data
# Fix the genomic distribution analysis
roh_genomic_distribution <- plink_roh_size_no_filter %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  group_by(CHR) %>%
  # Use middle position between POS1 and POS2 for binning
  mutate(POS = (POS1 + POS2)/2,
         position_bin = cut(POS, breaks = seq(min(POS1), max(POS2), length.out = 50), 
                            include.lowest = TRUE)) %>%
  group_by(CHR, position_bin, Country) %>%
  summarise(
    roh_count = n(),
    mean_kb = mean(KB),
    .groups = 'drop'
  )

# Plot ROH density across chromosomes
roh_genomic_distribution <- ggplot(roh_genomic_distribution, aes(x = position_bin, y = roh_count, fill = Country)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~CHR, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Genomic Distribution of ROH",
       x = "Chromosome Position",
       y = "Number of ROH Segments")

# Save to PDF
ggsave("roh_genomic_distribution.pdf", roh_genomic_distribution, width = 12, height = 8, dpi = 300)



# Step 1: Filter for ROH segments > 5 Mb
roh_gt5mb <- plink_roh_size_no_filter %>%
  filter(KB > 5000)

roh_genomic_distribution <- roh_gt5mb %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  group_by(CHR) %>%
  # Use middle position between POS1 and POS2 for binning
  mutate(POS = (POS1 + POS2)/2,
         position_bin = cut(POS, breaks = seq(min(POS1), max(POS2), length.out = 50), 
                            include.lowest = TRUE)) %>%
  group_by(CHR, position_bin, Country) %>%
  summarise(
    roh_count = n(),
    mean_kb = mean(KB),
    .groups = 'drop'
  )


library(dplyr)
library(forcats)

roh_genomic_distribution_gt5mb <- roh_genomic_distribution_gt5mb %>%
  mutate(
    bin_start = bin_5Mb,  # keep as numeric
    bin_label = paste0(bin_start / 1e6, "Mb")
  ) %>%
  group_by(CHR) %>%
  mutate(
    bin_label = fct_reorder(bin_label, bin_start, .fun = min)  # sort x-axis correctly within each facet
  ) %>%
  ungroup()

roh_genomic_distribution_gt5mb <- ggplot(roh_genomic_distribution_gt5mb, aes(x = bin_label, y = roh_count, fill = Country)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.9) +
  facet_wrap(~CHR, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Distribution of Long ROHs (>5 Mb)",
       x = "Genomic Bin (Mb)",
       y = "Number of ROH Segments")

ggsave("Long_ROH_Distribution_5MbBins.pdf", roh_genomic_distribution_gt5mb, width = 12, height = 10, dpi = 300)

# Load the USA selection candidates: 184 SNPs
selection_candidates <- read.csv("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/185overlap_prda_pcadapt_full_info.csv")

library(GenomicRanges)

# Subset your plink_roh_size_no_filter to long ROHs
roh_long <- plink_roh_size_no_filter %>%
  filter(KB > 5000)

# Convert to GRanges
roh_gr <- GRanges(
  seqnames = roh_long$CHR,
  ranges = IRanges(start = roh_long$POS1, end = roh_long$POS2),
  FID = roh_long$FID
)

# Convert SNPs to GRanges
snp_gr <- GRanges(
  seqnames = selection_candidates$CHROM,
  ranges = IRanges(start = selection_candidates$POS, width = 1),
  SNP_ID = selection_candidates$SNP
)

# Find SNPs that fall within any long ROH
overlaps <- findOverlaps(snp_gr, roh_gr)

# Extract overlapping SNPs
snp_in_roh <- selection_candidates[queryHits(overlaps), ]


# Count how many overlap
nrow(snp_in_roh)  # How many of the 185 fall in long ROHs

# Optional: Add a new column to flag SNPs in long ROHs
selection_candidates$in_long_roh <- selection_candidates$SNP %in% snp_in_roh$SNP
table(selection_candidates$in_long_roh)

write.csv(snp_in_roh, "Selection_SNPs_in_Long_ROHs.csv", row.names = FALSE)

length(unique(snp_in_roh$SNP))  # should return 46

snp_in_roh_unique <- snp_in_roh %>%
  distinct(SNP, .keep_all = TRUE)

write.csv(snp_in_roh_unique, "Unique_Selection_SNPs_in_Long_ROHs.csv", row.names = FALSE)


# Assign each SNP to its 5 Mb bin
snp_in_long_roh_binned <- selection_candidates %>%
  filter(in_long_roh) %>%
  mutate(
    bin_start = floor(POS / 5e6) * 5e6,
    bin_label = paste0(bin_start / 1e6, "Mb"),
    CHR = as.character(CHROM)  # match to ROH data's CHR
  )

# Ensure bin_label is ordered the same way
snp_in_long_roh_binned <- snp_in_long_roh_binned %>%
  group_by(CHR) %>%
  mutate(bin_label = forcats::fct_reorder(bin_label, bin_start)) %>%
  ungroup()

# Recreate the data frame (you can copy from earlier code)
roh_genomic_distribution_gt5mb <- roh_gt5mb %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  mutate(POS = (POS1 + POS2) / 2,
         bin_5Mb = floor(POS / 5e6) * 5e6) %>%
  group_by(CHR, bin_5Mb, Country) %>%
  summarise(
    roh_count = n(),
    mean_kb = mean(KB),
    .groups = 'drop'
  ) %>%
  mutate(
    bin_start = bin_5Mb,
    bin_label = paste0(bin_start / 1e6, "Mb")
  ) %>%
  group_by(CHR) %>%
  mutate(bin_label = forcats::fct_reorder(bin_label, bin_start)) %>%
  ungroup()

# Step 1: Get all bins across both datasets
all_bin_levels <- unique(c(roh_genomic_distribution_gt5mb$bin_label,
                           snp_in_long_roh_binned$bin_label))

# Step 2: Apply the same bin factor levels to both
roh_genomic_distribution_gt5mb$bin_label <- factor(roh_genomic_distribution_gt5mb$bin_label, levels = all_bin_levels)
snp_in_long_roh_binned$bin_label <- factor(snp_in_long_roh_binned$bin_label, levels = all_bin_levels)


# Get all combinations of CHR and bin_start from SNPs and ROHs
all_bins <- expand.grid(
  CHR = unique(c(roh_genomic_distribution_gt5mb$CHR, snp_in_long_roh_binned$CHR)),
  bin_5Mb = unique(c(roh_genomic_distribution_gt5mb$bin_start, snp_in_long_roh_binned$bin_start))
)

# Join with ROH counts, fill NAs with 0
roh_genomic_distribution_full <- all_bins %>%
  left_join(roh_genomic_distribution_gt5mb, by = c("CHR", "bin_5Mb" = "bin_start")) %>%
  mutate(
    roh_count = ifelse(is.na(roh_count), 0, roh_count),
    bin_label = paste0(bin_5Mb / 1e6, "Mb")
  ) %>%
  group_by(CHR) %>%
  mutate(bin_label = fct_reorder(bin_label, bin_5Mb)) %>%
  ungroup()

roh_genomic_distribution_full <- all_bins %>%
  left_join(roh_genomic_distribution_gt5mb, by = c("CHR", "bin_5Mb" = "bin_start")) %>%
  mutate(
    roh_count = ifelse(is.na(roh_count), 0, roh_count),
    mean_kb = ifelse(is.na(mean_kb), 0, mean_kb),
    Country = ifelse(is.na(Country), "USA", Country),  # Assuming all your data is USA
    bin_label = paste0(bin_5Mb / 1e6, "Mb")
  ) %>%
  group_by(CHR) %>%
  mutate(bin_label = fct_reorder(bin_label, bin_5Mb)) %>%
  ungroup()


ggplot(roh_genomic_distribution_full, aes(x = bin_label, y = roh_count, fill = Country)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.9) +
  # Highlight overlapping selection SNPs
  geom_point(data = snp_in_long_roh_binned,
             aes(x = bin_label, y = -0.3),  # offset point to appear below bars
             color = "red", size = 1.5, inherit.aes = FALSE) +
  facet_wrap(~CHR, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Distribution of Long ROHs (>5 Mb) with Candidate SNPs",
       x = "Genomic Bin (Mb)",
       y = "Number of ROH Segments")


library(ggplot2)
library(forcats)
library(dplyr)

# Apply consistent factor levels for bins
all_bin_levels <- unique(c(roh_genomic_distribution_full$bin_label,
                           snp_in_long_roh_binned$bin_label))

roh_genomic_distribution_full$bin_label <- factor(roh_genomic_distribution_full$bin_label,
                                                  levels = all_bin_levels)
snp_in_long_roh_binned$bin_label <- factor(snp_in_long_roh_binned$bin_label,
                                           levels = all_bin_levels)

# Plot: ROH bar height + red SNP dots
ggplot(roh_genomic_distribution_full, aes(x = bin_label, y = roh_count, fill = Country)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.9) +
  geom_point(data = snp_in_long_roh_binned,
             aes(x = bin_label, y = -0.3),
             color = "red", size = 1.5, inherit.aes = FALSE) +
  facet_wrap(~CHR, scales = "free_x") +
  scale_fill_manual(values = c("USA" = "#3690c0")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.position = "none"
  ) +
  labs(
    title = "Distribution of Long ROHs (>5 Mb) with Selection SNPs",
    x = "Genomic Bin (Mb)",
    y = "Number of ROH Segments"
  )

ggsave("Long_ROH_Distribution_with_SNPs_USAcolor.pdf", width = 12, height = 10, dpi = 300)
ggsave("Long_ROH_Distribution_with_SNPs_USAcolor.png", width = 12, height = 10, dpi = 300)


# Pull SNPs in PRKG1, HIVEP1, CLINT1
snp_targets <- snp_in_long_roh_binned %>%
  filter(Gene_Name %in% c("PRKG1", "HIVEP1", "CLINT1")) %>%
  select(SNP, CHR = CHROM, POS, Gene_Name)

# Join to ROHs where the SNP falls within the ROH's POS1–POS2 range on the same chromosome
snp_roh_match <- snp_targets %>%
  left_join(roh_gt5mb, by = "CHR", relationship = "many-to-many") %>%
  filter(POS >= POS1, POS <= POS2) %>%
  select(SNP, Gene_Name, CHR, POS, FID, POS1, POS2, KB)








#compute ROH metrics per chromosome (chr_1 to chr_13) for USA samples only, including normalization by chromosome length
# ===============================
# ROH PER CHROMOSOME (USA ONLY)
# ===============================

# Step 1: Define chromosome lengths in Mb
chrom_lengths <- tibble::tibble(
  CHR = paste0("chr_", 1:13),
  length_Mb = c(
    404.011691, 238.731119, 227.282743, 217.703105, 182.426808,
    165.984238, 149.132726, 142.605837, 133.192330, 123.181216,
    83.140728, 78.191594, 55.466916
  )
)

# Step 2: Filter to USA samples and chromosomes 1–13
roh_usa_chr <- plink_roh_size_no_filter %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  filter(Country == "USA", CHR %in% paste0("chr_", 1:13))

# Step 3: Summarize ROH count and length per chromosome per individual
roh_chr_summary <- roh_usa_chr %>%
  group_by(FID, CHR) %>%
  summarise(
    total_roh_count = n(),
    total_roh_length_Kb = sum(KB),
    .groups = 'drop'
  )

# Step 4: Merge with chromosome lengths and calculate normalized metrics
roh_chr_summary <- roh_chr_summary %>%
  left_join(chrom_lengths, by = "CHR") %>%
  mutate(
    total_roh_length_Mb = total_roh_length_Kb / 1000,
    FROH_chr = total_roh_length_Mb / length_Mb,
    roh_density = total_roh_count / length_Mb
  )

# Step 5: Summarize across individuals per chromosome
roh_chr_stats <- roh_chr_summary %>%
  group_by(CHR) %>%
  summarise(
    mean_count = mean(total_roh_count),
    sd_count = sd(total_roh_count),
    mean_length_Mb = mean(total_roh_length_Mb),
    sd_length_Mb = sd(total_roh_length_Mb),
    mean_FROH_chr = mean(FROH_chr),
    sd_FROH_chr = sd(FROH_chr),
    mean_density = mean(roh_density),
    .groups = 'drop'
  )

# View or save results
print(roh_chr_stats)
write_csv(roh_chr_stats, "USA_ROH_per_chromosome_normalized_chr1to13.csv")



# ===============================
# Step 0: Define chromosome lengths
# ===============================
chrom_lengths <- tibble::tibble(
  CHR = paste0("chr_", 1:13),
  length_Mb = c(
    404.011691, 238.731119, 227.282743, 217.703105, 182.426808,
    165.984238, 149.132726, 142.605837, 133.192330, 123.181216,
    83.140728, 78.191594, 55.466916
  )
)

# ===============================
# Step 1: Filter to long ROHs (>5 Mb) from USA samples
# ===============================
roh_long_usa_chr <- plink_roh_size_no_filter %>%
  filter(KB > 5000) %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  filter(Country == "USA", CHR %in% paste0("chr_", 1:13))

# ===============================
# Step 2: Summarize by individual and chromosome
# ===============================
roh_long_chr_summary <- roh_long_usa_chr %>%
  group_by(FID, CHR) %>%
  summarise(
    total_roh_count = n(),
    total_roh_length_Kb = sum(KB),
    .groups = 'drop'
  )

# ===============================
# Step 3: Normalize by chromosome length
# ===============================
roh_long_chr_summary <- roh_long_chr_summary %>%
  left_join(chrom_lengths, by = "CHR") %>%
  mutate(
    total_roh_length_Mb = total_roh_length_Kb / 1000,
    FROH_chr = total_roh_length_Mb / length_Mb,
    roh_density = total_roh_count / length_Mb
  )

# ===============================
# Step 4: Summarize across individuals per chromosome
# ===============================
roh_long_chr_stats <- roh_long_chr_summary %>%
  group_by(CHR) %>%
  summarise(
    mean_count = mean(total_roh_count),
    sd_count = sd(total_roh_count),
    mean_length_Mb = mean(total_roh_length_Mb),
    sd_length_Mb = sd(total_roh_length_Mb),
    mean_FROH_chr = mean(FROH_chr),
    sd_FROH_chr = sd(FROH_chr),
    mean_density = mean(roh_density),
    .groups = 'drop'
  )

# ===============================
# Step 5: View or Save
# ===============================
print(roh_long_chr_stats)
write_csv(roh_long_chr_stats, "USA_Long_ROH_per_chromosome_normalized.csv")


# Step 1: Filter to long ROHs (>5 Mb) and USA samples
roh_long_usa_chr <- plink_roh_size_no_filter %>%
  filter(KB > 5000) %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  filter(Country == "USA", CHR %in% paste0("chr_", 1:13))

# Step 2: Summarize by chromosome (not per individual)
roh_long_chr_raw <- roh_long_usa_chr %>%
  group_by(CHR) %>%
  summarise(
    num_long_roh = n(),
    total_length_Mb = sum(KB) / 1000,
    mean_length_Mb = mean(KB) / 1000,
    median_length_Mb = median(KB) / 1000,
    max_length_Mb = max(KB) / 1000,
    .groups = 'drop'
  ) %>%
  arrange(desc(num_long_roh))

# Step 3: View or save
print(roh_long_chr_raw)
write_csv(roh_long_chr_raw, "USA_Long_ROH_raw_by_chromosome.csv")



library(ggplot2)
library(dplyr)

# OPTIONAL: restrict to certain chromosomes (e.g., chr_1, chr_6, chr_8)
# roh_long_usa_chr <- roh_long_usa_chr %>% filter(CHR %in% c("chr_1", "chr_6", "chr_8"))

# Step 1: reorder FID for plotting
roh_long_usa_chr <- roh_long_usa_chr %>%
  mutate(FID = factor(FID, levels = unique(FID[order(CHR, POS1)])))

roh_long_usa_chr %>% filter(CHR == "chr_8") %>% arrange(desc(KB)) %>% head()

roh_position_plot <- ggplot(roh_long_usa_chr, 
                            aes(xmin = POS1, xmax = POS2, ymin = as.numeric(FID) - 0.4, ymax = as.numeric(FID) + 0.4)) +
  geom_rect(fill = "#c861a5", alpha = 0.8) +
  facet_wrap(~CHR, scales = "fixed")  # instead of "free_x"
  scale_x_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(x = "Genomic Position (bp)", y = "Individual", title = "Long ROH (>5Mb) Regions Across Chromosomes") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())


# Step 2: Plot ROH segments along the chromosome
roh_position_plot <- ggplot(roh_long_usa_chr, aes(xmin = POS1, xmax = POS2, y = FID)) +
  geom_segment(aes(x = POS1, xend = POS2, y = FID, yend = FID), color = "#c861a5", size = 1.2, alpha = 0.8) +
  facet_wrap(~CHR, scales = "free_x", ncol = 1) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Long ROH (>5Mb) Regions Across Chromosomes",
    x = "Genomic Position (bp)",
    y = "Individual"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# View or save the plot
print(roh_position_plot)
ggsave("ROH_chr_position_plot.pdf", roh_position_plot, width = 12, height = 10, dpi = 300)


# Filter long ROHs (>5Mb) in USA samples
roh_long_usa <- plink_roh_size_no_filter %>%
  filter(KB > 5000) %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  filter(Country == "USA")

# Count total number of long ROHs in USA
total_long_roh_usa <- nrow(roh_long_usa)
cat("Total long ROHs (>5Mb) in USA:", total_long_roh_usa, "\n")

# Count unique individuals with at least one long ROH
individuals_with_long_roh <- roh_long_usa %>%
  distinct(FID)

num_individuals <- nrow(individuals_with_long_roh)
cat("Number of USA individuals with long ROHs:", num_individuals, "\n")

# Print list of individuals
cat("Individuals with long ROHs:\n")
print(individuals_with_long_roh$FID)

# Optional: Export full list with details of each long ROH
write_csv(roh_long_usa, "USA_Long_ROHs_over_5Mb.csv")



library(dplyr)
library(ggplot2)

# Summarize total length (in Mb) and count of long ROHs (>5 Mb) per individual
roh_summary <- plink_roh_size_no_filter %>%
  left_join(popdata, by = c("FID" = "SLF_ID")) %>%
  group_by(FID, Country) %>%
  summarise(
    total_length_Mb = sum(KB) / 1000,
    roh_count = n(),
    .groups = "drop"
  )

# Clean scatter plot with strong axes
ggplot(roh_summary, aes(x = total_length_Mb, y = roh_count, color = Country)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("USA" = "#3690c0", "CN" = "#c861a5")) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(size = 0.6),
    legend.background = element_blank()
  ) +
  labs(
    x = "Total ROH Length(Mb)",
    y = "Number of ROH Segments",
    color = "Country"
  )

ggsave("scatter_plot_ROH_distribution.pdf", width = 8, height = 6)
ggsave("scatter_plot_ROH_distribution.png", width = 8, height = 6, dpi = 300)

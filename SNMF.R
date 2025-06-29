setwd("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated")

library(LEA)
library(RColorBrewer)

#convert ped to geno file type - only run once
ped2geno('/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Genetic_Distance/Combined.ped')
ped2geno('/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Shanghai.ped')
ped2geno('/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA.ped')

SH_snmf <- snmf("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Shanghai.geno",K = 1:10, 
                 entropy = TRUE, ploidy = 2, project = "new", repetitions = 10)

pdf("SH_snmf_cross_entropy_plot.pdf", width = 8, height = 6)
# plot ancestry suggests 2 ancestral populations
plot(SH_snmf, col = "blue", pch = 19, cex = 1.2) 
dev.off()


#needs entropy=T in snmf Cross-entropy quantifies how well the probability distribution q predicts p. A smaller cross-entropy indicates that q is a better approximation of p.
entropy=sapply(1:10,FUN=function(i)
  cross.entropy(SH_snmf,K = i))

# Compute the mean cross-entropy for each K
mean_cross_entropy <- colMeans(entropy)

# Find the K with the lowest mean cross-entropy
lowest_mean_value <- min(mean_cross_entropy)
best_k <- which.min(mean_cross_entropy)

# Print the results
cat("The K with the lowest mean cross-entropy is:", best_k, "\n")
cat("The lowest mean cross-entropy value is:", lowest_mean_value, "\n")


USA_snmf <- snmf("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA.geno",K = 1:20, 
                entropy = TRUE, ploidy = 2, project = "new")

pdf("K20_USA_snmf_cross_entropy_plot.pdf", width = 8, height = 6)
# plot ancestry suggests 2 ancestral populations
plot(USA_snmf , col = "blue", pch = 19, cex = 1.2) 
dev.off()

#needs entropy=T in snmf Cross-entropy quantifies how well the probability distribution q predicts p. A smaller cross-entropy indicates that q is a better approximation of p.
entropy=sapply(1:10,FUN=function(i)
  cross.entropy(SH_snmf,K = i))

# Compute the mean cross-entropy for each K
mean_cross_entropy <- colMeans(entropy)

# Find the K with the lowest mean cross-entropy
lowest_mean_value <- min(mean_cross_entropy)
best_k <- which.min(mean_cross_entropy)

# Print the results
cat("The K with the lowest mean cross-entropy is:", best_k, "\n")
cat("The lowest mean cross-entropy value is:", lowest_mean_value, "\n")



country_colors <- c(
  "USA" = "#3690c0",
  "CN" = "#c861a5"
)

SH_colors <- c(
  "Jinqiao" = "#3690c0",
  "Sheshan" = "#984ea3",
  "Tianma" = "#4daf4a"
)

USA_colors <- c(
  "Bronx" = "#252525", 
  "Brooklyn" =  "#e6ab02",
  "Camden, NJ" = "#a65628",
  "Connecticut" = "#e41a1c",
  "Governor Island" = "#08589e",
  "Manhattan" = "#c861a5",
  "Staten island" = "#3690c0",
  "West Orange, NJ" = "#009E73"
)


# Load population data
popdata <- read.csv("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_SLF_MS_Info.csv", header=T)   # Population metadata

popbed <- popdata[,c("SLF_ID","Location")]
poplist.names <- popbed[,2]

# plot ancestry suggests 2 ancestral populations
plot(SLF_snmf, col = "blue", pch = 19, cex = 1.2) 

#needs entropy=T in snmf Cross-entropy quantifies how well the probability distribution q predicts p. A smaller cross-entropy indicates that q is a better approximation of p.
entropy=sapply(1:10,FUN=function(i)
  cross.entropy(SLF_snmf,K = i))

# Compute the mean cross-entropy for each K
mean_cross_entropy <- colMeans(entropy)

# Find the K with the lowest mean cross-entropy
lowest_mean_value <- min(mean_cross_entropy)
best_k <- which.min(mean_cross_entropy)

# Print the results
cat("The K with the lowest mean cross-entropy is:", best_k, "\n")
cat("The lowest mean cross-entropy value is:", lowest_mean_value, "\n")

# Create a data frame for visualization
df <- data.frame(
  K = 1:length(mean_cross_entropy),
  Mean_Cross_Entropy = mean_cross_entropy
)

write_csv(df, 'snmf_cv_mean.csv')

ce <-  cross.entropy(SLF_snmf, K =2)
ce

############# Explore Shanghai Data ######


# Load Shanghai-specific .geno file and run SNMF for K = 1 to 10
SH_snmf <- snmf("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Genetic_Distance/Shanghai.geno",
                K = 1:10, 
                entropy = TRUE, 
                ploidy = 2, 
                project = "new", 
                repetitions = 10)

# Plot cross-entropy to choose the optimal K
plot(SH_snmf, col = "blue", pch = 19, cex = 1.2)

# Cross-entropy for each K
entropy <- sapply(1:10, FUN = function(i) cross.entropy(SH_snmf, K = i))
print(entropy)

# Load population data
popdata <- read.csv("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_SLF_MS_Info.csv", header = TRUE)

popdata <- popdata %>% filter(Location %in% c("Jinqiao", "Sheshan", "Tianma"))

# Prepare Shanghai-specific metadata
popbed_SH <- popdata[, c("SLF_ID", "Location")]

# USA sNMF
usa_snmf <- snmf("Genetic_Distance/USA.geno", K = 1:10, entropy = TRUE, ploidy = 2, project = "new", repetitions = 10)
plot(usa_snmf, col = "blue", pch = 19)
pdf("usa_snmf_cross_entropy_plot.pdf", width = 8, height = 6)
plot(usa_snmf, col = "blue", pch = 19)
dev.off()

usa_entropy <- sapply(1:10, function(k) cross.entropy(usa_snmf, K = k))
usa_best_k <- which.min(colMeans(usa_entropy))
usa_best_run_k2 <- which.min(cross.entropy(usa_snmf, K = 2))
usa_best_run_k3 <- which.min(cross.entropy(usa_snmf, K = 3))

# Load USA popdata
popdata <- read.csv("Final_SLF_MS_Info.csv")
popdata_usa <- popdata %>% filter(Country == "USA")
popbed_usa <- popdata_usa[, c("SLF_ID", "Location")]

# ---- K = 2 ----
qmat_usa_k2 <- LEA::Q(usa_snmf, K = 2, run = usa_best_run_k2)
colnames(qmat_usa_k2) <- c("P1", "P2")

qdf_usa_k2 <- qmat_usa_k2 %>%
  as_tibble() %>%
  mutate(individual = popbed_usa$SLF_ID, region = popbed_usa$Location) %>%
  pivot_longer(cols = P1:P2, names_to = "ancestry", values_to = "proportion")

qdf_usa_k2 <- qdf_usa_k2 %>%
  group_by(region, individual) %>%
  mutate(likely_assignment = ancestry[which.max(proportion)],
         assignment_prob = max(proportion)) %>%
  arrange(region, likely_assignment, desc(assignment_prob)) %>%
  ungroup() %>%
  mutate(individual = forcats::fct_inorder(factor(individual)))

region_order_usa <- popdata_usa %>%
  group_by(Location) %>%
  summarize(mean_lat = mean(lat, na.rm = TRUE)) %>%
  arrange(mean_lat) %>%
  pull(Location)

qdf_usa_k2 <- qdf_usa_k2 %>%
  mutate(region = factor(region, levels = region_order_usa))

ancestry_colors_k2 <- c("P1" = "#1f78b4", "P2" = "#33a02c")

ggplot(qdf_usa_k2, aes(x = individual, y = proportion, fill = ancestry)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = ancestry_colors_k2) +
  facet_wrap(~ region, nrow = 1, scales = "free_x", strip.position = "bottom") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12), panel.grid = element_blank()) +
  ylab("Proportion") +
  ggsave("snmf_K2_USA_plot.pdf", width = 16, height = 8)


# Shanghai sNMF
sh_snmf <- snmf("Genetic_Distance/Shanghai.geno", K = 1:10, entropy = TRUE, ploidy = 2, project = "new", repetitions = 10)
plot(sh_snmf, col = "blue", pch = 19)

pdf("SH_snmf_cross_entropy_plot.pdf", width = 8, height = 6)
plot(sh_snmf, col = "blue", pch = 19)
dev.off()

sh_entropy <- sapply(1:10, function(k) cross.entropy(sh_snmf, K = k))
sh_best_k <- which.min(colMeans(sh_entropy))
sh_best_run_k2 <- which.min(cross.entropy(sh_snmf, K = 2))
sh_best_run_k3 <- which.min(cross.entropy(sh_snmf, K = 3))

# Shanghai metadata
popdata <- read.csv("Final_SLF_MS_Info.csv")
popdata_sh <- popdata %>% filter(Location %in% c("Jinqiao", "Sheshan", "Tianma"))
popbed_sh <- popdata_sh[, c("SLF_ID", "Location")]


# Load required libraries
library(LEA)
library(tidyverse)
library(RColorBrewer)

# Set working directory
setwd("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated")

# Load population data
popdata <- read.csv("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_SLF_MS_Info.csv", header = TRUE)

#================================================================================
# USA K=1 ANALYSIS AND PLOTTING
#================================================================================
setwd("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Genetic_Distance")

# Load the existing USA SNMF project (assuming it's already run)
usa_snmf <- load.snmfProject("USA.snmfProject")

# Get cross-entropy for K=1
ce_k1_usa <- cross.entropy(usa_snmf, K = 1)
best_run_k1_usa <- which.min(ce_k1_usa)

cat("USA K=1 - Best run:", best_run_k1_usa, "Cross-entropy:", min(ce_k1_usa), "\n")

# Extract Q matrix for K=1 (will be all 1s since there's only one population)
q_mat_k1_usa <- LEA::Q(usa_snmf, K = 1, run = best_run_k1_usa)
colnames(q_mat_k1_usa) <- "P1"

# Prepare USA metadata
popdata_usa <- popdata %>% filter(Country == "USA")
popbed_usa <- popdata_usa[, c("SLF_ID", "Location")]

# Convert Q matrix to data frame
q_df_k1_usa <- q_mat_k1_usa %>%
  as_tibble() %>%
  mutate(individual = popbed_usa$SLF_ID,
         region = popbed_usa$Location,
         ancestry = "P1",
         proportion = P1) %>%
  select(-P1)

# Order regions by latitude (south to north)
region_order_usa <- popdata_usa %>%
  group_by(Location) %>%
  summarize(mean_lat = mean(lat, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_lat) %>%
  pull(Location)

# Order individuals within regions
q_df_ordered_k1_usa <- q_df_k1_usa %>%
  mutate(region = factor(region, levels = region_order_usa)) %>%
  arrange(region, individual) %>%
  mutate(individual = forcats::fct_inorder(factor(individual)))

# Define color for USA (single population)
usa_color <- "#3690c0"

# Create K=1 plot for USA
plot_k1_usa <- q_df_ordered_k1_usa %>%
  ggplot(aes(x = individual, y = proportion, fill = ancestry)) +
  geom_col(width = 0.9, fill = usa_color) +
  labs(title = "USA Spotted Lanternfly Population Structure (K=1)",
       subtitle = "No population subdivision detected - single panmictic population",
       y = "Ancestry Proportion") +
  theme_minimal() +
  facet_wrap(~ region, nrow = 1, strip.position = "bottom", scales = "free_x") +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "transparent", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
    legend.position = "none"
  ) +
  ylim(0, 1)

# Save USA K=1 plot
ggsave(filename = "snmf_K1_USA_ancestry_plot.pdf", 
       plot = plot_k1_usa, 
       device = "pdf", 
       width = 16, 
       height = 8, 
       units = "in")

ggsave(filename = "snmf_K1_USA_ancestry_plot.png", 
       plot = plot_k1_usa, 
       device = "png", 
       width = 16, 
       height = 8, 
       units = "in",
       dpi = 300)

cat("USA K=1 plot saved: snmf_K1_USA_ancestry_plot.pdf\n")

#================================================================================
# SHANGHAI K=1 ANALYSIS AND PLOTTING
#================================================================================

# Load the existing Shanghai SNMF project
sh_snmf <- load.snmfProject("Shanghai.snmfProject")

# Get cross-entropy for K=1
ce_k1_sh <- cross.entropy(sh_snmf, K = 1)
best_run_k1_sh <- which.min(ce_k1_sh)

cat("Shanghai K=1 - Best run:", best_run_k1_sh, "Cross-entropy:", min(ce_k1_sh), "\n")

# Extract Q matrix for K=1
q_mat_k1_sh <- LEA::Q(sh_snmf, K = 1, run = best_run_k1_sh)
colnames(q_mat_k1_sh) <- "P1"

# Prepare Shanghai metadata
popdata_sh <- popdata %>% filter(Location %in% c("Jinqiao", "Sheshan", "Tianma"))
popbed_sh <- popdata_sh[, c("SLF_ID", "Location")]

# Convert Q matrix to data frame
q_df_k1_sh <- q_mat_k1_sh %>%
  as_tibble() %>%
  mutate(individual = popbed_sh$SLF_ID,
         region = popbed_sh$Location,
         ancestry = "P1",
         proportion = P1) %>%
  select(-P1)

# Order Shanghai regions by latitude
region_order_sh <- popdata_sh %>%
  group_by(Location) %>%
  summarize(mean_lat = mean(lat, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_lat) %>%
  pull(Location)

# Order individuals within regions
q_df_ordered_k1_sh <- q_df_k1_sh %>%
  mutate(region = factor(region, levels = region_order_sh)) %>%
  arrange(region, individual) %>%
  mutate(individual = forcats::fct_inorder(factor(individual)))

# Define color for Shanghai (single population)
sh_color <- "#c861a5"

# Create K=1 plot for Shanghai
plot_k1_sh <- q_df_ordered_k1_sh %>%
  ggplot(aes(x = individual, y = proportion, fill = ancestry)) +
  geom_col(width = 0.9, fill = sh_color) +
  labs(title = "Shanghai Spotted Lanternfly Population Structure (K=1)",
       subtitle = "No population subdivision detected - single panmictic population",
       y = "Ancestry Proportion") +
  theme_minimal() +
  facet_wrap(~ region, nrow = 1, strip.position = "bottom", scales = "free_x") +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "transparent", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
    legend.position = "none"
  ) +
  ylim(0, 1)

# Save Shanghai K=1 plot
ggsave(filename = "snmf_K1_Shanghai_ancestry_plot.pdf", 
       plot = plot_k1_sh, 
       device = "pdf", 
       width = 12, 
       height = 8, 
       units = "in")

ggsave(filename = "snmf_K1_Shanghai_ancestry_plot.png", 
       plot = plot_k1_sh, 
       device = "png", 
       width = 12, 
       height = 8, 
       units = "in",
       dpi = 300)

cat("Shanghai K=1 plot saved: snmf_K1_Shanghai_ancestry_plot.pdf\n")

#================================================================================
# COMBINED SUMMARY PLOT (OPTIONAL)
#================================================================================

# Create a combined summary showing both regions
combined_data <- bind_rows(
  q_df_ordered_k1_usa %>% mutate(country = "USA"),
  q_df_ordered_k1_sh %>% mutate(country = "Shanghai")
)

# Create combined plot
plot_combined_k1 <- combined_data %>%
  mutate(fill_color = ifelse(country == "USA", usa_color, sh_color)) %>%
  ggplot(aes(x = individual, y = proportion)) +
  geom_col(aes(fill = fill_color), width = 0.9) +
  scale_fill_identity() +
  labs(title = "Spotted Lanternfly Population Structure (K=1)",
       subtitle = "Both USA and Shanghai show panmictic population structure",
       y = "Ancestry Proportion") +
  theme_minimal() +
  facet_wrap(~ paste(country, region), nrow = 2, strip.position = "bottom", scales = "free_x") +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "transparent", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
    legend.position = "none"
  ) +
  ylim(0, 1)

# Save combined plot
ggsave(filename = "snmf_K1_Combined_ancestry_plot.pdf", 
       plot = plot_combined_k1, 
       device = "pdf", 
       width = 20, 
       height = 10, 
       units = "in")

#================================================================================
# CROSS-ENTROPY SUMMARY
#================================================================================

# Create summary of cross-entropy results
summary_data <- data.frame(
  Region = c("USA", "Shanghai"),
  Optimal_K = c(1, 1),
  Cross_Entropy_K1 = c(min(ce_k1_usa), min(ce_k1_sh)),
  N_Samples = c(nrow(popbed_usa), nrow(popbed_sh)),
  N_Locations = c(length(unique(popbed_usa$Location)), length(unique(popbed_sh$Location)))
)

write.csv(summary_data, "SNMF_K1_Summary.csv", row.names = FALSE)

cat("\n=== SNMF K=1 ANALYSIS SUMMARY ===\n")
print(summary_data)



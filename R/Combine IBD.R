# 1. Isolation by geographic distance
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")

library(hierfstat)
library(geodist)
library(ecodist)
library(adespatial)
library(vegan)
library(ggplot2)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")

# Enhanced visualization
IBD_plot <- ggplot(combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance)) +
  geom_point(alpha = 0.8, color = "#3690c0") +
  geom_smooth(method = "lm", color = "black") +   
  theme_minimal() +
  labs(x = "Log Geographic Distance (km)", y = "Genetic Distance (1-IBS)")

SH_IBD_plot <- ggplot(SH_combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance)) +
  geom_point(alpha = 0.8, color = "#3690c0") +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(x = "Log Geographic Distance (km)", y = "Genetic Distance (1-IBS)")


p <- ggplot(combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance)) +
  geom_point(alpha = 0.8, color = "#3690c0") +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  labs(x = "Log Geographic Distance (km)", y = "Pairwise Genetic Distance (1-IBS)")+
  geom_point(data=SH_combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance),
             alpha = 0.8, color = "#c861a5") +
  geom_smooth(data=SH_combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance),
              method = "lm", color = "#e41a1c") +
  theme_classic(base_size = 14) +  # Increase overall base font size
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
  )

# Save as PDF
ggsave("Combined_SLF_IBD_Plot.pdf", plot = p, width = 7, height = 6, units = "in", dpi = 300)

#Part 5: Landscape Genetics
## 5.1 IBD (Isolation By Distance)
# Genetic distances were calculated as 1 minus the identity-by-state (IBS) values using PLINK v1.09 (Purcell et al., 2007), providing a measure of genetic divergence, with values closer to 1 indicating greater dissimilarity. This was done on hpc.

#Isolation by geographic distance


# 1. Isolation by geographic distance
 setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")

library(hierfstat)
library(geodist)
library(ecodist)
library(adespatial)
library(vegan)
library(ggplot2)

 setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")
 
# Read in geo data
xySLF <- read.delim("USA_SLF_MS_GIS.txt")
names(xySLF) <- c("ID", "latitude", "longitude")

# Calculate the geographic distances using the Vincenty method and convert to kilometers
geographic_dist_matrix_km <- geodist(xySLF, measure = "Vincenty") / 1000

# Convert the geographic distance matrix to a dist object
geographic_dist_km <- as.dist(geographic_dist_matrix_km)

# Assuming geo_dist_matrix is a matrix with row and column names that match the IDs in xySLF$ID
geo_dist_matrix <- as.matrix(geographic_dist_matrix_km)
rownames(geographic_dist_matrix_km) <- xySLF$ID
colnames(geographic_dist_matrix_km) <- xySLF$ID

max_value <- max(geo_dist_matrix)
print(paste("Maximum geo distance:", max_value))

# Calculate the mean value of genetic distances
mean_value <- mean(geo_dist_matrix)
print(paste("Mean geo distance:", mean_value))
min(geo_dist_matrix)
max(geo_dist_matrix)

# Change to long format
geo_dist_long <- reshape2::melt(geographic_dist_matrix_km)
names(geo_dist_long) <- c("ID1", "ID2", "GeoDistance")

# Load the genetic distance calculated by plink
# Convert the genetic distance matrix to a dist object
genetic_dist <- read.delim('~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance/USA_1ibs.mdist', header = F)
genetic_dist <- as.matrix(genetic_dist)

max_value <- max(genetic_dist)
print(paste("Maximum genetic distance:", max_value))

# Calculate the mean value of genetic distances
mean_value <- mean(genetic_dist)
print(paste("Mean genetic distance:", mean_value))

# Set the column and row names of the matrix
colnames(genetic_dist) <- xySLF$ID
rownames(genetic_dist) <- xySLF$ID

# Change it to long format
ibs_long <- reshape2::melt(as.matrix(genetic_dist))
names(ibs_long) <- c("ID1", "ID2", "IBSDistance")

# Make sure they are character
geo_dist_long$ID1 <- as.character(geo_dist_long$ID1)
geo_dist_long$ID2 <- as.character(geo_dist_long$ID2)
ibs_long$ID1 <- as.character(ibs_long$ID1)
ibs_long$ID2 <- as.character(ibs_long$ID2)

# Merge the two datasets
combined_data_ibs <- merge(geo_dist_long, ibs_long, by = c("ID1", "ID2"))
combined_data_ibs <- combined_data_ibs[combined_data_ibs$GeoDistance > 0, ] # removing zero distances if any

# Log-transform the distances
combined_data_ibs$GeoDistance_log <- log(combined_data_ibs$GeoDistance + 1)
combined_data_ibs$IBSDistance_log <- log(combined_data_ibs$IBSDistance + 1)

# Mantel test with appropriate transformation
mantel_result_log <- mantel(genetic_dist, as.dist(geo_dist_matrix_log), 
                           method = "pearson", permutations = 1000)
#Mantel statistic r: 0.01081, # Significance: 0.41876 

# Enhanced visualization
IBD_plot <- ggplot(combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance)) +
  geom_point(alpha = 0.8, color = "#3690c0") +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(x = "Log Geographic Distance (km)", y = "Genetic Distance (1-IBS)")

# Saving the plot. Figure 5A.
ggsave("USA_Geo_IBS_log.pdf", plot = IBD_plot, width = 10, height = 8)
ggsave("USA_Geo_IBS_log.png", plot = IBD_plot, dpi = 300)

# Calculate the standard deviation of the geographic distances
sd_geo_distances <- sd(geo_dist_matrix)
sd_geo_distances

# Calculate the mean of the geographic distances
mean_geo_distances <- mean(geo_dist_matrix)
mean_geo_distances

# Function to calculate confidence interval
confidence_interval <- function(data, confidence = 0.95) {
  n <- length(data)
  se <- sd(data) / sqrt(n)
  alpha <- 1 - confidence
  error_margin <- qt(1 - alpha / 2, df = n - 1) * se
  lower_bound <- mean(data) - error_margin
  upper_bound <- mean(data) + error_margin
  return(c(lower_bound, upper_bound))
}

# Flatten the geographic distance matrix to a vector for CI calculation
geo_distance_vector <- as.vector(geo_dist_matrix)

# Calculate the 95% confidence interval for the geographic distances
ci_geo_distances <- confidence_interval(geo_distance_vector)

# Calculate the standard deviation of the genetic distances
sd_gen_distances <- sd(genetic_dist)

# Calculate the mean of the genetic distances
mean_gen_distances <- mean(genetic_dist)

# Flatten the genetic distance matrix to a vector for CI calculation
gen_distance_vector <- as.vector(genetic_dist)

# Calculate the 95% confidence interval for the genetic distances
ci_gen_distances <- confidence_interval(gen_distance_vector)

# Output the results for geographic distances
list(
  max_geo_distance = max(geo_dist_matrix),
  mean_geo_distance = mean_geo_distances,
  min_geo_distance = min(geo_dist_matrix),
  sd_geo_distance = sd_geo_distances,
  ci_geo_lower_bound = ci_geo_distances[1],
  ci_geo_upper_bound = ci_geo_distances[2]
)

# Output the results for genetic distances
list(
  max_gen_distance = max(genetic_dist),
  mean_gen_distance = mean_gen_distances,
  min_gen_distance = min(genetic_dist),
  sd_gen_distance = sd_gen_distances,
  ci_gen_lower_bound = ci_gen_distances[1],
  ci_gen_upper_bound = ci_gen_distances[2]
)



## SH IBD
```{r}
# 1. Isolation by geographic distance
 setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")

library(hierfstat)
library(geodist)
library(ecodist)
library(adespatial)
library(vegan)
library(ggplot2)

 setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance")
 
# Read in geo data
SH_xySLF <- read.delim("SH_SLF_MS_GIS.txt")
names(SH_xySLF) <- c("ID", "latitude", "longitude")

# Calculate the geographic distances using the Vincenty method and convert to kilometers
SH_geographic_dist_matrix_km <- geodist(SH_xySLF, measure = "Vincenty") / 1000

# Convert the geographic distance matrix to a dist object
SH_geographic_dist_km <- as.dist(SH_geographic_dist_matrix_km)

# Assuming geo_dist_matrix is a matrix with row and column names that match the IDs in xySLF$ID
SH_geo_dist_matrix <- as.matrix(SH_geographic_dist_matrix_km)
rownames(SH_geographic_dist_matrix_km) <- SH_xySLF$ID
colnames(SH_geographic_dist_matrix_km) <- SH_xySLF$ID

max_value <- max(geo_dist_matrix)
print(paste("Maximum geo distance:", max_value))

# Calculate the mean value of genetic distances
mean_value <- mean(geo_dist_matrix)
print(paste("Mean geo distance:", mean_value))
min(geo_dist_matrix)
max(geo_dist_matrix)

# Change to long format
SH_geo_dist_long <- reshape2::melt(SH_geographic_dist_matrix_km)
names(SH_geo_dist_long) <- c("ID1", "ID2", "GeoDistance")

# Load the genetic distance calculated by plink
# Convert the genetic distance matrix to a dist object
SH_genetic_dist <- read.delim('~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Genetic_Distance/Shanghai_1ibs.mdist', header = F)
SH_genetic_dist <- as.matrix(SH_genetic_dist)

SH_max_value <- max(SH_genetic_dist)
print(paste("Maximum genetic distance:", SH_max_value))

# Calculate the mean value of genetic distances
mean_value <- mean(genetic_dist)
print(paste("Mean genetic distance:", mean_value))

# Set the column and row names of the matrix
colnames(SH_genetic_dist) <- SH_xySLF$ID
rownames(SH_genetic_dist) <- SH_xySLF$ID

# Change it to long format
SH_ibs_long <- reshape2::melt(as.matrix(SH_genetic_dist))
names(SH_ibs_long) <- c("ID1", "ID2", "IBSDistance")

# Make sure they are character
SH_geo_dist_long$ID1 <- as.character(SH_geo_dist_long$ID1)
SH_geo_dist_long$ID2 <- as.character(SH_geo_dist_long$ID2)
SH_ibs_long$ID1 <- as.character(SH_ibs_long$ID1)
SH_ibs_long$ID2 <- as.character(SH_ibs_long$ID2)

# Merge the two datasets
SH_combined_data_ibs <- merge(SH_geo_dist_long, SH_ibs_long, by = c("ID1", "ID2"))
SH_combined_data_ibs <- SH_combined_data_ibs[SH_combined_data_ibs$GeoDistance > 0, ] # removing zero distances if any

# Log-transform the distances
SH_combined_data_ibs$GeoDistance_log <- log(SH_combined_data_ibs$GeoDistance + 1)
SH_combined_data_ibs$IBSDistance_log <- log(SH_combined_data_ibs$IBSDistance + 1)

# Mantel test with appropriate transformation
SH_mantel_result_sh <- mantel(SH_genetic_dist, SH_geo_dist_matrix, 
                           method = "pearson", permutations = 1000)
#Mantel statistic r: 0.91, # Significance: 0.000999  

# Enhanced visualization
SH_IBD_plot <- ggplot(SH_combined_data_ibs, aes(x = GeoDistance_log, y = IBSDistance)) +
  geom_point(alpha = 0.8, color = "#3690c0") +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(x = "Log Geographic Distance (km)", y = "Genetic Distance (1-IBS)")

# Saving the plot. Figure 5A.
ggsave("SH_Geo_IBS_log.pdf", plot = IBD_plot, width = 8, height = 6)
ggsave("SH_Geo_IBS_log.png", plot = IBD_plot, dpi = 300)

# Calculate the standard deviation of the geographic distances
sd_geo_distances <- sd(geo_dist_matrix)
sd_geo_distances

# Calculate the mean of the geographic distances
mean_geo_distances <- mean(geo_dist_matrix)
mean_geo_distances

# Function to calculate confidence interval
confidence_interval <- function(data, confidence = 0.95) {
  n <- length(data)
  se <- sd(data) / sqrt(n)
  alpha <- 1 - confidence
  error_margin <- qt(1 - alpha / 2, df = n - 1) * se
  lower_bound <- mean(data) - error_margin
  upper_bound <- mean(data) + error_margin
  return(c(lower_bound, upper_bound))
}

# Flatten the geographic distance matrix to a vector for CI calculation
geo_distance_vector <- as.vector(geo_dist_matrix)

# Calculate the 95% confidence interval for the geographic distances
ci_geo_distances <- confidence_interval(geo_distance_vector)

# Calculate the standard deviation of the genetic distances
sd_gen_distances <- sd(genetic_dist)

# Calculate the mean of the genetic distances
mean_gen_distances <- mean(genetic_dist)

# Flatten the genetic distance matrix to a vector for CI calculation
gen_distance_vector <- as.vector(genetic_dist)

# Calculate the 95% confidence interval for the genetic distances
ci_gen_distances <- confidence_interval(gen_distance_vector)

# Output the results for geographic distances
list(
  max_geo_distance = max(geo_dist_matrix),
  mean_geo_distance = mean_geo_distances,
  min_geo_distance = min(geo_dist_matrix),
  sd_geo_distance = sd_geo_distances,
  ci_geo_lower_bound = ci_geo_distances[1],
  ci_geo_upper_bound = ci_geo_distances[2]
)

# Output the results for genetic distances
list(
  max_gen_distance = max(genetic_dist),
  mean_gen_distance = mean_gen_distances,
  min_gen_distance = min(genetic_dist),
  sd_gen_distance = sd_gen_distances,
  ci_gen_lower_bound = ci_gen_distances[1],
  ci_gen_upper_bound = ci_gen_distances[2]
)

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
ggsave("Combined_SLF_IBD_Plot.pdf", plot = p, width = 8, height = 6, units = "in", dpi = 300)
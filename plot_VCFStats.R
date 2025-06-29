# **1. Set Working Directories**
setwd("/scratch/work/winchell/Lanternfly/WGS/15_NoDownsample_VCF/VCF_Stats")


# **2. Load Required Libraries**
library(data.table)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyverse)

# **3. Read and Explore VCF Data**
vcf_data <- fread('/scratch/work/winchell/Lanternfly/Aria/11_joint_genotyping/SLF.SNPs.vcf.table', header = TRUE, sep = "\t")
colnames(vcf_data)

# **4. Plot Feature Distributions**
output_pdf <- "VCF_Feature_Distributions.pdf"
features <- c("QUAL", "MQ", "SOR", "QD", "FS", "MQRankSum", "ReadPosRankSum")

pdf(output_pdf, width = 8, height = 6)
for (feature in features) {
  p <- ggplot(vcf_data, aes_string(x = feature)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    ggtitle(paste("Density Distribution of", feature)) +
    xlab(feature) +
    ylab("Density") +
    theme_minimal()
  print(p)
}
dev.off()
cat("Plots have been saved to", output_pdf, "\n")

# **5. Plot Specific Features with Filtering**
filtered_data <- vcf_data[vcf_data$QUAL < 1000, ]
qual_plot <- ggplot(filtered_data, aes(x = QUAL)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  ggtitle("Density Distribution of QUAL (QUAL < 1000)") +
  xlab("QUAL") + ylab("Density") + theme_minimal()
ggsave("QUAL_Density_Filtered.pdf", plot = qual_plot, width = 12, height = 10)

MQ_plot <- ggplot(vcf_data, aes(x = MQ)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  ggtitle("Density Distribution of MQ") +
  xlab("MQ") + ylab("Density") + theme_minimal() + xlim(0, 100)
ggsave("MQ_Density_Filtered.pdf", plot = MQ_plot, width = 12, height = 10)

# **6. Analyze VCF Statistics**

# **6.1 Depth Analysis**
var_depth <- read_delim("./SLF_Stats_site_mean_depth.ldepth.mean", 
                        delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

depth_plot <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  xlim(0, 20) + theme_light()
ggsave("mean_depth_density_plot.pdf", plot = depth_plot, width = 8, height = 6)
summary(var_depth$mean_depth)

# depth per individual
# read the data
dp <- read.table("./SLF_Stats_indv_depth.idepth", header =T)

# Plot the distribution of MEAN_COVERAGE with specified x-axis range and binwidth
MEAN_dp <- ggplot(dp, aes(x = INDV, y = MEAN_DEPTH)) +
  geom_bar(stat = "identity",fill = "lightblue") +
  theme_minimal() +
  labs(x = "Individual", y = "mean_depth") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x labels for better readability

ggsave(filename = "MEAN_dp_IND.pdf", plot = MEAN_dp, width = 10, height = 6)
summary(dp$MEAN_DEPTH)



# **6.2 Missing Data per site**
missingness <- read.table("SLF_Stats_missing_site.lmiss", header = TRUE)

missingness_plot <- ggplot(missingness, aes(F_MISS)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "missingness_density.pdf", plot = missingness_plot, width = 10, height = 6)

summary(missingness$F_MISS)

# missing data per individual
missingness_IND <- read.table("SLF_Stats_missing_indv.imiss", header = TRUE)
missingness_IND_plot <- ggplot(missingness_IND, aes(x = INDV, y = F_MISS)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Individual", y = "Frequency of Missing Data") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "missingness_IND.pdf", plot = missingness_IND_plot, width = 10, height = 6)
summary(missingness_IND$F_MISS)


# **6.3 Minor allele frequency**
var_freq <- read_delim("./SLF_Stats_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

maf_plot <- ggplot(var_freq, aes(maf)) + 
  geom_histogram(binwidth = 0.01, fill = "dodgerblue1", colour = "black", alpha = 0.7) +
  theme_light() + labs(title = "Distribution of Minor Allele Frequencies",
                       x = "Minor Allele Frequency (MAF)", y = "Count")

ggsave("SLF_combined_MAF.pdf", plot = maf_plot)

summary(var_freq$maf)



# **6.3 Heterozygosity and Inbreeding Coefficient**
hetstats <- read.table("Combined_het.het", header = TRUE)
popdata <- read.csv("2022 SLF Info new.csv", header = TRUE)
hetstats$O.HET <- hetstats$N_SITES - hetstats$O.HOM.
hetstats$E.HET <- hetstats$N_SITES - hetstats$E.HOM.
hetstats$Ho <- hetstats$O.HET / hetstats$N_SITES
hetstats$He <- hetstats$E.HET / hetstats$N_SITES
hetstats <- merge(hetstats, popdata, by.x = "INDV", by.y = "SLF_ID")

# **6.4 Minor Allele Frequency (MAF) Analysis**
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Function to process and plot MAF
process_maf <- function(file_name, dataset_name, output_pdf) {
  # Read the frequency data
  var_freq <- read_delim(file_name, 
                         delim = "\t", 
                         col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), 
                         skip = 1)
  
  # Calculate Minor Allele Frequency (MAF)
  var_freq$maf <- apply(var_freq %>% select(a1, a2), 1, function(z) min(z))
  
  # Plot the MAF distribution
  maf_plot <- ggplot(var_freq, aes(maf)) + 
    geom_histogram(binwidth = 0.01, fill = "dodgerblue1", colour = "black", alpha = 0.7) +
    theme_light() + 
    labs(
      title = paste("MAF Distribution:", dataset_name),
      x = "Minor Allele Frequency (MAF)", 
      y = "Count"
    )
  
  # Save the plot as a PDF
  ggsave(output_pdf, plot = maf_plot, width = 8, height = 6)
  
  # Print MAF summary statistics
  print(paste("Summary for", dataset_name))
  print(summary(var_freq$maf))
}

# Process and plot MAF for Shanghai and USA datasets
process_maf("Shanghai_samples_Stats_freq.frq", "Shanghai", "Shanghai_MAF_plot.pdf")
process_maf("USA_samples_Stats_freq.frq", "USA", "USA_MAF_plot.pdf")



var_freq <- read_delim("./USA_samples_Stats_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

# Ensure MAF values are > 0
var_freq_filtered <- var_freq %>% filter(maf > 0)

# Define updated MAF ranges including 10-20 and >20 alleles
maf_ranges <- c(0, 1 / (2 * sample_size), 2 / (2 * sample_size), 4 / (2 * sample_size), 
                10 / (2 * sample_size), 20 / (2 * sample_size), 1.0)

# Update the var_freq_filtered dataframe with new MAF groups
var_freq_filtered <- var_freq_filtered %>%
  mutate(
    maf_group = cut(
      maf, 
      breaks = maf_ranges, 
      include.lowest = TRUE, 
      labels = c("Singleton", "1-2 Alleles", "2-4 Alleles", "4-10 Alleles", "10-20 Alleles", ">20 Alleles"),
      right = TRUE
    )
  )

# Count the number of sites in each group
maf_counts <- var_freq_filtered %>%
  group_by(maf_group) %>%
  summarise(site_count = n(), .groups = "drop")

# Total number of SNPs
total_snps <- sum(maf_counts$site_count)

# Add a percentage column to maf_counts
maf_counts <- maf_counts %>%
  mutate(percentage = (site_count / total_snps) * 100)

# Print the updated maf_counts with percentages
print(maf_counts)

# Write the results to a CSV file
write_csv(maf_counts, "USA_maf_counts.csv")

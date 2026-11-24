# **1. Set Working Directories**
setwd('/scratch/work/winchell/Lanternfly/Aria/12_Final_VCF')
setwd('/scratch/work/winchell/Lanternfly/Aria/11_joint_genotyping')

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
var_depth <- read_delim("./Shanghai_site_mean_depth.ldepth.mean", 
                        delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

depth_plot <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  xlim(0, 20) + theme_light()
ggsave("mean_depth_density_plot.pdf", plot = depth_plot, width = 8, height = 6)
summary(var_depth$mean_depth)

# **6.2 Missing Data**
missingness <- read.table("Combined_missing.imiss", header = TRUE)
missingness_IND <- ggplot(missingness, aes(x = INDV, y = F_MISS)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Individual", y = "Frequency of Missing Data") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "missingness_IND.pdf", plot = missingness_IND, width = 10, height = 6)

# **6.3 Heterozygosity and Inbreeding Coefficient**
hetstats <- read.table("Combined_het.het", header = TRUE)
popdata <- read.csv("2022 SLF Info new.csv", header = TRUE)
hetstats$O.HET <- hetstats$N_SITES - hetstats$O.HOM.
hetstats$E.HET <- hetstats$N_SITES - hetstats$E.HOM.
hetstats$Ho <- hetstats$O.HET / hetstats$N_SITES
hetstats$He <- hetstats$E.HET / hetstats$N_SITES
hetstats <- merge(hetstats, popdata, by.x = "INDV", by.y = "SLF_ID")

# **6.4 Minor Allele Frequency (MAF) Analysis**
process_maf <- function(file_name, output_pdf) {
  var_freq <- read_delim(file_name, 
                         delim = "\t", 
                         col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), 
                         skip = 1)
  var_freq$maf <- apply(var_freq %>% select(a1, a2), 1, function(z) min(z))
  maf_plot <- ggplot(var_freq, aes(maf)) + 
    geom_histogram(binwidth = 0.01, fill = "dodgerblue1", colour = "black", alpha = 0.7) +
    theme_light() + labs(title = "Distribution of Minor Allele Frequencies",
                         x = "Minor Allele Frequency (MAF)", y = "Count")
  ggsave(output_pdf, plot = maf_plot, width = 8, height = 6)
  print(summary(var_freq$maf))
}

process_maf("Shanghai_freq2.frq", "Shanghai_MAF_plot.pdf")
process_maf("USA_freq2.frq", "USA_MAF_plot.pdf")




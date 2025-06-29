####### Look for shared 185 SNPs Stats
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/Overlap")

# Load heterozygosity file: SH_rural
het_rural <- read.table("SH_rural_185_het.het", 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        row.names = NULL, 
                        check.names = FALSE)

# Load SH_urban
het_urban <- read.table("SH_urban_185_het.het", 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        row.names = NULL, 
                        check.names = FALSE)

# Load USA
het_usa <- read.table("USA_185_het.het", 
                      header = TRUE, 
                      stringsAsFactors = FALSE, 
                      row.names = NULL, 
                      check.names = FALSE)

# -----------------------------
# Define function to calculate heterozygosity stats
# -----------------------------
calculate_het_stats <- function(df, pop_name) {
  df$Population <- pop_name
  df$O.HET <- df$N_SITES - df$`O(HOM)`  # Observed heterozygous sites
  df$E.HET <- df$N_SITES - df$`E(HOM)`  # Expected heterozygous sites
  df$Ho <- df$O.HET / df$N_SITES        # Proportion observed
  df$He <- df$E.HET / df$N_SITES        # Proportion expected
  return(df)
}

# -----------------------------
# Apply to each dataset
# -----------------------------
het_rural <- calculate_het_stats(het_rural, "SH_rural")
het_urban <- calculate_het_stats(het_urban, "SH_urban")
het_usa   <- calculate_het_stats(het_usa,   "USA")

# -----------------------------
# Combine all into one data frame
# -----------------------------
het_all <- bind_rows(het_rural, het_urban, het_usa)
write.csv()
# -----------------------------
# Summary of heterozygosity by population
# -----------------------------
het_summary <- het_all %>%
  group_by(Population) %>%
  summarise(
    mean_Ho = mean(Ho, na.rm = TRUE),
    mean_He = mean(He, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    sd_He = sd(He, na.rm = TRUE),
    n = n()
  )

print(het_summary)


######### Try to plot
# Convert the combined dataframe for plotting
het_plot_data <- het_all %>%
  select(INDV, Population, Ho) %>%
  rename(Group = Population, Heterozygosity = Ho)

# Re-level Group to get desired order on x-axis
het_plot_data$Group <- factor(het_plot_data$Group, levels = c("SH_rural", "SH_urban", "USA"))

# Load mean_se function if not already available
if(!exists("mean_se")) {
  mean_se <- function(x) {
    mean_x <- mean(x, na.rm = TRUE)
    se_x <- sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
    return(data.frame(y = mean_x, ymin = mean_x - se_x, ymax = mean_x + se_x))
  }
}

# Create violin plot with matching point colors
het_violin_plot <- ggplot(het_plot_data, aes(x = Group, y = Heterozygosity, fill = Group)) +
  geom_violin(alpha = 0.7) +
  geom_point(aes(color = Group), alpha = 0.8, size = 3, position = position_jitter(width = 0.15)) +  # Individual points with jitter
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 95, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("SH_rural" = "#009E73", "SH_urban" = "#c861a5", "USA" = "#3690c0"), guide = FALSE) +
  scale_color_manual(values = c("SH_rural" = "#009E73", "SH_urban" = "#c861a5", "USA" = "#3690c0"), guide = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Observed Heterozygosity (Ho)\nat Selected SNPs") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

# Print the plot
print(het_violin_plot)

# Optional: Save the plot
ggsave("185_heterozygosity_at_selected_SNPs.pdf", het_violin_plot, width = 7, height = 5)

# Random subsampling approach
set.seed(123) # For reproducibility
n_subsamples <- 100 # Number of iterations
sample_size <- 9 # Match the smallest group

# Create storage for results
subsample_results <- data.frame(
  iteration = integer(),
  mean_ho_usa = numeric()
)

# Run subsampling
for(i in 1:n_subsamples) {
  # Randomly sample from USA population
  usa_sample <- het_usa %>% 
    sample_n(sample_size) 
  
  # Calculate mean Ho for this sample
  mean_ho <- mean(usa_sample$Ho, na.rm = TRUE)
  
  # Store result
  subsample_results <- subsample_results %>%
    add_row(iteration = i, mean_ho_usa = mean_ho)
}

# Calculate confidence intervals
ci_usa <- quantile(subsample_results$mean_ho_usa, probs = c(0.025, 0.5, 0.975))

library(boot)

# Bootstrap function
boot_mean <- function(data, indices) {
  mean(data[indices])
}

# Apply to USA Ho values
usa_ho_values <- het_usa$Ho
usa_boot <- boot(usa_ho_values, boot_mean, R = 1000)
usa_ci <- boot.ci(usa_boot, type = "perc")

# Now compare with actual means from SH populations
mean_ho_rural <- mean(het_rural$Ho, na.rm = TRUE)
mean_ho_urban <- mean(het_urban$Ho, na.rm = TRUE)

# Kruskal-Wallis test since data may not be normally distributed
kruskal_test <- kruskal.test(Ho ~ Population, data = het_all)

# Post-hoc pairwise comparisons
if(kruskal_test$p.value < 0.05) {
  pairwise_results <- pairwise.wilcox.test(
    het_all$Ho, 
    het_all$Population,
    p.adjust.method = "bonferroni"
  )
}

# Alternatively, using subsampled data for USA
# Create multiple datasets with equal sample sizes
balanced_comparisons <- list()

for(i in 1:30) { # Perform 30 balanced comparisons
  usa_sample <- het_usa %>% sample_n(9)
  balanced_data <- bind_rows(
    het_rural,
    het_urban,
    usa_sample %>% mutate(Population = "USA_subsample")
  )
  
  # Store p-value from Kruskal-Wallis test
  kw_result <- kruskal.test(Ho ~ Population, data = balanced_data)
  balanced_comparisons[[i]] <- kw_result$p.value
}

# Check how many tests were significant
sum(balanced_comparisons < 0.05) / length(balanced_comparisons)

library(dplyr)
library(readr)
library(data.table)

# Read allele frequencies
SH_urban <- fread("SH_urban_185_freq.frq", header = F)
SH_rural <- fread("SH_rural_185_freq.frq", header = F)
USA <- fread("USA_185_freq.frq", header = F)

# Rename columns for clarity
colnames(SH_urban) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ")
colnames(SH_rural) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ")
colnames(USA) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ")

# Extract numeric frequencies from character strings
# Function to extract frequency values from strings like "T:0.636364"
extract_freq <- function(freq_str) {
  # Split by colon and take the second part
  as.numeric(strsplit(freq_str, ":")[[1]][2])
}

# Apply to each dataset - extract ALT frequencies
SH_urban$ALT_FREQ_VALUE <- sapply(SH_urban$ALT_FREQ, extract_freq)
SH_rural$ALT_FREQ_VALUE <- sapply(SH_rural$ALT_FREQ, extract_freq)
USA$ALT_FREQ_VALUE <- sapply(USA$ALT_FREQ, extract_freq)

# Create SNP_ID for joining
SH_urban$SNP_ID <- paste(SH_urban$CHROM, SH_urban$POS, sep = "_")
SH_rural$SNP_ID <- paste(SH_rural$CHROM, SH_rural$POS, sep = "_")
USA$SNP_ID <- paste(USA$CHROM, USA$POS, sep = "_")

# Create a merge-ready dataframe with just SNP_ID and frequency
SH_urban_freqs <- SH_urban %>% select(SNP_ID, ALT_FREQ_VALUE)
SH_rural_freqs <- SH_rural %>% select(SNP_ID, ALT_FREQ_VALUE)
USA_freqs <- USA %>% select(SNP_ID, ALT_FREQ_VALUE)

# Rename columns for clarity
colnames(SH_urban_freqs)[2] <- "SH_urban_freq"
colnames(SH_rural_freqs)[2] <- "SH_rural_freq"
colnames(USA_freqs)[2] <- "USA_freq"

# Merge datasets
merged_freqs <- SH_urban_freqs %>%
  inner_join(SH_rural_freqs, by = "SNP_ID") %>%
  inner_join(USA_freqs, by = "SNP_ID")

# Check merged data
head(merged_freqs)
cat("Total SNPs after merging:", nrow(merged_freqs), "\n")

annot_combined$SNP_ID <- gsub("_[ACGT]$", "", annot_combined$SNP) 

freqs_merged <- merge(merged_freqs, annot_combined, by = "SNP_ID", all.x = TRUE)
write.csv(freqs_merged , "Final185_freq.csv")

# Calculate distances between USA and each Shanghai population
merged_freqs$USA_urban_dist <- abs(merged_freqs$USA_freq - merged_freqs$SH_urban_freq)
merged_freqs$USA_rural_dist <- abs(merged_freqs$USA_freq - merged_freqs$SH_rural_freq)
merged_freqs$urban_rural_diff <- abs(merged_freqs$SH_urban_freq - merged_freqs$SH_rural_freq)

# Calculate which population is closer to USA for each SNP
merged_freqs$closer_to_USA <- ifelse(
  merged_freqs$USA_urban_dist < merged_freqs$USA_rural_dist, 
  "SH_urban", "SH_rural"
)

# Overall summary
cat("Mean USA-Urban Distance:", mean(merged_freqs$USA_urban_dist), "\n")
cat("Mean USA-Rural Distance:", mean(merged_freqs$USA_rural_dist), "\n")
cat("Percentage of SNPs where USA is closer to Urban:", 
    mean(merged_freqs$closer_to_USA == "SH_urban") * 100, "%\n")

# Statistical test - paired t-test
t_test_result <- t.test(
  merged_freqs$USA_urban_dist, 
  merged_freqs$USA_rural_dist, 
  paired = TRUE, 
  alternative = "less"  # Testing if USA-urban distance is less than USA-rural
)
print(t_test_result)


library(dplyr)

# Check total SNPs before deduplication
cat("Original number of rows:", nrow(pi_wide_merged), "\n")

# Deduplicate: keep only one row per SNP_ID (or POS)
pi_unique <- pi_wide_merged %>%
  group_by(SNP_ID) %>%
  slice(1) %>%  # Keep only the first row for each SNP_ID
  ungroup()

# Verify number of rows after deduplication
cat("Unique SNPs after deduplication:", nrow(pi_unique), "\n")

# Optional: keep only relevant columns if desired
pi_unique_clean <- pi_unique %>%
  select(SNP_ID, POS, CHROM.x, SH_rural, USA, SH_urban,
         diff_urban_rural, diff_usa_SHrual, diff_usa_SHurban,
         Gene_Name, Annotation, REF, ALT,
         predictor, correlation)

# Preview
head(pi_unique_clean)

# Step 1: Calculate distances between π values
pi_unique_clean <- pi_unique_clean %>%
  mutate(
    pi_dist_USA_urban = abs(USA - SH_urban),
    pi_dist_USA_rural = abs(USA - SH_rural),
    pi_dist_urban_rural = abs(SH_urban - SH_rural),
    closer_pi_to_USA = ifelse(pi_dist_USA_urban < pi_dist_USA_rural, "SH_urban", "SH_rural")
  )

# Step 2: Summary statistics
cat("Mean π Distance (USA vs SH_urban):", mean(pi_unique_clean$pi_dist_USA_urban, na.rm = TRUE), "\n")
cat("Mean π Distance (USA vs SH_rural):", mean(pi_unique_clean$pi_dist_USA_rural, na.rm = TRUE), "\n")
cat("Percentage of SNPs where USA π is closer to Urban:",
    mean(pi_unique_clean$closer_pi_to_USA == "SH_urban", na.rm = TRUE) * 100, "%\n")

# Step 3: Paired t-test — is USA π closer to urban than rural?
pi_ttest <- t.test(
  pi_unique_clean$pi_dist_USA_urban,
  pi_unique_clean$pi_dist_USA_rural,
  paired = TRUE,
  alternative = "less"  # Test if urban is closer to USA than rural
)
print(pi_ttest)

write.csv(pi_unique_clean, "Final185_pi_unique_clean.csv")


##### Look for the overlap of the genes of PCAdpt of USA and CN pacapt
# Function to extract unique functional genes (not LDel_XXX and not intergenic regions)
extract_functional_genes <- function(df) {
  # Create a vector to store functional gene names
  functional_genes <- character()
  
  for (i in 1:nrow(df)) {
    # Get the gene name
    gene_name <- df$Gene_Name[i]
    
    # Skip if the gene is in an intergenic region (contains a hyphen)
    if (grepl("-", gene_name)) {
      next
    }
    
    # Skip if the gene name starts with LDel_
    if (grepl("^LDel_", gene_name)) {
      next
    }
    
    # Add the gene to our list
    functional_genes <- c(functional_genes, gene_name)
  }
  
  # Return unique genes
  return(unique(functional_genes))
}

# Extract unique functional genes from each dataset
cn_genes <- extract_functional_genes(fdr_merged)
usa_genes <- extract_functional_genes(fdr_USA_Only)

# Find the overlap
overlapping_genes <- intersect(cn_genes, usa_genes)

shared_genes <- extract_functional_genes(shared_SNPs)

shared_genes_df <- merge


# Print results
cat("Number of unique functional genes in China dataset:", length(cn_genes), "\n")
cat("Number of unique functional genes in USA dataset:", length(usa_genes), "\n")
cat("Number of overlapping genes:", length(overlapping_genes), "\n")

# Print the overlapping genes
cat("\nOverlapping genes:\n")
sort(overlapping_genes)

shared_SNPs <- data.frame(SNP=intersect(fdr_merged$SNP, fdr_USA_Only$SNP))
shared_SNPs <- merge(fdr_USA_Only,shared_SNPs, all.x=F, all.y=F, by="SNP" )
sort(unique(shared_SNPs$Gene_Name))

# Save the results to files
write.csv(data.frame(Gene = cn_genes), "cn_functional_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = usa_genes), "usa_functional_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = overlapping_genes), "overlapping_genes.csv", row.names = FALSE)

# Create a Venn diagram (if VennDiagram package is installed)
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  library(VennDiagram)
  library(grid)
  
  # Create a Venn diagram
  venn_plot <- draw.pairwise.venn(
    area1 = length(cn_genes),
    area2 = length(usa_genes),
    cross.area = length(overlapping_genes),
    category = c("China", "USA"),
    fill = c("lightblue", "lightgreen"),
    alpha = 0.5,
    cat.pos = c(0, 0),
    cat.dist = c(0.025, 0.025)
  )
  
  # Save the Venn diagram
  pdf("pcadapt_gene_overlap.pdf", width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
  
  # For those without VennDiagram package, create a simple plot
} else {
  # Create a simple barplot
  counts <- c(length(cn_genes), length(usa_genes), length(overlapping_genes))
  labels <- c("China only", "USA only", "Overlap")
  
  pdf("pcadapt_gene_counts.pdf", width = 8, height = 6)
  barplot(counts, names.arg = labels, col = c("lightblue", "lightgreen", "coral"),
          main = "Functional Genes Under Selection", ylab = "Count")
  dev.off()
}

# Create a data frame to summarize results
summary_df <- data.frame(
  Dataset = c("China", "USA", "Overlap"),
  Gene_Count = c(length(cn_genes), length(usa_genes), length(overlapping_genes))
)

# Save summary
write.csv(summary_df, "pcadapt_gene_summary.csv", row.names = FALSE)


# Find genes unique to each dataset
unique_cn_genes <- setdiff(cn_genes, usa_genes)
unique_usa_genes <- setdiff(usa_genes, cn_genes)

# Print results
cat("Number of genes unique to China dataset:", length(unique_cn_genes), "\n")
cat("Number of genes unique to USA dataset:", length(unique_usa_genes), "\n")

# Save the unique gene lists
write.csv(data.frame(Gene = unique_cn_genes), "unique_cn_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = unique_usa_genes), "unique_usa_genes.csv", row.names = FALSE)


# Additional analysis: Check for functional patterns in gene names
# This is a simple approach to see if certain gene families are overrepresented

# Function to extract gene family prefixes (e.g., "CYP" from "CYP4V2_1")
extract_gene_family <- function(gene_vector) {
  # Extract prefixes like CYP, ABC, etc.
  families <- sapply(gene_vector, function(x) {
    # Extract letters at the beginning of the gene name
    match <- regexpr("^[A-Z]+", x)
    if (match > 0) {
      return(substr(x, match, match + attr(match, "match.length") - 1))
    } else {
      return(NA)
    }
  })
  return(families[!is.na(families)])
}

# Analyze gene families in each dataset
cn_families <- extract_gene_family(unique_cn_genes)
usa_families <- extract_gene_family(unique_usa_genes)

# Count frequencies
cn_family_counts <- table(cn_families)
usa_family_counts <- table(usa_families)

# Print top gene families
cat("\nTop gene families in China dataset:\n")
print(sort(cn_family_counts, decreasing = TRUE)[1:min(5, length(cn_family_counts))])

cat("\nTop gene families in USA dataset:\n")
print(sort(usa_family_counts, decreasing = TRUE)[1:min(5, length(usa_family_counts))])

# Save gene family analysis
write.csv(data.frame(Family = names(cn_family_counts), 
                     Count = as.numeric(cn_family_counts)),
          "cn_gene_families.csv", row.names = FALSE)

write.csv(data.frame(Family = names(usa_family_counts),
                     Count = as.numeric(usa_family_counts)),
          "usa_gene_families.csv", row.names = FALSE)


####### SH pcadpt (760 SNPs) and USA pcadpt (7343 SNPs), only has 19 SNPs overlap
shared_SNPs$SNP

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(grid)
library(gridExtra)

# Filter merged_freqs to only include the shared SNPs
shared_snps_data <- merged_freqs %>%
  filter(SNP_ID %in% shared_SNPs$SNP) %>%
  # Order by USA_freq in descending order
  arrange(desc(USA_freq)) %>%
  select(SNP_ID, SH_rural_freq, SH_urban_freq, USA_freq)

# Create a matrix for the heatmap
shared_matrix <- as.matrix(shared_snps_data[, c("SH_rural_freq", "SH_urban_freq", "USA_freq")])
rownames(shared_matrix) <- shared_snps_data$SNP_ID
colnames(shared_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Convert data to long format for ggplot
shared_long <- shared_snps_data %>%
  pivot_longer(
    cols = c(SH_rural_freq, SH_urban_freq, USA_freq),
    names_to = "Population",
    values_to = "Frequency"
  ) %>%
  mutate(
    Population = factor(Population, 
                        levels = c("SH_rural_freq", "SH_urban_freq", "USA_freq"),
                        labels = c("SH_Rural", "SH_Urban", "USA"))
  )

# Alternative: Custom pheatmap with better text contrast
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Calculate text color (white for dark cells, black for light cells)
  cell_value <- shared_matrix[i, j]
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  
  # Add text with appropriate color
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 8))
}

# Create pheatmap with custom cell rendering
shared_heatmap <- pheatmap(
  shared_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,  # We'll add custom text with custom_cell_fun
  fontsize_row = 8,
  cellwidth = 40,
  cellheight = 15,
  border_color = NA,
  main = "Allele Frequencies of 19 Shared SNPs\nOrdered by USA frequency (highest to lowest)",
  fontsize = 10,
  cell.fun = custom_cell_fun  # Custom function for cell rendering
)

# Save to PDF
pdf("shared_SNPs_pheatmap.pdf", width = 8, height = 10)
shared_heatmap
dev.off()


# Define custom cell function to display values with optimal contrast
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Get the cell value
  cell_value <- shared_matrix[i, j]
  
  # Determine text color based on background darkness
  # Dark backgrounds get white text, light backgrounds get black text
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  
  # Add text with appropriate color and formatting
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Create pheatmap with custom cell rendering
shared_heatmap <- pheatmap(
  shared_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,  # We'll use custom_cell_fun instead
  fontsize_row = 9,
  cellwidth = 50,
  cellheight = 18,
  border_color = "white",  # White borders between cells
  main = "Allele Frequencies of 19 Shared SNPs\nOrdered by USA frequency (highest to lowest)",
  fontsize = 12,
  cell.fun = custom_cell_fun  # Custom function for cell rendering
)

# Save to PDF
pdf("shared_SNPs_pheatmap.pdf", width = 8, height = 12)
shared_heatmap
dev.off()


######## Look the overlap of shanghai 760 SNPs and 2675 RDA SNPs in the USA, especially urban ones
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/AlleleFrequencies")
prda_SLF_annotation <- read_csv("prda_SLF_annotation.csv")

# Merge prda_SLF_annotation and merged_freqs_annoted datasets by SNP ID

# First, standardize the SNP ID format in prda_SLF_annotation
# Convert scientific notation positions to regular integers and create consistent IDs
standardize_prda_SNP_IDs <- function(data) {
  # Extract chromosome and position from SNP column
  data$SNP_ID <- paste0(
    data$CHROM,
    "_",
    as.integer(data$POS)  # Convert scientific notation to integer
  )
  return(data)
}

# Standardize the prda_SLF_annotation dataset
prda_SLF_standardized <- standardize_prda_SNP_IDs(prda_SLF_annotation)

# Check if SNP column in merged_freqs_annoted needs processing
# From the sample data, it appears SNP_ID is already in the right format

# Perform the merge by SNP ID
merged_data <- merge(
  prda_SLF_standardized,
  merged_freqs_annoted,
  by = "SNP_ID",
  all = FALSE  # Keep only matching rows
)

# Alternative approach if above doesn't work - create a common key first
if (nrow(merged_data) == 0) {
  # Extract parts from SNP in prda_SLF_annotation
  prda_SLF_standardized <- prda_SLF_annotation %>%
    mutate(
      chr_num = gsub("chr_", "", CHROM),
      pos_int = as.integer(POS),
      SNP_ID = paste0("chr_", chr_num, "_", pos_int)
    )
  
  # Try to match the format in merged_freqs_annoted
  merged_freqs_processed <- merged_freqs_annoted %>%
    mutate(SNP_ID_simple = sub("_[ACGT]$", "", SNP))  # Remove allele suffix if present
  
  # Merge using the processed IDs
  merged_data <- merge(
    prda_SLF_standardized,
    merged_freqs_processed,
    by.x = "SNP_ID",
    by.y = "SNP_ID_simple",
    all = FALSE
  )
}

# Check results
print(paste("Number of rows in merged dataset:", nrow(merged_data)))
print(paste("Number of columns in merged dataset:", ncol(merged_data)))

# Display the first few rows to verify
head(merged_data)

# Save the merged dataset
write.csv(merged_data, "SLF_merged_complete_dataset.csv", row.names = FALSE)

# Summary of data in merged dataset
cat("\nSummary of merged data:\n")

# Count of SNPs by chromosome
cat("\nSNPs per chromosome:\n")
table(merged_data$CHROM.x)  # Using .x suffix from the merge

# Count by annotation type
cat("\nDistribution by annotation type:\n")
table(merged_data$Annotation.x)  # Using .x suffix from the merge

# Frequency comparison between populations
cat("\nFrequency statistics by population:\n")
summary(merged_data[c("SH_urban_freq", "SH_rural_freq", "USA_freq")])



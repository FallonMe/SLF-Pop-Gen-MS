# Load required libraries
library(ggplot2)
library(dplyr)

# Define colors for each population
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

country_colors <- c(
  "USA" = "#3690c0",
  "CN" = "#c861a5"
)

habitat_colors <- c(
  "urban" = "#c861a5",
  "nonurban" = "#009E73"
)

SH_colors <- c(
  "Jinqiao" = "#3690c0",
  "Sheshan" = "#984ea3",
  "Tianma" = "#4daf4a"
)


# Set working directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")

# Load population data
popdata_SH <- read.csv("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/SH_SLF_MS_Info.csv", header = TRUE)  # Population metadata

# Load genetic data in .bed format
SLF_bed_SH <- read.pcadapt("./SH_final_annotated.bed", type = "bed")  # Load BED file for PCA

# Extract population identifiers for plotting
popbed_SH <- popdata_SH[, c("SLF_ID", "Location")]  # Keep only ID and location
poplist.names_SH <- popbed_SH[, 2]  # Extract population names (locations)

# Perform PCA with a high number of PCs to start
x_SH <- pcadapt(input = SLF_bed_SH, K = 10, method = "mahalanobis", ploidy = 2)  # Start with K=10 PCs to determine how many to retain

# Generate a scree plot to identify the optimal number of PCs
pdf("SH_pcadapt_Screeplot.pdf", width = 8, height = 6)  # Save scree plot as PDF
scree_SH <- plot(x_SH , option = "screeplot")  # Plot scree plot
dev.off()  # Close the PDF device

# Perform PCA with the optimal number of PCs based on the scree plot
x_SH <- pcadapt(input = SLF_bed_SH, K = 2)  # Retain 2 PCs based on the scree plot

# Extract PCA results for plotting
x2df_SH <- data.frame(x_SH$scores)  # PCA scores for K=2
x2df_SH$pops <- poplist.names_SH  # Add population information
x2df_SH$pops <- as.factor(x2df_SH$pops)  # Convert population to factor

# Plot PC1 vs PC2 with updated PCA (K=3) and legend for locations
pdf("SH_PC1_vs_PC2_K3_pcadapt.pdf", width = 10, height = 6)
ggplot(x2df_SH, aes(X1, X2, color = pops, fill = pops)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = SH_colors, name = "Location") + # Add legend for color-location mapping
  xlab("PC1") + ylab("PC2") +
  stat_ellipse() +
  theme(legend.position = "right") + 
  scale_fill_manual(values = SH_colors, name = "Location") 
dev.off()

# Proportion of variance explained
var_explained_SH <- (x_SH$singular.values)^2 / sum((x_SH$singular.values)^2)

# Run pcadapt with K = 2, using Mahalanobis distance and diploid organisms (ploidy=2)
slf.pcadapt <- pcadapt(input = SLF_bed_SH, K = 2, method = "mahalanobis", ploidy = 2)

# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space.
# "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.

names(slf.pcadapt)

# The Genomic Inflation Factor is stored in "gif" (more on that below)
# The (rescaled = stat/gif) test statistics are stored in "chi2.stat"
# The recalibrated p-values are stored in "pvalues"

plot(slf.pcadapt, option="manhattan")
plot(slf.pcadapt, option="qqplot")

# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.

slf.pcadapt$gif
#1.38

# GIF of 1=well calibrated, >1=liberal (too many small p-values), <1=conservative (too few small p-values)
# Note: GIFs > 2 indicate poor test calibration

hist(slf.pcadapt$pvalues, xlab="p-values", main=NULL, breaks=50)
# Note: these p-values are already rescaled with the GIF...

# Load annotation data
annot_combined <- read.csv("Final_SLF_Annotation_Info.csv")  

# Ensure SNP ID format is consistent
#annot_combined$SNP <- gsub("_[ACGT]$", "", annot_combined$SNP)  # Remove allele suffix if present

# False Discovery Rate (FDR) Correction using Q-values
qv <- qvalue::qvalue(slf.pcadapt$pvalues)
qv_adj <- qv$qvalues  # Extract the FDR-adjusted p-values

# Select only SNPs with FDR-adjusted p-values < 0.01
fdr_snps_filtered <- which(qv_adj < 0.01)  
fdr_snps <- annot_combined$SNP[fdr_snps_filtered]
fdr_pvals <- qv_adj[fdr_snps_filtered]  # Ensure correct p-values
#760 Significant SNPs out of 158335

# Create FDR DataFrame with Corrected Filtering
fdr_snps_df <- data.frame(SNP = fdr_snps, PValue = fdr_pvals, Method = "FDR")

# Ensure SNP IDs match the format in annot_combined
bonferroni_snps <- annot_combined$SNP[outliers_bonferroni]
fdr_snps <- annot_combined$SNP[fdr_snps_filtered]

# save with SNP ID and FDR-adjusted p-values
bonferroni_snps_df <- data.frame(SNP = bonferroni_snps, PValue = bonferroni_pvals, Method = "Bonferroni")

fdr_snps_df <- data.frame(SNP = fdr_snps, PValue = fdr_pvals, Method = "FDR")

# Merge FDR results with annotation using correct filtering
fdr_merged <- merge(fdr_snps_df, annot_combined, by = "SNP", all.x = TRUE)

fdr_merged$SNP <- gsub("_[ACGT]$", "", fdr_merged$SNP) 

# Check if incorrect values still persist
summary(fdr_merged$PValue)

# Save the fixed FDR results
write.csv(fdr_merged, "SH_pcadapt_SigFDR__results.csv", row.names = FALSE)

# Merge Bonferroni results with annotation
bonferroni_merged <- merge(bonferroni_snps_df, annot_combined, by = "SNP", all.x = TRUE)
write.csv(bonferroni_merged, "SH_pcadapt_Sigbonferroni_snps.csv", row.names = FALSE)

####genes with valid annotation and -log10(P) > 10
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load data
SH_data <- read.csv("SH_pcadapt_SigFDR__results.csv")

# Extract chromosome numbers
SH_data$chr_num <- as.numeric(gsub("chr_", "", SH_data$CHROM))

# Order data by chromosome and position
SH_data <- SH_data %>%
  arrange(chr_num, POS)

# Create alternating colors for chromosomes
SH_data$chr_color <- factor(SH_data$chr_num %% 2)

# Calculate chromosome boundaries for x-axis
chr_pos <- SH_data %>%
  group_by(CHROM) %>%
  summarize(
    chr_len = max(POS),
    chr_center = mean(POS)
  ) %>%
  arrange(chr_num = as.numeric(gsub("chr_", "", CHROM)))

# Create a cumulative position for each SNP
chr_cum_pos <- 0
SH_data$pos_cum <- 0

for (i in 1:nrow(chr_pos)) {
  # Filter SNPs for current chromosome
  chr_snps <- SH_data$CHROM == chr_pos$CHROM[i]
  
  # Assign cumulative position
  SH_data$pos_cum[chr_snps] <- SH_data$POS[chr_snps] + chr_cum_pos
  
  # Update chr_pos with cumulative position for chromosome center
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  
  # Update cumulative position for next chromosome
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5000000  # 5Mb gap between chromosomes
}
# Highlight: genes with valid annotation and -log10(P) > 5
highlight_snps <- SH_data %>%
  filter(!is.na(Gene_Name),
         Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name),
         -log10(PValue) > 5) %>%
  group_by(Gene_Name) %>%
  slice_min(order_by = PValue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Manhattan plot
p <- ggplot(SH_data, aes(x = pos_cum, y = -log10(PValue), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#CC6677")) +
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = NULL,
    y = "-log10(p-value)",
    title = "Selection Signals in SH SLF Population",
    subtitle = paste0("Yellow: ", nrow(highlight_snps), " genes above -log10(p) > 10")
  )

# Add yellow highlight points
p <- p + 
  geom_point(data = highlight_snps, 
             aes(x = pos_cum, y = -log10(PValue)), 
             color = "black", fill = "yellow", shape = 21, size = 2.2, stroke = 0.4) +
  geom_text_repel(data = highlight_snps,
                  aes(x = pos_cum, y = -log10(PValue), label = Gene_Name),
                  size = 3.2, fontface = "bold", color = "black",
                  box.padding = 0.6, point.padding = 0.5, max.overlaps = 40,
                  segment.size = 0.3, segment.color = "grey50")

# Show plot
print(p)

# Save the plot as PDF
ggsave("SLF_Top_SH_Manhattan_Plot.pdf", p, width = 12, height = 8, dpi = 300)
ggsave("SLF_Top_SH_Manhattan_Plot.png", p, width = 12, height = 8, dpi = 300)

# GO Term Analysis for SH selection candidates
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(goseq)

# Extract unique gene names using explicit dplyr reference
sig_genes_SH <- SH_data %>%
  filter(!is.na(Gene_Name)) %>%
  dplyr::select(Gene_Name) %>%
  distinct()

# Create a vector of significant gene names
sig_gene_list_SH <- sig_genes_SH$Gene_Name

# Extract all genes that have GO annotations as background
# Use the same gene_id_go_mapping from the USA analysis
gene_id_go_mapping <- read.csv("GeneID_GeneName_GOterms.csv")

all_genes <- gene_id_go_mapping %>%
  dplyr::select(GeneID, GeneName) %>%
  distinct()

# Get all unique gene identifiers
all_gene_ids <- unique(c(all_genes$GeneID, all_genes$GeneName))

# Create a binary vector where 1 = significant gene, 0 = not significant
gene_vector_SH <- numeric(length(all_gene_ids))
names(gene_vector_SH) <- all_gene_ids

# Populate the gene vector
for (i in 1:length(all_gene_ids)) {
  gene_id <- all_gene_ids[i]
  # Check if this gene is in our significant gene list
  gene_vector_SH[i] <- ifelse(gene_id %in% sig_gene_list_SH, 1, 0)
}

# Create gene2go dataframe
gene2go <- data.frame(
  gene = gene_id_go_mapping$GeneName,
  category = gene_id_go_mapping$GOterms
) %>%
  dplyr::distinct()

# Create bias values with sufficient uniqueness for goseq
set.seed(123)
bias <- rank(as.numeric(factor(names(gene_vector_SH)))) + runif(length(gene_vector_SH), min=0, max=0.1)
names(bias) <- names(gene_vector_SH)

# Create probability weighting function
pwf <- nullp(gene_vector_SH, bias.data = bias, plot.fit = FALSE)

# Run GO enrichment analysis
GO_results_SH <- goseq(pwf, gene2cat = gene2go, use_genes_without_cat = TRUE)

# Add FDR correction and process results
GO_results_all_SH <- GO_results_SH %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))

# Save all results
write.csv(GO_results_all_SH, "SH_Selection_GO_enrichment_all_results.csv")

# Now filter at a higher threshold to see if we get any results
GO_results_filtered_SH <- GO_results_all_SH %>%
  filter(FDR < 0.01) %>%  
  arrange(over_represented_pvalue)
# 0

# 9. Generate visualization even if no terms meet strict significance
# Select the top 20 most significant GO terms
top_GO_SH <- GO_results_filtered %>%
  arrange(over_represented_pvalue) %>%
  head(20)

# Create the bar plot
enrichment_plot <- ggplot(top_GO, aes(x = reorder(term, -over_represented_pvalue), 
                                      y = -log10(over_represented_pvalue), fill = ontology)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
  labs(title = "Over_represented GO Terms",
       x = NULL,
       y = "-log10(p-value)",
       fill = "Ontology") +
  theme_minimal()

# Save the plot
ggsave("USA_Selection_GO_top_terms.pdf", enrichment_plot, width = 12, height = 8)

# 10. Create a bubble plot with enrichment information
top_GO <- top_GO %>%
  mutate(enrichment_ratio = numDEInCat / 
           (numInCat * sum(gene_vector) / length(gene_vector)))

bubble_plot <- ggplot(top_GO, 
                      aes(x = enrichment_ratio, 
                          y = reorder(term, enrichment_ratio),
                          size = numDEInCat,
                          color = -log10(over_represented_pvalue))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "GO Term Enrichment in Selection Candidates",
       x = "Enrichment Ratio (>1 = enriched)",
       y = "GO Term",
       size = "Gene Count",
       color = "-log10(p-value)") +
  theme_minimal()

ggsave("USA_Selection_GO_enrichment_bubble_plot.pdf", bubble_plot, width = 12, height = 8)

# Print a summary
cat("Analysis complete. Top GO terms saved to files.\n")
cat("Number of terms with FDR < 0.01:", sum(GO_results_all_SH$FDR < 0.01), "\n")

save.image("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Pcadpt.RData")


#===============================================================================
# SPOTTED LANTERNFLY SELECTION ANALYSIS PIPELINE
#===============================================================================
# This script performs selection analysis on Spotted Lanternfly (SLF) genomic data
# using PCAdapt to identify loci under selection and GO term enrichment analysis 
# to understand the functional significance of these loci.
#===============================================================================

#-------------------------------------------------------------------------------
# PART 1: SETUP AND DATA LOADING
#-------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(pcadapt)
library(qvalue)

# Define color schemes for visualization
location_colors <- c(
  "Bronx" = "#252525", 
  "Brooklyn" = "#e6ab02",
  "Camden" = "#a65628",
  "Connecticut" = "#e41a1c",
  "Governor Island" = "#08589e",
  "Manhattan" = "#c861a5",
  "Staten Island" = "#3690c0", 
  "West Orange" = "#009E73"
)

country_colors <- c("USA" = "#3690c0")

habitat_colors <- c(
  "urban" = "#c861a5",
  "nonurban" = "#009E73"
)

# Set working directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")

# Load population metadata and genomic data
popdata_USA <- read.csv("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/USA_SLF_MS_Info.csv", header = TRUE)
SLF_bed_USA <- read.pcadapt("./USA_final_annotated.bed", type = "bed")

# Extract population identifiers for plotting
popbed_USA <- popdata_USA[, c("SLF_ID", "Location")]
poplist.names_USA <- popbed_USA[, 2]

#-------------------------------------------------------------------------------
# PART 2: INITIAL PCA EXPLORATION
#-------------------------------------------------------------------------------

# Perform PCA with a high number of PCs to evaluate population structure
x_USA <- pcadapt(input = SLF_bed_USA, K = 20, method = "mahalanobis", ploidy = 2)

# Generate scree plot to identify optimal number of PCs
pdf("USA_pcadapt_Screeplot.pdf", width = 8, height = 6)
scree_USA <- plot(x_USA, option = "screeplot")
dev.off()

# Perform PCA with optimal number of PCs (K=1 based on scree plot)
x_USA <- pcadapt(input = SLF_bed_USA, K = 1, method = "mahalanobis", ploidy = 2)

# Visualize PCA results by location
x2df_USA <- data.frame(x_USA$scores)
x2df_USA$pops <- poplist.names_USA
x2df_USA$pops <- as.factor(x2df_USA$pops)

# Create jitter plot for K=1
ggplot(x2df_USA, aes(x = pops, y = x_USA.scores, color = pops)) +
  geom_jitter(width = 0.2, size = 2) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  xlab("Location") +
  ylab("PC1") +
  theme(legend.position = "none")

# Save PCA visualization
pdf("USA_PC1_vs_PC2_K1_pcadapt.pdf", width = 10, height = 6)
ggplot(x2df_USA, aes(X1, X2, color = pops, fill = pops)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  xlab("PC1") + ylab("PC2") +
  stat_ellipse() +
  theme(legend.position = "right") + 
  scale_fill_manual(values = location_colors, name = "Location") 
dev.off()

# Calculate proportion of variance explained
var_explained_USA <- (x_USA$singular.values)^2 / sum((x_USA$singular.values)^2)

#-------------------------------------------------------------------------------
# PART 3: PCADAPT SELECTION ANALYSIS
#-------------------------------------------------------------------------------

# Run pcadapt with K = 5 for selection analysis
slf.pcadapt_USA <- pcadapt(input = SLF_bed_USA, K = 5, method = "mahalanobis", ploidy = 2)

# Examine PCAdapt output
names(slf.pcadapt_USA)

# Generate diagnostic plots
plot(slf.pcadapt_USA, option="manhattan")
plot(slf.pcadapt_USA, option="qqplot")

# Check genomic inflation factor (GIF)
slf.pcadapt_USA$gif
#1.3
# Note: GIF of 1=well calibrated, >1=liberal (too many small p-values), <1=conservative

# Histogram of p-values to examine distribution
hist(slf.pcadapt_USA$pvalues, xlab="p-values", main=NULL, breaks=50)

# Visualize updated PCA results
x2df_USA <- data.frame(slf.pcadapt_USA$scores)
x2df_USA$pops <- poplist.names_USA
x2df_USA$pops <- as.factor(x2df_USA$pops)

pdf("USA_PC1_vs_PC2_K1_pcadapt.pdf", width = 10, height = 6)
ggplot(x2df_USA, aes(X1, X2, color = pops, fill = pops)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  xlab("PC1") + ylab("PC2") +
  stat_ellipse() +
  theme(legend.position = "right") + 
  scale_fill_manual(values = location_colors, name = "Location") 
dev.off()

#-------------------------------------------------------------------------------
# PART 4: MULTI-K PCADAPT ANALYSIS
#-------------------------------------------------------------------------------

# Function to run pcadapt across multiple K values and identify consistent outliers
run_pcadapt_multiple_K <- function(bed_file, pop_data, K_values = c(2, 3, 4, 5), annot_data = NULL, 
                                   fdr_threshold = 0.01, output_dir = "./pcadapt_results/") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load required libraries
  require(pcadapt)
  require(ggplot2)
  require(qvalue)
  require(dplyr)
  require(tidyr)
  require(VennDiagram)
  
  # Load genetic data in .bed format if provided as string path
  if (is.character(bed_file)) {
    genetic_data <- read.pcadapt(bed_file, type = "bed")
  } else {
    genetic_data <- bed_file # Assume it's already loaded
  }
  
  # Extract population identifiers for plotting
  if (is.character(pop_data)) {
    popdata <- read.csv(pop_data, header = TRUE)
  } else {
    popdata <- pop_data # Assume it's already loaded
  }
  
  popbed <- popdata[, c("SLF_ID", "Location")]  # Adjust column names if needed
  poplist.names <- popbed[, 2]  # Extract population names (locations)
  
  # Store results for each K value
  results_list <- list()
  candidate_loci_list <- list()
  
  # Assuming you have a color vector defined - if not, create one
  if (!exists("location_colors")) {
    n_locations <- length(unique(poplist.names))
    location_colors <- rainbow(n_locations)
    names(location_colors) <- unique(poplist.names)
  }
  
  # Run pcadapt for each K value
  for (k in K_values) {
    cat(paste0("Running pcadapt with K = ", k, "\n"))
    
    # Run pcadapt
    pcadapt_results <- pcadapt(input = genetic_data, K = k, method = "mahalanobis", ploidy = 2)
    
    # Store results
    results_list[[paste0("K", k)]] <- pcadapt_results
    
    # Extract PCA results for plotting
    pca_df <- data.frame(pcadapt_results$scores)
    pca_df$pops <- poplist.names
    pca_df$pops <- as.factor(pca_df$pops)
    
    # Generate plot filenames
    pca_plot_file <- paste0(output_dir, "PCA_plot_K", k, ".pdf")
    manhattan_plot_file <- paste0(output_dir, "Manhattan_plot_K", k, ".pdf")
    qqplot_file <- paste0(output_dir, "QQ_plot_K", k, ".pdf")
    
    # Plot PCA results
    pdf(pca_plot_file, width = 10, height = 6)
    if (k >= 2) {
      p <- ggplot(pca_df, aes(X1, X2, color = pops, fill = pops)) +
        geom_point(size = 3) +
        theme_bw() +
        scale_color_manual(values = location_colors, name = "Location") +
        xlab("PC1") + ylab("PC2") +
        stat_ellipse() +
        theme(legend.position = "right") +
        scale_fill_manual(values = location_colors, name = "Location") +
        ggtitle(paste0("PCA Plot - K = ", k))
      print(p)
    } else {
      # For K=1, plot differently since there's only one PC
      p <- ggplot(pca_df, aes(x = pops, y = X1, color = pops)) +
        geom_jitter(width = 0.2, size = 2) +
        theme_bw() +
        scale_color_manual(values = location_colors, name = "Location") +
        xlab("Location") +
        ylab("PC1") +
        theme(legend.position = "right") +
        ggtitle(paste0("PCA Plot - K = ", k))
      print(p)
    }
    dev.off()
    
    # Plot diagnostics
    pdf(manhattan_plot_file, width = 10, height = 6)
    plot(pcadapt_results, option = "manhattan")
    dev.off()
    
    pdf(qqplot_file, width = 10, height = 6)
    plot(pcadapt_results, option = "qqplot")
    dev.off()
    
    # Calculate FDR-adjusted p-values
    qv <- qvalue(pcadapt_results$pvalues)
    qv_adj <- qv$qvalues
    
    # Select SNPs with FDR-adjusted p-values < threshold
    fdr_snps_indices <- which(qv_adj < fdr_threshold)
    
    # Store candidate loci information
    if (!is.null(annot_data)) {
      fdr_snps <- annot_data$SNP[fdr_snps_indices]
      fdr_pvals <- qv_adj[fdr_snps_indices]
      
      # Create a dataframe with candidate SNPs
      candidate_df <- data.frame(
        SNP = fdr_snps,
        p_value = pcadapt_results$pvalues[fdr_snps_indices],
        q_value = fdr_pvals
      )
      
      # Save candidate SNPs to file
      write.csv(candidate_df, 
                file = paste0(output_dir, "candidate_loci_K", k, ".csv"),
                row.names = FALSE)
      
      # Store candidate loci for overlap analysis
      candidate_loci_list[[paste0("K", k)]] <- fdr_snps
    } else {
      # If no annotation data is provided, just save the indices
      candidate_loci_list[[paste0("K", k)]] <- fdr_snps_indices
      write.csv(data.frame(SNP_index = fdr_snps_indices), 
                file = paste0(output_dir, "candidate_loci_indices_K", k, ".csv"),
                row.names = FALSE)
    }
    
    # Print some summary statistics
    cat(paste0("K = ", k, " summary:\n"))
    cat(paste0("  Genomic Inflation Factor (GIF): ", round(pcadapt_results$gif, 3), "\n"))
    cat(paste0("  Number of candidate loci: ", length(fdr_snps_indices), "\n\n"))
  }
  
  # Find overlapping candidate loci between different K values
  if (length(K_values) > 1) {
    # Create Venn diagram for overlapping candidates
    venn_file <- paste0(output_dir, "venn_diagram_overlapping_loci.pdf")
    
    venn_colors <- c("red", "blue", "green", "purple", "orange")[1:length(K_values)]
    venn_labels <- paste0("K=", K_values)
    
    pdf(venn_file, width = 8, height = 8)
    venn.diagram(
      x = candidate_loci_list,
      category.names = venn_labels,
      filename = NULL,
      output = TRUE,
      imagetype = "pdf",
      scaled = TRUE,
      col = "black",
      fill = venn_colors,
      alpha = 0.5,
      main = "Overlapping candidate loci"
    )
    dev.off()
    
    # Find SNPs that overlap in all K values
    overlap_all <- Reduce(intersect, candidate_loci_list)
    
    # Save overlapping SNPs
    if (!is.null(annot_data) && length(overlap_all) > 0) {
      if (is.character(overlap_all[1])) {
        # If we have SNP IDs
        overlap_df <- data.frame(SNP = overlap_all)
      } else {
        # If we only have indices
        overlap_df <- data.frame(SNP_index = overlap_all)
        if (!is.null(annot_data)) {
          overlap_df$SNP <- annot_data$SNP[overlap_all]
        }
      }
      
      write.csv(overlap_df, 
                file = paste0(output_dir, "overlapping_candidate_loci_all_K.csv"),
                row.names = FALSE)
    }
    
    # Create a table showing the number of SNPs for each K and the overlap
    overlap_summary <- data.frame(
      K_value = paste0("K=", K_values),
      Candidates = sapply(candidate_loci_list, length)
    )
    
    # Add pairwise overlaps
    if (length(K_values) > 1) {
      for (i in 1:(length(K_values)-1)) {
        for (j in (i+1):length(K_values)) {
          k1 <- paste0("K", K_values[i])
          k2 <- paste0("K", K_values[j])
          overlap_label <- paste0("Overlap_", K_values[i], "_", K_values[j])
          overlap_count <- length(intersect(candidate_loci_list[[k1]], candidate_loci_list[[k2]]))
          overlap_summary[[overlap_label]] <- c(rep(NA, i-1), overlap_count, rep(NA, length(K_values)-i-1))
        }
      }
    }
    
    overlap_summary$Overlap_all <- c(length(overlap_all), rep(NA, length(K_values)-1))
    
    write.csv(overlap_summary, 
              file = paste0(output_dir, "overlap_summary.csv"),
              row.names = FALSE)
    
    cat("Overlap analysis complete.\n")
    cat(paste0("Number of loci appearing in all K values: ", length(overlap_all), "\n"))
  }
  
  # Return results
  return(list(
    pcadapt_results = results_list,
    candidate_loci = candidate_loci_list,
    overlap = if(exists("overlap_all")) overlap_all else NULL
  ))
}

#-------------------------------------------------------------------------------
# PART 5: EXECUTE MULTI-K PCADAPT ANALYSIS
#-------------------------------------------------------------------------------

# Load additional required packages
library(VennDiagram)
library(grid)
library(gridExtra)

# Execute the multi-K pcadapt analysis
results <- run_pcadapt_multiple_K(
  bed_file = "./USA_final_annotated.bed",
  pop_data = "~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/USA_SLF_MS_Info.csv",
  K_values = c(2, 3, 4, 5),
  annot_data = annot_combined,
  fdr_threshold = 0.01,
  output_dir = "./pcadapt_results/"
)

# Compare number of candidates at different K values
k2_candidates <- read.csv("./pcadapt_results/candidate_loci_K2.csv")
k3_candidates <- read.csv("./pcadapt_results/candidate_loci_K3.csv")
k4_candidates <- read.csv("./pcadapt_results/candidate_loci_K4.csv")
k5_candidates <- read.csv("./pcadapt_results/candidate_loci_K5.csv")

# Number of candidate loci at each K value
nrow(k2_candidates)
nrow(k3_candidates)
nrow(k4_candidates)
nrow(k5_candidates)

# Verify all q-values are below 0.01 in each dataset
all(k2_candidates$q_value < 0.01)
all(k3_candidates$q_value < 0.01)
all(k4_candidates$q_value < 0.01)
all(k5_candidates$q_value < 0.01)

# Get summary statistics of q-values
summary(k2_candidates$q_value)
summary(k3_candidates$q_value)
summary(k4_candidates$q_value)
summary(k5_candidates$q_value)

# Create Venn diagram visualization of overlapping candidates
venn <- venn.diagram(
  x = candidate_lists,
  category.names = paste0("K=", c(2, 3, 4, 5)),
  filename = NULL,
  output = TRUE,
  scaled = TRUE,
  col = "black",
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  main = "Overlapping candidate loci"
)

# Save Venn diagram
pdf("./custom_venn_diagram.pdf", width = 8, height = 8)
grid.newpage()
grid.draw(venn)
dev.off()

#-------------------------------------------------------------------------------
# PART 6: PROCESS AND MERGE CANDIDATE LOCI
#-------------------------------------------------------------------------------

# Load the overlapping SNPs across all K values
overlap_all_k <- read.csv("./pcadapt_results/overlapping_candidate_loci_all_K.csv")
head(overlap_all_k)
nrow(overlap_all_k)  # Number of SNPs identified across all K values - 7343 out of 158335

# Load annotation data and merge with candidate SNPs
annot_combined <- read.csv("Final_SLF_Annotation_Info.csv")  
fdr_USA_Only <- merge(overlap_all_k, annot_combined, by = "SNP", all.x = TRUE)
write.csv(fdr_USA_Only, "USA_Only_pcadapt_SigFDR_results.csv", row.names = FALSE)

# Create a comprehensive dataframe with statistics from all K values
comprehensive_df <- data.frame(SNP = overlap_all_k$SNP)

# Add statistics for each K value
for (k in c(2, 3, 4, 5)) {
  # Get the candidate data for this K
  k_data <- get(paste0("k", k, "_candidates"))
  
  # Match the overlapping SNPs with their statistics
  matched_indices <- match(comprehensive_df$SNP, k_data$SNP)
  
  # Add p-values and q-values to the comprehensive dataframe
  comprehensive_df[paste0("p_value_K", k)] <- k_data$p_value[matched_indices]
  comprehensive_df[paste0("q_value_K", k)] <- k_data$q_value[matched_indices]
}

# Add annotation information
if ("gene" %in% colnames(annot_combined) || "function" %in% colnames(annot_combined)) {
  matched_indices <- match(comprehensive_df$SNP, annot_combined$SNP)
  
  # Add relevant annotation columns
  annotation_cols <- c()
  if ("gene" %in% colnames(annot_combined)) annotation_cols <- c(annotation_cols, "gene")
  if ("function" %in% colnames(annot_combined)) annotation_cols <- c(annotation_cols, "function")
  if ("scaffold" %in% colnames(annot_combined)) annotation_cols <- c(annotation_cols, "scaffold")
  if ("position" %in% colnames(annot_combined)) annotation_cols <- c(annotation_cols, "position")
  
  for (col in annotation_cols) {
    if (col %in% colnames(annot_combined)) {
      comprehensive_df[col] <- annot_combined[matched_indices, col]
    }
  }
}

# Calculate average p-value and q-value across all K values
comprehensive_df$avg_p_value <- rowMeans(comprehensive_df[, grep("^p_value_K", colnames(comprehensive_df))])
comprehensive_df$avg_q_value <- rowMeans(comprehensive_df[, grep("^q_value_K", colnames(comprehensive_df))])

# Sort by average q-value
comprehensive_df <- comprehensive_df[order(comprehensive_df$avg_q_value), ]

# Merge with annotation data and save
comprehensive_fdr_USA_Only <- merge(comprehensive_df, annot_combined, by = "SNP", all.x = TRUE)
write.csv(comprehensive_fdr_USA_Only, "comprehensive_fdr_USA_Only_SigFDR.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# PART 7: MANHATTAN PLOT VISUALIZATION
#-------------------------------------------------------------------------------

library(ggrepel)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection")
# Load data
snp_data <- read.csv("comprehensive_fdr_USA_Only_SigFDR.csv")

# Extract chromosome numbers
snp_data$chr_num <- as.numeric(gsub("chr_", "", snp_data$CHROM))

# Order data by chromosome and position
snp_data <- snp_data %>%
  arrange(chr_num, POS)

# Create alternating colors for chromosomes
snp_data$chr_color <- factor(snp_data$chr_num %% 2)

# Calculate chromosome boundaries for x-axis
chr_pos <- snp_data %>%
  group_by(CHROM) %>%
  summarize(
    chr_len = max(POS),
    chr_center = mean(POS)
  ) %>%
  arrange(chr_num = as.numeric(gsub("chr_", "", CHROM)))

# Create a cumulative position for each SNP
chr_cum_pos <- 0
snp_data$pos_cum <- 0

for (i in 1:nrow(chr_pos)) {
  # Filter SNPs for current chromosome
  chr_snps <- snp_data$CHROM == chr_pos$CHROM[i]
  
  # Assign cumulative position
  snp_data$pos_cum[chr_snps] <- snp_data$POS[chr_snps] + chr_cum_pos
  
  # Update chr_pos with cumulative position for chromosome center
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  
  # Update cumulative position for next chromosome
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5000000  # 5Mb gap between chromosomes
}

# Define p-value threshold for highlighting genes
p_value_threshold <- 1e-10
custom_threshold <- -log10(p_value_threshold)

# Mark significant SNPs
snp_data$is_significant <- snp_data$avg_p_value < p_value_threshold

# Find genes with valid names above threshold
valid_genes <- snp_data %>%
  filter(
    avg_p_value < p_value_threshold,
    !grepl("^LDel_", Gene_Name),
    !grepl("-", Gene_Name),
    !is.na(Gene_Name),
    Gene_Name != ""
  )

# For each valid gene, keep only the most significant SNP
significant_snps <- valid_genes %>%
  group_by(Gene_Name) %>%
  filter(avg_p_value == min(avg_p_value)) %>%
  ungroup()

# Add manually specified genes of interest
additional_genes <- c("MDH1")  
for (gene in additional_genes) {
  if (!gene %in% significant_snps$Gene_Name) {
    # Find the best SNP for this gene
    gene_snps <- snp_data %>% 
      filter(Gene_Name == gene) %>%
      arrange(avg_p_value) %>%
      head(1)
    
    if (nrow(gene_snps) > 0) {
      significant_snps <- rbind(significant_snps, gene_snps)
    }
  }
}

# Create Manhattan plot
p <- ggplot(snp_data, aes(x = pos_cum, y = -log10(avg_p_value), color = chr_color)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#CC6677")) +
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Chromosome",
    y = "-log10(p-value)",
    title = "Manhattan Plot of SLF Selection Signals"
  )

# Add labels for valid genes above threshold
p <- p + geom_text_repel(
  data = significant_snps,
  aes(label = Gene_Name),
  size = 3.5,
  fontface = "bold",
  color = "black",
  box.padding = 0.5,
  point.padding = 0.5,
  max.overlaps = 25,
  segment.size = 0.2,
  segment.color = "grey50"
)

# Add points with different appearance for highlighted SNPs
p <- p + geom_point(
  data = significant_snps,
  aes(x = pos_cum, y = -log10(avg_p_value)),
  color = "black",
  size = 1.5,
  shape = 21,
  fill = "yellow",
  alpha = 0.8
)

# Save the plot
ggsave("SLF_USA_Manhattan_Plot.pdf", p, width = 12, height = 7, dpi = 300)
ggsave("SLF_USA_Manhattan_Plot.png", p, width = 12, height = 7, dpi = 300)
print(p)


#-------------------------------------------------------------------------------
# PART 8: GO TERM ENRICHMENT ANALYSIS
#-------------------------------------------------------------------------------

# Create a binary vector where 1 = significant gene, 0 = not significant
gene_vector <- numeric(length(all_gene_ids))
names(gene_vector) <- all_gene_ids

# Populate the gene vector
for (i in 1:length(all_gene_ids)) {
  gene_id <- all_gene_ids[i]
  # Check if this gene is in our significant gene list
  gene_vector[i] <- ifelse(gene_id %in% sig_gene_list, 1, 0)
}

# Create gene2go dataframe for goseq
gene2go <- data.frame(
  gene = gene_id_go_mapping$GeneName,
  category = gene_id_go_mapping$GOterms
) %>%
  dplyr::distinct()

# Create bias values with sufficient uniqueness for goseq
# Using rank of gene IDs plus small random noise for bias
set.seed(123) # For reproducibility
bias <- rank(as.numeric(factor(names(gene_vector)))) + runif(length(gene_vector), min=0, max=0.1)
names(bias) <- names(gene_vector)

# Create probability weighting function with permissive parameters
pwf <- nullp(gene_vector, bias.data = bias, plot.fit = FALSE)

# Run GO enrichment analysis with higher p-value threshold
GO_results <- goseq(pwf, gene2cat = gene2go, use_genes_without_cat = TRUE)

# Add FDR correction and process results
GO_results_all <- GO_results %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))

# Save all results
write.csv(GO_results_all, "USAOnly_Selection_GO_enrichment_all_results.csv", row.names = FALSE)

# Filter at a higher threshold to capture more results
GO_results_filtered <- GO_results_all %>%
  filter(FDR < 0.05) %>%  # More lenient threshold
  arrange(over_represented_pvalue)

# Generate visualization of GO enrichment results
# Select the top significant GO terms
top_GO <- GO_results_filtered %>%
  arrange(over_represented_pvalue) %>%
  head(31)

# Create bar plot of enriched GO terms
enrichment_plot <- ggplot(top_GO, aes(x = reorder(term, -over_represented_pvalue), 
                                      y = -log10(over_represented_pvalue), fill = ontology)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
  labs(title = "Over_represented GO Terms",
       x = NULL,
       y = "-log10(p-value)",
       fill = "Ontology") +
  theme_minimal()

# Save the GO term enrichment plot
ggsave("USA_Selection_GO_top_terms.pdf", enrichment_plot, width = 12, height = 8)

# Create a bubble plot with enrichment information
top_GO <- top_GO %>%
  mutate(enrichment_ratio = numDEInCat / 
           (numInCat * sum(gene_vector) / length(gene_vector)))

bubble_plot <- ggplot(top_GO, 
                      aes(x = enrichment_ratio, 
                          y = reorder(term, enrichment_ratio),
                          size = numDEInCat,
                          color = -log10(over_represented_pvalue))) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "GO Term Enrichment in Selection Candidates",
       x = "Enrichment Ratio (>1 = enriched)",
       y = "GO Term",
       size = "Gene Count",
       color = "-log10(p-value)") +
  theme_minimal()

# Save the bubble plot
ggsave("USA_Selection_GO_enrichment_bubble_plot.pdf", bubble_plot, width = 12, height = 8)

# Print summary of GO enrichment analysis
cat("Analysis complete. Top GO terms saved to files.\n")
cat("Number of terms with FDR < 0.05:", sum(GO_results_all$FDR < 0.05), "\n")

# Save R workspace
save.image("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Pcadpt.RData")

#-------------------------------------------------------------------------------
# PART 9: MANHATTAN PLOT WITH OVERLAP HIGHLIGHTING
#-------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load data
snp_data <- read.csv("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/comprehensive_fdr_USA_Only_SigFDR.csv")

# Extract chromosome numbers
snp_data$chr_num <- as.numeric(gsub("chr_", "", snp_data$CHROM))

# Order data by chromosome and position
snp_data <- snp_data %>%
  arrange(chr_num, POS)

# Create alternating colors for chromosomes
snp_data$chr_color <- factor(snp_data$chr_num %% 2)

# Calculate chromosome boundaries for x-axis
chr_pos <- snp_data %>%
  group_by(CHROM) %>%
  summarize(
    chr_len = max(POS),
    chr_center = mean(POS),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

# Initialize center_cum column
chr_pos$center_cum <- 0

# Create a cumulative position for each SNP
chr_cum_pos <- 0
snp_data$pos_cum <- 0
for (i in 1:nrow(chr_pos)) {
  # Filter SNPs for current chromosome
  chr_snps <- snp_data$CHROM == chr_pos$CHROM[i]
  
  # Assign cumulative position
  snp_data$pos_cum[chr_snps] <- snp_data$POS[chr_snps] + chr_cum_pos
  
  # Update chr_pos with cumulative position for chromosome center
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  
  # Update cumulative position for next chromosome
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5000000  # 5Mb gap between chromosomes
}

# Create unique identifier for SNPs
snp_data$SNP_ID <- paste(snp_data$CHROM, snp_data$POS, sep="_")

# Process overlap SNPs
if(is.data.frame(prda_pcadapt)) {
  overlap_ids <- paste(prda_pcadapt$CHROM, prda_pcadapt$POS, sep="_")
} else {
  # If prda_pcadapt is a vector of already formatted IDs, use directly
  overlap_ids <- prda_pcadapt
}

# Mark overlap SNPs
snp_data$is_overlap <- snp_data$SNP_ID %in% overlap_ids

# Find genes in the overlap set with valid names
overlap_genes <- snp_data %>%
  filter(
    is_overlap == TRUE,
    !grepl("^LDel_", Gene_Name),
    !grepl("-", Gene_Name),
    !is.na(Gene_Name),
    Gene_Name != ""
  )

# Extract unique gene names from overlap
unique(overlap_genes$Gene_Name)
write.csv(overlap_genes, "overlap_genes185_unqiuegenenames.csv")

# Use all overlap genes without filtering by significance
significant_snps <- overlap_genes

# Count how many overlap SNPs were found
overlap_count <- sum(snp_data$is_overlap, na.rm = TRUE)
print(paste("Found", overlap_count, "overlap SNPs out of expected 185"))

#-------------------------------------------------------------------------------
# PART 10: MANHATTAN PLOT VARIANTS
#-------------------------------------------------------------------------------
# Identify unique functional genes (not LDel, not empty, not containing '-')
highlight_snps_usa <- significant_snps %>%
  filter(!is.na(Gene_Name),
         Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  group_by(Gene_Name) %>%
  slice_min(order_by = avg_p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(snp_data, aes(x = pos_cum, y = -log10(as.numeric(avg_p_value)), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#3690c0")) +
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  ) +
  labs(
    x = "Chromosome",
    y = "-log10(p-value)",
    title = "Manhattan Plot of SLF Selection Signals",
    subtitle = paste0("All points represent SNPs significant at FDR < 0.01; pink highlights show ", 
                      sum(snp_data$is_overlap), " SNPs identified by both pRDA and PCAdapt")
  )

# IMPORTANT: Add highlighted points AFTER the main plot is created
# Make sure the fill color is applied correctly and alpha is appropriate
p <- p + geom_point(
  data = snp_data %>% filter(is_overlap == TRUE),
  aes(x = pos_cum, y = -log10(as.numeric(avg_p_value))),
  color = "black",         # outline color
  fill = "#c861a5",        # fill color (pink)
  size = 2,                # make slightly larger than regular points
  shape = 21,              # filled circle with outline
  stroke = 0.4,            # outline thickness
  alpha = 1                # full opacity to ensure visibility
)

p <- p + geom_text_repel(data = highlight_snps_usa,
                aes(x = pos_cum, y = -log10(avg_p_value), label = Gene_Name),
                size = 3.2, fontface = "bold", color = "black",
                box.padding = 0.6, point.padding = 0.5, max.overlaps = 40,
                segment.size = 0.3, segment.color = "grey50")


# Update the subtitle to match the actual color used
p <- p + labs(
  subtitle = paste0("All points represent SNPs significant at FDR < 0.01; pink highlights show ", 
                    sum(snp_data$is_overlap), " SNPs identified by both pRDA and PCAdapt")
)

print(p)

# Save the plot with alternating chromosome colors
ggsave("SLF_USA_Manhattan_Plot_Overlap.pdf", p, width = 14, height = 8, dpi = 300)
ggsave("SLF_USA_Manhattan_Plot_Overlap.png", p, width = 14, height = 8, dpi = 300)










######### GO Term Analysis for overlap and method-specific SNPs #######
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(goseq)  # For GO enrichment analysis
library(ggrepel)
library(VennDiagram)
library(ggVennDiagram)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection")

# Load gene-GO term mapping
gene_id_go_mapping <- read.csv("GeneID_GeneName_GOterms.csv")

# Function to perform GO enrichment analysis for a specific gene set
perform_go_analysis <- function(annotation_df, file_prefix, set_name) {
  # Extract gene names from the annotation dataframe
  if("Gene_Name" %in% colnames(annotation_df)) {
    gene_col <- "Gene_Name"
  } else if("GeneName" %in% colnames(annotation_df)) {
    gene_col <- "GeneName"
  } else {
    cat("Warning: No gene name column found in", set_name, "dataset\n")
    return(NULL)
  }
  
  # Extract unique gene names
  sig_genes <- annotation_df %>%
    filter(!is.na(!!sym(gene_col))) %>%
    select(!!sym(gene_col)) %>%
    distinct()
  
  # Create a vector of significant gene names
  sig_gene_list <- sig_genes[[gene_col]]
  
  # Print the number of genes with annotations
  cat("Number of genes with annotations in", set_name, ":", length(sig_gene_list), "\n")
  
  # Extract all genes that have GO annotations as background
  all_genes <- gene_id_go_mapping %>%
    select(GeneID, GeneName) %>%
    distinct()
  
  # Get all unique gene identifiers as background
  all_gene_ids <- unique(c(all_genes$GeneID, all_genes$GeneName))
  
  # Create a binary vector where 1 = gene from set, 0 = not in set
  gene_vector <- numeric(length(all_gene_ids))
  names(gene_vector) <- all_gene_ids
  
  # Populate the gene vector
  for (i in 1:length(all_gene_ids)) {
    gene_id <- all_gene_ids[i]
    gene_vector[i] <- ifelse(gene_id %in% sig_gene_list, 1, 0)
  }
  
  # Create the gene2go dataframe for mapping genes to GO terms
  gene2go <- data.frame(
    gene = gene_id_go_mapping$GeneName,
    category = gene_id_go_mapping$GOterms
  ) %>%
    dplyr::distinct()
  
  # Create bias values for goseq
  set.seed(123) # For reproducibility
  bias <- rank(as.numeric(factor(names(gene_vector)))) + runif(length(gene_vector), min=0, max=0.1)
  names(bias) <- names(gene_vector)
  
  # Create probability weighting function with permissive parameters
  pwf <- nullp(gene_vector, bias.data = bias, plot.fit = FALSE)
  
  # Run GO enrichment analysis with a higher p-value threshold initially
  GO_results <- goseq(pwf, gene2cat = gene2go, use_genes_without_cat = TRUE)
  
  # Add FDR correction and process results
  GO_results_all <- GO_results %>%
    as.data.frame() %>%
    mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))
  
  GO_results_all <- GO_results_all %>%
    filter(!is.na(term) & term != "NA")
  
  # Save all results
  write.csv(GO_results_all, paste0(file_prefix, "_GO_enrichment_all_results.csv"), row.names = FALSE)
  
  # Filter at a higher threshold to see if we get any results
  GO_results_filtered <- GO_results_all %>%
    filter(FDR < 0.05) %>%  # Using threshold of FDR < 0.05
    arrange(over_represented_pvalue)
  
  # Print the number of significant terms
  cat("Number of GO terms with FDR < 0.05 in", set_name, ":", nrow(GO_results_filtered), "\n")
  
  # Define top_GO for plotting - either the significant terms or top 15 by p-value
  top_GO <- if(nrow(GO_results_filtered) > 0) {
    GO_results_filtered %>% 
      head(min(15, nrow(GO_results_filtered)))
  } else {
    # If no significant terms, use the top 15 from all results
    GO_results_all %>%
      arrange(over_represented_pvalue) %>%
      head(15)
  }
  
  # Only create plots if we have data
  if(nrow(top_GO) > 0) {
    # Create the bar plot
    enrichment_plot <- ggplot(top_GO, aes(x = reorder(term, -over_represented_pvalue), 
                                          y = -log10(over_represented_pvalue), fill = ontology)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
      labs(title = paste("GO Terms Enrichment for", set_name),
           subtitle = paste0(nrow(GO_results_filtered), " terms with FDR < 0.05"),
           x = NULL,
           y = "-log10(p-value)",
           fill = "Ontology") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
    
    # Save the plot
    ggsave(paste0(file_prefix, "_GO_top_terms.pdf"), enrichment_plot, width = 12, height = 8)
    
    # Create a bubble plot with enrichment information
    top_GO <- top_GO %>%
      mutate(enrichment_ratio = numDEInCat / 
               (numInCat * sum(gene_vector) / length(gene_vector)))
    
    bubble_plot <- ggplot(top_GO, 
                          aes(x = enrichment_ratio, 
                              y = reorder(term, enrichment_ratio),
                              size = numDEInCat,
                              color = -log10(over_represented_pvalue))) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      labs(title = paste("GO Term Enrichment in", set_name),
           x = "Enrichment Ratio (>1 = enriched)",
           y = "GO Term",
           size = "Gene Count",
           color = "-log10(p-value)") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
    
    ggsave(paste0(file_prefix, "_GO_enrichment_bubble_plot.pdf"), bubble_plot, width = 12, height = 8)
  }
  
  # Create a table with detailed information for the top GO terms
  top_GO_table <- top_GO %>%
    select(term, ontology, numDEInCat, numInCat, over_represented_pvalue, FDR, enrichment_ratio) %>%
    arrange(over_represented_pvalue)
  
  # Save the detailed table
  write.csv(top_GO_table, paste0(file_prefix, "_top_GO_terms_details.csv"), row.names = FALSE)
  
  # Return all results for comparison, but mark which ones are significant
  return(GO_results_all %>% mutate(is_significant = FDR < 0.05))
}

# Perform GO analysis for each set
overlap_go <- perform_go_analysis(overlap_annotation, "overlap", "Overlap SNPs")
prda_only_go <- perform_go_analysis(prda_only_annotation, "prda_only", "pRDA-only SNPs")
pcadapt_only_go <- perform_go_analysis(pcadapt_only_annotation, "pcadapt_only", "PCAdapt-only SNPs")

# Create a list of unique genes for each category to verify the overlap
extract_genes <- function(annotation_df) {
  if("Gene_Name" %in% colnames(annotation_df)) {
    gene_col <- "Gene_Name"
  } else if("GeneName" %in% colnames(annotation_df)) {
    gene_col <- "GeneName"
  } else {
    return(character(0))
  }
  
  genes <- annotation_df %>%
    filter(!is.na(!!sym(gene_col))) %>%
    select(!!sym(gene_col)) %>%
    distinct() %>%
    pull()
  
  return(genes)
}

overlap_genes <- extract_genes(overlap_annotation)
prda_only_genes <- extract_genes(prda_only_annotation)
pcadapt_only_genes <- extract_genes(pcadapt_only_annotation)

# Create a summary of gene counts
gene_summary <- data.frame(
  Category = c("Overlap SNPs", "pRDA-only SNPs", "PCAdapt-only SNPs"),
  Gene_Count = c(length(overlap_genes), length(prda_only_genes), length(pcadapt_only_genes))
)

# Save gene lists
write.csv(data.frame(Gene = overlap_genes), "overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = prda_only_genes), "prda_only_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pcadapt_only_genes), "pcadapt_only_genes.csv", row.names = FALSE)
write.csv(gene_summary, "gene_count_summary.csv", row.names = FALSE)

# Extract the significant GO terms from each analysis and remove NAs
overlap_sig_terms <- overlap_go %>% 
  filter(is_significant) %>% 
  pull(term) %>%
  .[!is.na(.)]  # Remove NAs from the vector

prda_sig_terms <- prda_only_go %>% 
  filter(is_significant) %>% 
  pull(term) %>%
  .[!is.na(.)]  # Remove NAs from the vector

pcadapt_sig_terms <- pcadapt_only_go %>% 
  filter(is_significant) %>% 
  pull(term) %>%
  .[!is.na(.)]  # Remove NAs from the vector

# Print the length of each list to check
cat("Number of overlap GO terms with FDR < 0.05:", length(overlap_sig_terms), "\n")
cat("Number of PCAdapt-only GO terms with FDR < 0.05:", length(pcadapt_sig_terms), "\n")
cat("Number of pRDA-only GO terms with FDR < 0.05:", length(prda_sig_terms), "\n")

# Compare top GO terms across categories - only if all results are available
if(!is.null(overlap_go) && !is.null(prda_only_go) && !is.null(pcadapt_only_go)) {
  # Create a combined dataset of significant terms
  sig_terms <- bind_rows(
    overlap_go %>% filter(is_significant) %>% mutate(Category = "Overlap"),
    prda_only_go %>% filter(is_significant) %>% mutate(Category = "pRDA-only"),
    pcadapt_only_go %>% filter(is_significant) %>% mutate(Category = "PCAdapt-only")
  )
  
  # If we have significant terms, create visualizations
  if(nrow(sig_terms) > 0) {
    # Create a visualization comparing significant terms
    comparison_plot <- ggplot(sig_terms, 
                              aes(x = reorder(paste(term, " (", Category, ")", sep=""), 
                                              -over_represented_pvalue),
                                  y = -log10(over_represented_pvalue), 
                                  fill = Category)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("Overlap" = "#FFA500", 
                                   "pRDA-only" = "#4682B4", 
                                   "PCAdapt-only" = "#228B22")) +
      labs(title = "Comparison of Significant GO Terms Across Categories",
           subtitle = "GO terms with FDR < 0.05",
           x = NULL,
           y = "-log10(p-value)",
           fill = "Category") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
    ggsave("GO_terms_category_comparison.pdf", comparison_plot, width = 14, height = 10)
    
    # Create a dotplot showing significant terms
    dotplot <- ggplot(sig_terms, 
                      aes(x = Category, y = term, 
                          size = numDEInCat, 
                          color = -log10(over_represented_pvalue))) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = "Significant GO Terms Across Categories",
           subtitle = "GO terms with FDR < 0.05",
           x = NULL,
           y = NULL,
           size = "Gene Count",
           color = "-log10(p-value)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("GO_terms_dotplot_comparison.pdf", dotplot, width = 12, height = 8)
  }
  
  # If we don't have enough significant terms, create a visualization of the top terms regardless of significance
  if(nrow(sig_terms) < 10) {
    # Get top terms by p-value from each category
    top_terms <- bind_rows(
      overlap_go %>% arrange(over_represented_pvalue) %>% head(10) %>% mutate(Category = "Overlap"),
      prda_only_go %>% arrange(over_represented_pvalue) %>% head(10) %>% mutate(Category = "pRDA-only"),
      pcadapt_only_go %>% arrange(over_represented_pvalue) %>% head(10) %>% mutate(Category = "PCAdapt-only")
    )
    
    # Create dotplot of top terms
    top_dotplot <- ggplot(top_terms, 
                          aes(x = Category, y = term, 
                              size = numDEInCat, 
                              color = -log10(over_represented_pvalue))) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = "Top GO Terms Across Categories",
           subtitle = "Top 10 GO terms by p-value (regardless of significance)",
           x = NULL,
           y = NULL,
           size = "Gene Count",
           color = "-log10(p-value)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("GO_terms_top_dotplot_comparison.pdf", top_dotplot, width = 12, height = 8)
  }
}

# Create a list for venn diagram with proper empty list handling
go_term_list <- list(
  "Overlap GO Terms" = if(length(overlap_sig_terms) > 0) overlap_sig_terms else character(0),
  "PCAdapt-only GO Terms" = if(length(pcadapt_sig_terms) > 0) pcadapt_sig_terms else character(0),
  "pRDA-only GO Terms" = if(length(prda_sig_terms) > 0) prda_sig_terms else character(0)
)

# Create directory for Venn diagram output
dir.create("venn_temp", showWarnings = FALSE)

# Try the ggVennDiagram approach first (handles empty sets better)
library(ggVennDiagram)

# Create and save the ggVenn plot
p <- ggVennDiagram(go_term_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("GO Term Overlap Between Methods (FDR < 0.05)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )

ggsave("GO_terms_ggvenn_diagram.png", p, width = 8, height = 6, dpi = 300)

# Also try the traditional VennDiagram approach with error handling
tryCatch({
  venn.plot <- venn.diagram(
    x = go_term_list,
    category.names = c("Overlap\nGO Terms", "PCAdapt-only\nGO Terms", "pRDA-only\nGO Terms"),
    filename = "venn_temp/GO_terms_venn_diagram.png",
    output = TRUE,
    
    # Output features
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'solid',
    fill = c('#b2df8a', '#FF8000', '#1f78b4'),
    alpha = 0.5,
    
    # Numbers
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    cat.dist = 0.07,
    cat.pos = c(0, 0, 180)
  )
}, error = function(e) {
  cat("VennDiagram error:", e$message, "\n")
  cat("Using ggVennDiagram output instead.\n")
})

# Create a gene-level Venn diagram too
gene_list <- list(
  "Overlap Genes" = overlap_genes,
  "PCAdapt-only Genes" = pcadapt_only_genes,
  "pRDA-only Genes" = prda_only_genes
)

# Create gene Venn diagram
gene_venn <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#E69F00") +
  ggtitle("Gene Overlap Between Methods") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16)
  )

ggsave("Gene_venn_diagram.png", gene_venn, width = 8, height = 6, dpi = 300)

# Save workspace
save.image("GO_term_analysis.RData")

# Print summary information
cat("\n===== ANALYSIS SUMMARY =====\n")
cat("Total genes analyzed:\n")
print(gene_summary)
cat("\nSignificant GO terms:\n")
cat("- Overlap: ", length(overlap_sig_terms), "\n")
cat("- PCAdapt-only: ", length(pcadapt_sig_terms), "\n")
cat("- pRDA-only: ", length(prda_sig_terms), "\n")


# Create a more comprehensive dotplot showing all significant GO terms from all categories
# First, prepare the data
all_sig_terms <- bind_rows(
  overlap_go %>%
    filter(FDR < 0.05) %>%
    mutate(Category = "Overlap"),
  pcadapt_only_go %>%
    filter(FDR < 0.05) %>%
    mutate(Category = "PCAdapt-only"),
  prda_only_go %>%
    filter(FDR < 0.05) %>%
    mutate(Category = "pRDA-only")
)

# Create a copy of the terms ordered by p-value for y-axis ordering (most significant at top)
ordered_terms <- all_sig_terms %>%
  arrange(over_represented_pvalue) %>%  # Smallest p-values first = most significant
  pull(term) %>%
  unique() %>%
  rev()  # This reverses the order so the most significant appears at the top

# Make term a factor to preserve the p-value ordering across all categories
all_sig_terms$term <- factor(all_sig_terms$term, levels = ordered_terms)

# Set Category as a factor with the desired order
all_sig_terms$Category <- factor(all_sig_terms$Category, 
                                 levels = c("Overlap", "PCAdapt-only", "pRDA-only"))

write.csv(all_sig_terms, "all_sig_terms.csv")

# Create the plot with terms ordered by p-value
combined_plot <- ggplot(all_sig_terms,
                        aes(x = Category, y = term,
                            size = numDEInCat,
                            color = -log10(over_represented_pvalue))) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Significant GO Terms Across Methods",
       subtitle = paste("Overlap:", sum(all_sig_terms$Category == "Overlap"),
                        "terms | PCAdapt-only:", sum(all_sig_terms$Category == "PCAdapt-only"),
                        "terms | pRDA-only:", sum(all_sig_terms$Category == "pRDA-only"), "terms"),
       x = NULL,
       y = NULL,
       size = "Gene Count",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

combined_plot

ggsave("FDR005_pvalue_ordered_GO_terms_combined_dotplot.pdf", combined_plot, width = 12, height = max(8, nrow(all_sig_terms) * 0.3))
ggsave("FDR005_pvalue_ordered_GO_terms_combined_dotplot.png", combined_plot, width = 12, height = max(8, nrow(all_sig_terms) * 0.3), dpi =400)







############## PBE analysis ################
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/PBS")
# --- Setup ---
# Load FST results
fst_results <- read.table("combined_fst_results.txt", header = TRUE)

# Filter rows with NaN or FST >= 1 in any FST column
fst_cols <- c("FST_Rural_Urban", "FST_Rural_US", "FST_Urban_US")
valid_rows <- complete.cases(fst_results[, fst_cols]) &
  apply(fst_results[, fst_cols], 1, function(row) all(row < 1))
filtered_fst <- fst_results[valid_rows, ]
cat("SNPs after removing NaN and FST >= 1:", nrow(filtered_fst), "\n")

# Set negative FST values to 0
filtered_fst[, fst_cols] <- lapply(filtered_fst[, fst_cols], function(x) ifelse(x < 0, 0, x))

# --- Calculate Branch Lengths (T) ---
filtered_fst$T_Rural_Urban <- -log(1 - filtered_fst$FST_Rural_Urban)
filtered_fst$T_Rural_US    <- -log(1 - filtered_fst$FST_Rural_US)
filtered_fst$T_Urban_US    <- -log(1 - filtered_fst$FST_Urban_US)

# --- Calculate Standard PBS ---
filtered_fst$PBS_Urban <- (filtered_fst$T_Rural_Urban + filtered_fst$T_Urban_US - filtered_fst$T_Rural_US) / 2
filtered_fst$PBS_Rural <- (filtered_fst$T_Rural_Urban + filtered_fst$T_Rural_US - filtered_fst$T_Urban_US) / 2
filtered_fst$PBS_USA   <- (filtered_fst$T_Rural_US + filtered_fst$T_Urban_US - filtered_fst$T_Rural_Urban) / 2

# --- Calculate PBE (Population Branch Excess) ---
# Calculate required median values
med_PBS_Urban <- median(filtered_fst$PBS_Urban, na.rm = TRUE)
med_PBS_Rural <- median(filtered_fst$PBS_Rural, na.rm = TRUE)
med_T_Rural_US <- median(filtered_fst$T_Rural_US, na.rm = TRUE) # T_BC for PBE_Urban
med_T_Urban_US <- median(filtered_fst$T_Urban_US, na.rm = TRUE) # T_BC for PBE_Rural

# Check if median T values are zero to avoid division by zero
if (med_T_Rural_US <= 0) {
  warning("Median T_Rural_US is zero or negative. PBE_Urban cannot be reliably calculated.")
  filtered_fst$PBE_Urban <- NA
} else {
  filtered_fst$PBE_Urban <- filtered_fst$PBS_Urban - (filtered_fst$T_Rural_US * med_PBS_Urban / med_T_Rural_US)
}

if (med_T_Urban_US <= 0) {
  warning("Median T_Urban_US is zero or negative. PBE_Rural cannot be reliably calculated.")
  filtered_fst$PBE_Rural <- NA
} else {
  filtered_fst$PBE_Rural <- filtered_fst$PBS_Rural - (filtered_fst$T_Urban_US * med_PBS_Rural / med_T_Urban_US)
}

# --- Sort and Output Results ---
urban_pbe_sorted <- filtered_fst[order(-filtered_fst$PBE_Urban), ]
rural_pbe_sorted <- filtered_fst[order(-filtered_fst$PBE_Rural), ]

# Print the top results for each population based on PBE
cat("\nTop SNPs with highest PBE for Urban Shanghai:\n")
print(head(urban_pbe_sorted[, c("CHROM", "POS", "PBE_Urban", "PBS_Urban")]))

cat("\nTop SNPs with highest PBE for Rural Shanghai:\n")
print(head(rural_pbe_sorted[, c("CHROM", "POS", "PBE_Rural", "PBS_Rural")]))

# Write results to files
write.table(filtered_fst, "pbe_results_all.txt", sep="\t", quote=FALSE, row.names=FALSE)

# --- Load Libraries for Visualization ---
library(ggplot2)
library(dplyr)
library(ggrepel)

# Merge with gene annotation data
filtered_fst_merged <- merge(filtered_fst, annot_combined, by = "POS", all.x = TRUE)

# Function to filter gene names
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name) &
        Gene_Name != "" &
        !grepl("^LDel_", Gene_Name) & 
        !grepl("-", Gene_Name)
    )
}

# Count total SNPs in each population
urban_count <- nrow(urban_data_pbe)
rural_count <- nrow(rural_data_pbe)

# Create data frame for plotting
snp_count_data <- data.frame(
  Population = c("Urban", "Rural"),
  Count = c(urban_count, rural_count)
)


# 1. SNPs with positive PBE_Urban AND negative PBE_Rural (Urban-specific selection)
urban_specific <- filtered_fst_merged %>%
  filter(PBE_Urban > 0 & PBE_Rural < 0) %>%
  mutate(Selection_Pattern = "Urban-specific")

# 2. SNPs with negative PBE_Urban AND positive PBE_Rural (Rural-specific selection)
rural_specific <- filtered_fst_merged %>%
  filter(PBE_Urban < 0 & PBE_Rural > 0) %>%
  mutate(Selection_Pattern = "Rural-specific")

# 3. SNPs with positive values in BOTH PBE_Urban AND PBE_Rural (Divergent selection)
both_positive <- filtered_fst_merged %>%
  filter(PBE_Urban > 0 & PBE_Rural > 0) %>%
  mutate(Selection_Pattern = "Both_positive")

# 4. SNPs with negative values in BOTH PBE_Urban AND PBE_Rural (Conserved/Background)
both_negative <- filtered_fst_merged %>%
  filter(PBE_Urban < 0 & PBE_Rural < 0) %>%
  mutate(Selection_Pattern = "Both_negative")

# Create a summary of counts for each category
category_counts <- data.frame(
  Category = c("Urban-specific", "Rural-specific", "Both positive", "Both negative"),
  Count = c(nrow(urban_specific), nrow(rural_specific), 
            nrow(both_positive), nrow(both_negative)),
  Percentage = c(nrow(urban_specific), nrow(rural_specific), 
                 nrow(both_positive), nrow(both_negative)) / nrow(filtered_fst_merged) * 100
)

# Print the summary
print(category_counts)

# Optional: Combine all four datasets with their categories for potential unified analysis
all_categorized <- bind_rows(urban_specific, rural_specific, both_positive, both_negative)

library(dplyr)

# Step 1: Clean up column names and remove duplicated columns
filtered_fst_merged <- filtered_fst_merged %>%
  select(-matches("\\.y$")) %>%             # remove all ".y" columns
  rename_with(~ sub("\\.x$", "", .), ends_with(".x"))  # rename ".x" → no suffix


###### For plotting
# Urban-specific selection: positive PBE_Urban, negative PBE_Rural
urban_specific <- filtered_fst_merged %>%
  select(CHROM, POS, pos_cum, chr_color, Gene_Name, PBE_Urban, PBE_Rural) %>%
  filter(PBE_Urban > 0 & PBE_Rural < 0) %>%
  mutate(PBE_plot = PBE_Urban, Population = "Urban")

# Rural-specific selection: positive PBE_Rural, negative PBE_Urban
rural_specific <- filtered_fst_merged %>%
  select(CHROM, POS, pos_cum, chr_color, Gene_Name, PBE_Urban, PBE_Rural) %>%
  filter(PBE_Rural > 0 & PBE_Urban < 0) %>%
  mutate(PBE_plot = -PBE_Rural, Population = "Rural")

# Print counts to see distribution
cat("Urban-specific SNPs:", nrow(urban_specific), "\n")
cat("Rural-specific SNPs:", nrow(rural_specific), "\n")
cat("Total environment-specific SNPs:", nrow(combined_data_pbe), "\n")

# Create a unique SNP identifier for both datasets
urban_specific <- urban_specific %>%
  mutate(SNP = paste0(CHROM, "_", POS))

rural_specific <- rural_specific %>%
  mutate(SNP = paste0(CHROM, "_", POS))


# Update the combined data with the SNP identifier
combined_data_pbe <- bind_rows(urban_specific, rural_specific)

# Print the first few rows of the combined data
head(combined_data_pbe)

# Count total SNPs in each environment
urban_count <- nrow(urban_specific)
rural_count <- nrow(rural_specific)

mirror_data <- data.frame(
  Population = c("Urban", "Rural"),  # Shorter names for cleaner display
  Count = c(urban_count, -rural_count)  # Negative for rural to point downward
)

# Create the mirror barplot with closer bars
p_mirror_barplot <- ggplot(mirror_data, aes(x = Population, y = Count, fill = Population)) +
  # Make bars slimmer
  geom_bar(stat = "identity", width = 0.2) +
  scale_fill_manual(values = c("Urban" = "#3690c0", "Rural" = "#3690c0")) +
  
  # Add a zero line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  
  # Use absolute values for y-axis labels
  scale_y_continuous(labels = abs) +
  
  # Make the plot much narrower to bring bars closer together
  coord_fixed(ratio = 0.01) +  # This makes the plot very narrow
  
  # Clean theme
  labs(x = NULL,
       y = "Number of SNPs",
       title = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Save the narrow plot
ggsave("SLF_SNP_Mirror_Barplot_Narrow.pdf", p_mirror_barplot, 
       width = 3, height = 8, dpi = 300)


# Apply filtering
pbe_urban_valid_genes <- filter_valid_genes(urban_specific)
pbe_rural_valid_genes <- filter_valid_genes(rural_specific)

# Keep only the first occurrence for each gene
pbe_urban_dedup <- pbe_urban_valid_genes %>%
  group_by(Gene_Name) %>%
  slice(1) %>%
  ungroup()

pbe_rural_dedup <- pbe_rural_valid_genes %>%
  group_by(Gene_Name) %>%
  slice(1) %>%
  ungroup()

# --- Prepare PBE Data for Heatmap ---
combined_pbe_genes <- bind_rows(
  mutate(pbe_urban_dedup, Category = "Urban Specific (PBE)"),
  mutate(pbe_rural_dedup, Category = "Rural Specific (PBE)")
)

# Reshape for heatmap
heatmap_data_pbe <- combined_pbe_genes %>%
  select(Gene_Name, PBE_Urban, PBE_Rural, Category) %>%
  pivot_longer(
    cols = c(PBE_Urban, PBE_Rural),
    names_to = "Population_Metric",
    values_to = "PBE_Value"
  ) %>%
  mutate(
    Population = factor(Population_Metric,
                        levels = c("PBE_Urban", "PBE_Rural"),
                        labels = c("Shanghai Urban", "Shanghai Rural"))
  )





# --- Create PBE ggplot manhattan plot---
# Determine range for diverging color scale
pbe_range <- range(heatmap_data_pbe$PBE_Value, na.rm = TRUE)
pbe_midpoint <- 0
pbe_limit <- max(abs(pbe_range))

# --- Create Manhattan plot for PBE ---
# Extract chromosome numbers and prepare data
filtered_fst_merged$chr_num <- as.numeric(gsub("chr_", "", filtered_fst_merged$CHROM))
filtered_fst_merged <- filtered_fst_merged %>% arrange(chr_num, POS)

# Calculate chromosome boundaries
chr_pos <- filtered_fst_merged %>%
  group_by(CHROM) %>%
  summarize(
    chr_len = max(POS),
    chr_center = mean(POS)
  ) %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

# Create cumulative position
chr_cum_pos <- 0
filtered_fst_merged$pos_cum <- 0
for (i in 1:nrow(chr_pos)) {
  chr_snps <- filtered_fst_merged$CHROM.x == chr_pos$CHROM[i]
  filtered_fst_merged$pos_cum[chr_snps] <- filtered_fst_merged$POS[chr_snps] + chr_cum_pos
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5000000
}

# Threshold for highlighting
highlight_threshold_pbe <- 0

# STEP 1: Ensure pos_cum is correctly computed
combined_data_pbe <- combined_data_pbe %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(pos_cum = as.numeric(POS) + ifelse(row_number() == 1, 0, lag(cumsum(as.numeric(POS)), default = 0))) %>%
  ungroup()

# STEP 2: Define specific subsets again AFTER pos_cum
urban_specific <- combined_data_pbe %>%
  filter(PBE_Urban > highlight_threshold_pbe & PBE_Rural < 0) %>%
  mutate(PBE_plot = PBE_Urban, Population = "Urban")

rural_specific <- combined_data_pbe %>%
  filter(PBE_Rural > highlight_threshold_pbe & PBE_Urban < 0) %>%
  mutate(PBE_plot = -PBE_Rural, Population = "Rural")

# STEP 3: Highlight functional genes with pos_cum included
highlight_urban_pbe <- urban_specific %>%
  filter(!is.na(Gene_Name),
         Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  top_n(40, PBE_Urban)

highlight_rural_pbe <- rural_specific %>%
  filter(!is.na(Gene_Name),
         Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  top_n(40, PBE_Rural)

highlight_genes_pbe <- bind_rows(highlight_urban_pbe, highlight_rural_pbe)


# Plot: PBE Manhattan with only functional genes labeled
ggplot(combined_data_pbe, aes(x = pos_cum, y = PBE_plot)) +
  geom_point(color = "#3690c0", alpha = 0.6, size = 1) +
  
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_hline(yintercept = c(-highlight_threshold_pbe, highlight_threshold_pbe),
             linetype = "dashed", color = "darkred") +
  
  geom_point(data = highlight_genes_pbe,
             aes(x = pos_cum, y = PBE_plot),
             color = "black", fill = "#c861a5", shape = 21, size = 2, stroke = 0.4) +
  
  geom_text_repel(data = highlight_genes_pbe,
                  aes(x = pos_cum, y = PBE_plot, label = Gene_Name),
                  size = 3.2, fontface = "bold", color = "black",
                  box.padding = 0.6, point.padding = 0.5,
                  segment.size = 0.3, segment.color = "gray50",
                  max.overlaps = 15,
                  min.segment.length = 0,
                  force = 2) +
  
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  scale_y_continuous(labels = abs, limits = c(-2, 2)) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(x = NULL,
       y = "Population Branch Excess (PBE)",
       title = "Selection Signals in Spotted Lanternfly Populations")






# ─── Load Required Packages ───────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)

# ─── INPUT: Ensure `combined_data_pbe` contains CHROM, POS, PBE_Urban, PBE_Rural, Gene_Name ───

# ─── STEP 1: Compute cumulative genomic positions (pos_cum) ────────────────
combined_data_pbe$chr_num <- as.numeric(gsub("chr_", "", combined_data_pbe$CHROM))
combined_data_pbe <- combined_data_pbe %>% arrange(chr_num, POS)

# Calculate chromosome sizes and centers
chr_pos <- combined_data_pbe %>%
  group_by(CHROM) %>%
  summarize(
    chr_len = max(POS, na.rm = TRUE),
    chr_center = mean(POS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

# Add pos_cum to main dataset
chr_cum_pos <- 0
combined_data_pbe$pos_cum <- 0
for (i in seq_len(nrow(chr_pos))) {
  chr_name <- chr_pos$CHROM[i]
  idx <- combined_data_pbe$CHROM == chr_name
  combined_data_pbe$pos_cum[idx] <- combined_data_pbe$POS[idx] + chr_cum_pos
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5e6
}

# ─── STEP 2: Identify Urban- and Rural-specific SNPs ────────────────────────
highlight_threshold <- 0  # You can adjust to focus on stronger signals

urban_specific <- combined_data_pbe %>%
  filter(PBE_Urban > highlight_threshold, PBE_Rural < 0) %>%
  mutate(PBE_plot = PBE_Urban, Population = "Urban")

rural_specific <- combined_data_pbe %>%
  filter(PBE_Rural > highlight_threshold, PBE_Urban < 0) %>%
  mutate(PBE_plot = -PBE_Rural, Population = "Rural")

# Combine
combined_pbe_mirror <- bind_rows(urban_specific, rural_specific)

# ─── STEP 3: Filter and Highlight Top Functional Genes ──────────────────────
highlight_urban <- urban_specific %>%
  filter(!is.na(Gene_Name), Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  top_n(40, PBE_Urban)

highlight_rural <- rural_specific %>%
  filter(!is.na(Gene_Name), Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  top_n(40, PBE_Rural)

highlight_genes <- bind_rows(highlight_urban, highlight_rural)

# ─── STEP 4: Plot Mirror Manhattan Plot ─────────────────────────────────────
mirror_plot <- ggplot(combined_pbe_mirror, aes(x = pos_cum, y = PBE_plot)) +
  geom_point(color = "#3690c0", alpha = 0.6, size = 1) +
  
  # Highlight lines
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  geom_hline(yintercept = c(-highlight_threshold, highlight_threshold), 
             linetype = "dashed", color = "darkred", size = 0.3) +
  
  # Highlight top genes
  geom_point(data = highlight_genes, 
             aes(x = pos_cum, y = PBE_plot), 
             color = "black", fill = "#c861a5", shape = 21, size = 2, stroke = 0.4) +
  
  geom_text_repel(data = highlight_genes,
                  aes(label = Gene_Name),
                  size = 3.2, fontface = "bold", color = "black",
                  box.padding = 0.6, point.padding = 0.5,
                  segment.size = 0.3, segment.color = "gray50",
                  max.overlaps = 15, min.segment.length = 0, force = 2) +
  
  # Format chromosome axis
  scale_x_continuous(
    breaks = chr_pos$center_cum,
    labels = chr_pos$CHROM,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  scale_y_continuous(
    limits = c(-2, 2),  # Adjust based on your real PBE range
    breaks = seq(-2, 2, 0.5),
    labels = abs
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Population Branch Excess (PBE)",
    title = NULL
  )

# ─── STEP 5: Save Plot ──────────────────────────────────────────────────────
ggsave("Mirror_Manhattan_PBE_FunctionalGenes.png", mirror_plot, width = 12, height = 8, dpi = 300)
ggsave("Mirror_Manhattan_PBE_FunctionalGenes.pdf", mirror_plot, width = 12, height = 8)

# Create a cleaner Manhattan plot for PBE values
p_mirror_pbe <- ggplot(combined_data_pbe, aes(x = pos_cum, y = PBE_plot)) +
  # Use a single color for all points
  geom_point(color = "#3690c0", alpha = 0.6, size = 1) +
  
  # Keep essential reference lines
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_hline(yintercept = c(-highlight_threshold_pbe, highlight_threshold_pbe),
             linetype = "dashed", color = "darkred") +
  
  # Highlight significant genes with yellow points
  geom_point(data = highlight_genes_pbe,
             aes(x = pos_cum, y = PBE_plot),
             color = "black", fill = "#c861a5", shape = 21, size = 2, stroke = 0.4) +
  
  # Control label density and appearance
  geom_text_repel(data = highlight_genes_pbe,
                  aes(label = Gene_Name),
                  size = 3.2, fontface = "bold", color = "black",
                  box.padding = 0.6, point.padding = 0.5, 
                  segment.size = 0.3, segment.color = "gray50",
                  max.overlaps = 15,
                  min.segment.length = 0,  # Always draw lines regardless of distance
                  force = 2) +
  
  # Simplify axis formatting
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  theme_bw() +
  scale_y_continuous(labels = abs, limits = c(-2, 2)) +
  
  # Clean up theme elements
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  
  # Clear labels
  labs(x = NULL,
       y = "Population Branch Excess (PBE)",
       title = "Selection Signals in Spotted Lanternfly Populations")





print(p_mirror_pbe)

ggsave("SLF_PBE_Manhattan_Plot.pdf", p_mirror_pbe, width = 12, height = 8, dpi = 300)


# Combine plots using patchwork
combined_plot <- p_mirror_pbe  + p_mirror_barplot +
  plot_layout(widths = c(20, 3))

ggsave("SLF_Manhattan_with_Mirror_Barplot.pdf", combined_plot, 
       width = 14, height = 8, dpi = 300)




# Check the structure of both datasets
str(urban_specific)
str(freqs_wide)

# Create SNP column in urban_specific that matches the format of SNP_ID in freqs_merged_all
# Create a clean subset of columns from freqs_merged_all that we want to include
freqs_subset <- freqs_wide %>%
  select(SNP_ID, SH_urban, SH_rural, USA)

freqs_subset$SNP_ID <- gsub("_[ACGT]$", "", freqs_subset$SNP_ID) 

# STEP 1: Construct SNP_ID in urban_specific
urban_specific <- urban_specific %>%
  mutate(SNP_ID = paste0(CHROM, "_", POS))

# STEP 3: Left join using SNP_ID
urban_with_freqs <- left_join(urban_specific, freqs_subset, by = "SNP_ID")

# Extract annotation information
snp_annotations <- urban_with_freqs %>%
  select(SNP_ID, Gene_Name) %>%
  distinct()

# Create cleaner row labels
row_labels <- snp_annotations %>%
  mutate(
    # Use gene name for known genes, otherwise use shortened SNP_ID
    display_label = case_when(
      !is.na(Gene_Name) & !grepl("^LDel_", Gene_Name) & !grepl("-", Gene_Name) ~ Gene_Name,
      TRUE ~ ifelse(is.na(SNP_ID), 
                    paste0(CHROM, "_", POS), # Use CHROM_POS if SNP_ID is NA
                    str_replace(SNP_ID, "chr_", "")) # Otherwise clean the SNP_ID
    )
  ) %>%
  pull(display_label)
names(row_labels) <- snp_annotations$SNP_ID

# Create a vector for row label font face (bold for gene names)
is_gene <- !is.na(snp_annotations$Gene_Name) & 
  !grepl("^LDel_", snp_annotations$Gene_Name) & 
  !grepl("-", snp_annotations$Gene_Name)
row_fontface <- ifelse(is_gene, "bold", "plain")
names(row_fontface) <- snp_annotations$SNP_ID

### miss code for making urban matrix

# Define custom cell function to display values with optimal contrast
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Get the cell value
  cell_value <- urban_matrix[i, j]
  
  # Determine text color based on background darkness
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  
  # Add text with appropriate color and formatting
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}


# Create the heatmap
urban_heatmap <- pheatmap(
  urban_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,  # Don't cluster rows so they stay in SH_urban_freq order
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize_row = 10,
  cellwidth = 60,
  cellheight = 18,
  border_color = "white",
  main = "Allele Frequencies of Urban-Specific SNPs Under Selection\n(Sorted by Urban Frequency)",
  fontsize = 13,
  cell.fun = custom_cell_fun,
  labels_row = row_labels,
  fontface_row = row_fontface,  # Apply bold to gene names
  angle_col = 45
)

# Save to PDF
pdf("urban_specific_SNPs_sorted_by_urban_freq.pdf", width = 8, height = 40)
urban_heatmap
dev.off()





### What if we just focus on functional annotated, not showing with SNPID 

library(pheatmap)
library(grid)
library(dplyr)
library(stringr)


# Function to filter gene names - keeping only "valid" gene names
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name) &
        Gene_Name != "" &
        !grepl("^LDel_", Gene_Name) & 
        !grepl("-", Gene_Name)
    )
}

# Apply the filter to keep only records with valid gene names
urban_with_freqs_genes <- filter_valid_genes(urban_with_freqs)

# Sort by SH_urban_freq in descending order
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  arrange(desc(SH_urban_freq))

# Prepare the data matrix for the heatmap
urban_matrix <- as.matrix(urban_with_freqs_genes[, c("SH_rural_freq", "SH_urban_freq", "USA_freq")])
rownames(urban_matrix) <- urban_with_freqs_genes$SNP_ID
colnames(urban_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Use gene names as row labels
row_labels <- urban_with_freqs_genes$Gene_Name
names(row_labels) <- urban_with_freqs_genes$SNP_ID

# Define custom cell function to display values with optimal contrast
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Get the cell value
  cell_value <- urban_matrix[i, j]
  
  # Determine text color based on background darkness
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  
  # Add text with appropriate color and formatting
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Create the heatmap focused on functional genes
PBE_urban_heatmap <- pheatmap(
  urban_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,  # We'll add custom text
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "Genes under selection in SH Urban",
  fontsize = 10,
  cell.fun = custom_cell_fun,  # Custom function for cell rendering,
  labels_row = row_labels,
  angle_col = 45
)

# Save to PDF
pdf("urban_specific_functional_genes_heatmap.pdf", width = 8, height = 10)
PBE_urban_heatmap
dev.off()



library(pheatmap)
library(grid)
library(dplyr)
library(stringr)

# Apply the filter to keep only records with valid gene names
rural_with_freqs_genes <- filter_valid_genes(rural_with_freqs)

# Sort by SH_rural_freq in descending order (since this is rural data)
rural_with_freqs_genes <- rural_with_freqs_genes %>%
  arrange(desc(SH_rural_freq))

# Prepare the data matrix for the heatmap
rural_matrix <- as.matrix(rural_with_freqs_genes[, c("SH_rural_freq", "SH_urban_freq", "USA_freq")])
rownames(rural_matrix) <- rural_with_freqs_genes$SNP_ID
colnames(rural_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Use gene names as row labels
row_labels <- rural_with_freqs_genes$Gene_Name
names(row_labels) <- rural_with_freqs_genes$SNP_ID

# Define custom cell function to display values with optimal contrast
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  # Get the cell value
  cell_value <- rural_matrix[i, j]
  
  # Determine text color based on background darkness
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  
  # Add text with appropriate color and formatting
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Create the heatmap focused on functional genes
PBE_rural_heatmap <- pheatmap(
  rural_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,  # We'll add custom text
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "Genes under selection in SH Rural",
  fontsize = 10,
  cell.fun = custom_cell_fun,  # Custom function for cell rendering
  labels_row = row_labels,
  angle_col = 45
)

# Save to PDF
pdf("rural_specific_functional_genes_heatmap.pdf", width = 8, height = 10)
PBE_rural_heatmap
dev.off()




# Function to filter valid genes
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name) &
        Gene_Name != "" &
        !grepl("^LDel_", Gene_Name) & 
        !grepl("-", Gene_Name)
    )
}

# Apply filter to each dataframe
fdr_merged_filtered <- filter_valid_genes(fdr_merged)
fdr_USA_Only_filtered <- filter_valid_genes(fdr_USA_Only)
overlap_full_info_subset_filtered <- filter_valid_genes(freqs_merged)

write.table(overlap_full_info_subset_filtered$Gene_Name, file = "overlap_full_info_subset_filtered_genes.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)

intersect(fdr_merged$SNP, freqs_merged$SNP_ID)






library(dplyr)
library(pheatmap)
library(grid)

# Step 1: Fix SNP_ID format in freqs_wide (remove allele info like "_A")
freqs_subset <- freqs_wide %>%
  mutate(SNP_ID_clean = gsub("_[ACGT]$", "", SNP_ID)) %>%
  select(SNP_ID_clean, SH_urban, SH_rural, USA, SNP_ID)

# Step 2: Create SNP_ID to match in urban_specific
urban_specific <- urban_specific %>%
  mutate(SNP_ID_clean = paste0(CHROM, "_", POS))

# Step 3: Join with allele frequency data
urban_with_freqs <- left_join(urban_specific, freqs_subset, by = "SNP_ID_clean")

# Step 4: Filter to only SNPs with valid gene names
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name),
      Gene_Name != "",
      !grepl("^LDel_", Gene_Name),
      !grepl("-", Gene_Name)
    )
}

urban_with_freqs_genes <- filter_valid_genes(urban_with_freqs)

# Step 5: Sort by SH_urban frequency
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  arrange(desc(SH_urban))

# Step 6: Create better row labels: GeneName_SNPID
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  mutate(Row_Label = paste0(Gene_Name, "_", CHROM, "_", POS))

# Step 7: Create matrix and label rows
urban_matrix <- as.matrix(urban_with_freqs_genes[, c("SH_rural", "SH_urban", "USA")])
rownames(urban_matrix) <- urban_with_freqs_genes$Row_Label
colnames(urban_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Optional: custom cell text display (same as before)
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  cell_value <- urban_matrix[i, j]
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Step 8: Generate heatmap
PBE_urban_heatmap <- pheatmap(
  urban_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "SNPs under selection in SH Urban",
  fontsize = 10,
  cell.fun = custom_cell_fun,
  angle_col = 45
)

# Save to PDF
pdf("PBE_urban_specific_SNPs_heatmap.pdf", width = 8, height = 6)
PBE_urban_heatmap
dev.off()


library(dplyr)
library(pheatmap)
library(grid)

# Step 1: Fix SNP_ID format in freqs_wide (remove allele info like "_A")
freqs_subset <- freqs_wide %>%
  mutate(SNP_ID_clean = gsub("_[ACGT]$", "", SNP_ID)) %>%
  select(SNP_ID_clean, SH_urban, SH_rural, USA, SNP_ID)

# Step 2: Create SNP_ID to match in rural_specific
rural_specific <- rural_specific %>%
  mutate(SNP_ID_clean = paste0(CHROM, "_", POS))

# Step 3: Join with allele frequency data
rural_with_freqs <- left_join(rural_specific, freqs_subset, by = "SNP_ID_clean")

# Step 4: Filter to only SNPs with valid gene names
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name),
      Gene_Name != "",
      !grepl("^LDel_", Gene_Name),
      !grepl("-", Gene_Name)
    )
}

rural_with_freqs_genes <- filter_valid_genes(rural_with_freqs)

# Step 5: Sort by SH_rural frequency
rural_with_freqs_genes <- rural_with_freqs_genes %>%
  arrange(desc(SH_rural))

# Step 6: Create better row labels: GeneName_CHROM_POS
rural_with_freqs_genes <- rural_with_freqs_genes %>%
  mutate(Row_Label = paste0(Gene_Name, "_", CHROM, "_", POS))

# Step 7: Create matrix and label rows
rural_matrix <- as.matrix(rural_with_freqs_genes[, c("SH_rural", "SH_urban", "USA")])
rownames(rural_matrix) <- rural_with_freqs_genes$Row_Label
colnames(rural_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Optional: custom cell text display
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  cell_value <- rural_matrix[i, j]
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  grid.text(sprintf("%.2f", cell_value), x, y, 
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Step 8: Generate heatmap
PBE_rural_heatmap <- pheatmap(
  rural_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "SNPs under selection in SH Rural",
  fontsize = 10,
  cell.fun = custom_cell_fun,
  angle_col = 45
)

# Save to PDF
pdf("PBE_specific_functional_SNPs_heatmap.pdf", width = 8, height = 6)
PBE_rural_heatmap
dev.off()


# Count how many times each gene appears
gene_counts <- table(urban_with_freqs_genes$Gene_Name)

# Create updated gene labels based on count
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  mutate(Gene_Label = ifelse(gene_counts[Gene_Name] > 1,
                             paste0(Gene_Name, "_", SNP_ID),
                             Gene_Name))

# Count how many times each gene appears
gene_counts <- table(rural_with_freqs_genes$Gene_Name)



# Load required libraries
library(pheatmap)
library(grid)
library(dplyr)
library(stringr)

# ──────────── FUNCTION: Filter Valid Genes ────────────
filter_valid_genes <- function(df) {
  df %>%
    filter(
      !is.na(Gene_Name),
      Gene_Name != "",
      !grepl("^LDel_", Gene_Name),
      !grepl("-", Gene_Name)
    )
}

# ──────────── URBAN HEATMAP ────────────

# Filter valid genes
urban_with_freqs_genes <- filter_valid_genes(urban_with_freqs)

# Count gene occurrences
urban_gene_counts <- table(urban_with_freqs_genes$Gene_Name)

# Add label: Gene_Name_SNP_ID if gene occurs more than once
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  mutate(Gene_Label = ifelse(urban_gene_counts[Gene_Name] > 1,
                             paste0(Gene_Name, "_", SNP_ID),
                             Gene_Name))
urban_with_freqs_genes <- urban_with_freqs_genes %>%
  arrange(desc(SH_urban))  # <- use the correct column name


# Create matrix for heatmap
urban_matrix <- as.matrix(urban_with_freqs_genes[, c("SH_rural", "SH_urban", "USA")])
rownames(urban_matrix) <- urban_with_freqs_genes$SNP_ID
colnames(urban_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Assign row labels
urban_row_labels <- urban_with_freqs_genes$Gene_Label
names(urban_row_labels) <- urban_with_freqs_genes$SNP_ID

# Cell text contrast function
custom_cell_fun <- function(j, i, x, y, width, height, fill) {
  cell_value <- urban_matrix[i, j]
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  grid.text(sprintf("%.2f", cell_value), x, y,
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Create heatmap
PBE_urban_heatmap <- pheatmap(
  urban_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "Genes under selection in SH Urban",
  fontsize = 10,
  cell.fun = custom_cell_fun,
  labels_row = urban_row_labels,
  angle_col = 45
)

# Save plots
pdf("urban_specific_functional_genes_labelled.pdf", width = 6, height = 8)
print(PBE_urban_heatmap)
dev.off()

png("urban_specific_functional_genes_labelled.png", width = 600, height = 800, res = 300)
print(PBE_urban_heatmap)
dev.off()


# ──────────── RURAL HEATMAP ────────────

# Filter valid genes
rural_with_freqs_genes <- filter_valid_genes(rural_with_freqs)

# Sort by SH_rural_freq
rural_with_freqs_genes <- rural_with_freqs_genes %>%
  arrange(desc(SH_rural))

# Use gene name directly (no duplicate genes)
rural_row_labels <- rural_with_freqs_genes$Gene_Name
names(rural_row_labels) <- rural_with_freqs_genes$SNP_ID

# Create matrix for heatmap
rural_matrix <- as.matrix(rural_with_freqs_genes[, c("SH_rural", "SH_urban", "USA")])
rownames(rural_matrix) <- rural_with_freqs_genes$SNP_ID
colnames(rural_matrix) <- c("SH_Rural", "SH_Urban", "USA")

# Same cell function
custom_cell_fun_rural <- function(j, i, x, y, width, height, fill) {
  cell_value <- rural_matrix[i, j]
  text_col <- ifelse(cell_value > 0.5, "white", "black")
  grid.text(sprintf("%.2f", cell_value), x, y,
            gp = gpar(col = text_col, fontsize = 9, fontface = "bold"))
}

# Create heatmap
PBE_rural_heatmap <- pheatmap(
  rural_matrix,
  color = colorRampPalette(c("#c6dbef", "#4292c6", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize_row = 9,
  cellwidth = 40,
  cellheight = 12,
  border_color = NA,
  main = "Genes under selection in SH Rural",
  fontsize = 10,
  cell.fun = custom_cell_fun_rural,
  labels_row = rural_row_labels,
  angle_col = 45
)

# Save plots
pdf("rural_specific_functional_genes_labelled.pdf", width = 6, height = 8)
print(PBE_rural_heatmap)
dev.off()

png("rural_specific_functional_genes_labelled.png", width = 600, height = 800, res = 300)
print(PBE_rural_heatmap)
dev.off()


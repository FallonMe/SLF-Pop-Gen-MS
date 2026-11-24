
# Function to run pcadapt with different K values and find overlapping candidate loci
run_pcadapt_multiple_K <- function(bed_file, pop_data, K_values = c(2, 3, 4, 5, 6), annot_data = NULL, 
                                   fdr_threshold = 0.01, output_dir = "./USA_new/") {
  
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
    # Create Venn diagram for overlapping candidates (only show K2-K5 for visibility)
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
    
    # Find SNPs that overlap in all K values (ultra-conservative)
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

# Load required packages
library(pcadapt)
library(ggplot2)
library(qvalue)
library(dplyr)
library(tidyr)
library(VennDiagram)

# Run analysis with K2-K6
results <- run_pcadapt_multiple_K(
  bed_file = "./USA_final_annotated.bed",
  pop_data = "~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/USA_SLF_MS_Info.csv",
  K_values = c(2, 3, 4, 5, 6),  # Added K=6
  annot_data = annot_combined,
  fdr_threshold = 0.01,
  output_dir = "./USA_new/"  # Changed output directory
)

# Compare number of candidates at different K values
k2_candidates <- read.csv("./USA_new/candidate_loci_K2.csv")
k3_candidates <- read.csv("./USA_new/candidate_loci_K3.csv")
k4_candidates <- read.csv("./USA_new/candidate_loci_K4.csv")
k5_candidates <- read.csv("./USA_new/candidate_loci_K5.csv")
k6_candidates <- read.csv("./USA_new/candidate_loci_K6.csv")  # Added K=6

# Number of candidate loci at each K value
cat("Number of candidates:\n")
cat("K2:", nrow(k2_candidates), "\n")
cat("K3:", nrow(k3_candidates), "\n")
cat("K4:", nrow(k4_candidates), "\n")
cat("K5:", nrow(k5_candidates), "\n")
cat("K6:", nrow(k6_candidates), "\n")

# Load the overlapping SNPs across all K values (now K2-K6)
overlap_all_k <- read.csv("./USA_new/overlapping_candidate_loci_all_K.csv")
cat("Ultra-conservative consensus (K2-K6):", nrow(overlap_all_k), "SNPs\n")

# Create final results using K=6 p-values (as recommended)
k6_data <- k6_candidates
overlap_with_k6_pvals <- merge(overlap_all_k, k6_data[, c("SNP", "p_value", "q_value")], 
                               by = "SNP", suffixes = c("", "_K6"))

# Add annotation data
final_results <- merge(overlap_with_k6_pvals, annot_combined, by = "SNP", all.x = TRUE)
final_results <- final_results[order(final_results$p_value), ]

# Save final ultra-conservative results with K=6 p-values
write.csv(final_results, "./USA_new/Final_overlap_K2-K6_final.csv", row.names = FALSE)


# Load required libraries
library(VennDiagram)
library(grid)
library(gridExtra)

# Create candidate lists from the saved CSV files
k2_candidates <- read.csv("./USA_new/candidate_loci_K2.csv")
k3_candidates <- read.csv("./USA_new/candidate_loci_K3.csv")
k4_candidates <- read.csv("./USA_new/candidate_loci_K4.csv")
k5_candidates <- read.csv("./USA_new/candidate_loci_K5.csv")
k6_candidates <- read.csv("./USA_new/candidate_loci_K6.csv")

# Create the candidate_lists object
candidate_lists <- list(
  K2 = k2_candidates$SNP,
  K3 = k3_candidates$SNP,
  K4 = k4_candidates$SNP,
  K5 = k5_candidates$SNP,
  K6 = k6_candidates$SNP
)

# Create the plot object but don't save it yet
venn <- venn.diagram(
  x = candidate_lists,
  category.names = paste0("K=", c(2, 3, 4, 5, 6)),
  filename = NULL,
  output = TRUE,
  scaled = TRUE,
  col = "black",
  fill = c("red", "blue", "green", "purple", "orange"),  # Added orange for K6
  alpha = 0.5,
  main = "Overlapping candidate loci K2-K6"
)

# Open a PDF device
pdf("./USA_new/custom_venn_diagram_K2-K6.pdf", width = 10, height = 10)

# Draw the diagram
grid.newpage()
grid.draw(venn)

# Close the PDF device
dev.off()

# Also create a PNG version
png("./USA_new/custom_venn_diagram_K2-K6.png", width = 800, height = 800)
grid.newpage()
grid.draw(venn)
dev.off()

# Print summary of overlaps
cat("Summary of candidate SNPs:\n")
cat("K2:", length(candidate_lists$K2), "\n")
cat("K3:", length(candidate_lists$K3), "\n")
cat("K4:", length(candidate_lists$K4), "\n")
cat("K5:", length(candidate_lists$K5), "\n")
cat("K6:", length(candidate_lists$K6), "\n")


#######Go term analysis for USA 7318 SNPs
#########Mahattan plot
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/USA_new")
# Load data
snp_data_pcadpt <- read.csv("./Final_overlap_K2-K6_final.csv")
snp_data_prda <- read.csv("./prda_SLF_annotation.csv")

# Ensure both have the SNP_ID column
pcadapt_overlap <- snp_data_pcadpt %>%
  mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>%
  filter(SNP_ID %in% shared_ids)

prda_overlap <- snp_data_prda %>%
  mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>%
  filter(SNP_ID %in% shared_ids)

# Identify columns in prda_overlap that are NOT duplicated
cols_to_drop <- intersect(colnames(pcadapt_overlap), colnames(prda_overlap))
prda_overlap_clean <- prda_overlap %>%
  select(-all_of(cols_to_drop[cols_to_drop != "SNP_ID"]))  # keep SNP_ID for join

# Now merge without duplicate columns and without .x/.y suffixes
merged_overlap <- left_join(pcadapt_overlap, prda_overlap_clean, by = "SNP_ID")

# Save
write.csv(merged_overlap, "Overlap_SNPs_pcadapt_prda.csv", row.names = FALSE)
cat("✅ Merged without duplicated columns. Rows:", nrow(merged_overlap), "\n")


####### Manhattan plot with highligting the overlap SNPs with annotated names
#-------------------------------------------------------------------------------
# MANHATTAN PLOT WITH OVERLAP HIGHLIGHTING
#-------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load data
snp_data <- read.csv("./Final_overlap_K2-K6_final.csv")

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
write.csv(overlap_genes, "overlap_genes184_unqiuegenenames.csv")

# Use all overlap genes without filtering by significance
significant_snps <- overlap_genes

# Count how many overlap SNPs were found
overlap_count <- sum(snp_data$is_overlap, na.rm = TRUE)
print(paste("Found", overlap_count, "overlap SNPs out of expected 185"))

#-------------------------------------------------------------------------------
#  MANHATTAN PLOT for overlap VARIANTS
#-------------------------------------------------------------------------------
# Identify unique functional genes (not LDel, not empty, not containing '-')
highlight_snps_usa <- significant_snps %>%
  filter(!is.na(Gene_Name),
         Gene_Name != "",
         !grepl("^LDel_", Gene_Name),
         !grepl("-", Gene_Name)) %>%
  group_by(Gene_Name) %>%
  slice_min(order_by = FDR_p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Compute max log10(FDR) for scaling
max_log10p <- ceiling(max(-log10(snp_data$FDR_p_value), na.rm = TRUE))
# Find minimum FDR value
min_log10p <- floor(min(-log10(as.numeric(snp_data$FDR_p_value)), na.rm = TRUE))


# Manhattan plot
p <- ggplot(snp_data, aes(x = pos_cum, y = -log10(FDR_p_value), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#3690c0")) +
  scale_x_continuous(
    label = chr_pos$CHROM,
    breaks = chr_pos$center_cum,
    guide = guide_axis(check.overlap = TRUE)
  ) +
  scale_y_continuous(
    name = expression(-log[10](FDR[p-value])),
    limits = c(min_log10p - 1, max_log10p + 2),
    breaks = pretty(c(min_log10p -1 , max_log10p + 2)),
    expand = c(0, 0)
  ) +
  geom_point(
    data = snp_data %>% filter(is_overlap == TRUE),
    aes(x = pos_cum, y = -log10(FDR_p_value)),
    color = "black",
    fill = "#c861a5",
    size = 2,
    shape = 21,
    stroke = 0.4,
    alpha = 1
  ) +
  geom_text_repel(
    data = highlight_snps_usa,
    aes(x = pos_cum, y = -log10(FDR_p_value), label = Gene_Name),
    size = 3.2,
    fontface = "bold",
    color = "black",
    box.padding = 0.6,
    point.padding = 0.5,
    max.overlaps = 100,
    segment.size = 0.3,
    segment.color = "grey50"
  ) +
  labs(
    x = "Chromosome",
    title = "Manhattan Plot of SLF Selection Signals",
    subtitle = paste0("All points represent SNPs significant at FDR < 0.01; pink highlights show ", 
                      sum(snp_data$is_overlap), " SNPs identified by both pRDA and PCAdapt")
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
  )

# Print plot
print(p)

# Save the plot with alternating chromosome colors
ggsave("SLF_USA_Manhattan_Plot_Overlap.pdf", p, width = 14, height = 8, dpi = 300)
ggsave("SLF_USA_Manhattan_Plot_Overlap.png", p, width = 14, height = 8, dpi = 300)


library(ggplot2)
library(dplyr)
library(ggrepel)

# Load pcadapt SNPs
snp_data <- read.csv("./USA_new/Final_overlap_K2-K6_final.csv")

# Load overlap table with full annotations
merged_overlap <- read.csv("Overlap_SNPs_pcadapt_prda.csv")

# Generate consistent SNP_ID for matching
snp_data <- snp_data %>%
  mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = "_"))

# Extract chromosome number for plotting
snp_data$chr_num <- as.numeric(gsub("chr_", "", snp_data$CHROM))
snp_data <- snp_data %>% arrange(chr_num, POS)

# Alternating chromosome colors
snp_data$chr_color <- factor(snp_data$chr_num %% 2)

# Chromosome position and center
chr_pos <- snp_data %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS), chr_center = mean(POS), .groups = "drop") %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

# Add cumulative positions
chr_pos$center_cum <- 0
chr_cum_pos <- 0
snp_data$pos_cum <- 0

for (i in 1:nrow(chr_pos)) {
  chr_snps <- snp_data$CHROM == chr_pos$CHROM[i]
  snp_data$pos_cum[chr_snps] <- snp_data$POS[chr_snps] + chr_cum_pos
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5e6  # 5Mb gap
}

# Tag overlapping SNPs
snp_data$is_overlap <- snp_data$SNP_ID %in% merged_overlap$SNP_ID

# Pull highlightable SNPs with readable gene names
highlight_snps <- merged_overlap %>%
  filter(
    !is.na(Gene_Name),
    Gene_Name != "",
    !grepl("^LDel_", Gene_Name),
    !grepl("-", Gene_Name)
  ) %>%
  group_by(Gene_Name) %>%
  slice_min(order_by = FDR_p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Add position info for labeling
highlight_snps <- highlight_snps %>%
  mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = "_")) %>%
  left_join(snp_data[, c("SNP_ID", "pos_cum")], by = "SNP_ID")

# Build Manhattan plot
p <- ggplot(snp_data, aes(x = pos_cum, y = -log10(FDR_p_value), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#3690c0")) +
  scale_x_continuous(breaks = chr_pos$center_cum, labels = chr_pos$CHROM) +
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
    y = "-log10(FDR p-value)",
    title = "Manhattan Plot of SLF Selection Signals",
    subtitle = paste0("Pink points (", sum(snp_data$is_overlap), ") shared by pRDA and pcadapt")
  )

# Add highlight points (overlap)
p <- p + geom_point(
  data = snp_data %>% filter(is_overlap),
  aes(x = pos_cum, y = -log10(FDR_p_value)),
  color = "black", fill = "#c861a5", shape = 21,
  size = 2, stroke = 0.4, alpha = 1
)

# Add gene labels
p <- p + geom_text_repel(
  data = highlight_snps,
  aes(x = pos_cum, y = -log10(FDR_p_value), label = Gene_Name),
  size = 3.2, fontface = "bold", color = "black",
  box.padding = 0.6, point.padding = 0.5, max.overlaps = 40,
  segment.size = 0.3, segment.color = "grey50"
)

# Print and save
print(p)
ggsave("SLF_USA_Manhattan_Plot_Overlap.pdf", p, width = 14, height = 8, dpi = 300)
ggsave("SLF_USA_Manhattan_Plot_Overlap.png", p, width = 14, height = 8, dpi = 300)


#-------------------------------------------------------------------------------
# Go term analysis fro USA and overlap separetely
#-------------------------------------------------------------------------------
library(goseq)
library(dplyr)
library(ggplot2)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")
# Load your datasets
final_results <- read.csv("./USA_new/Final_overlap_K2-K6_final.csv")
merged_overlap <- read.csv("./USA_new/Overlap_SNPs_pcadapt_prda.csv")

# Extract significant gene list from PCAdapt results (your sig_gene_list equivalent)
sig_gene_list_pcadapt <- final_results %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

cat("PCAdapt significant genes:", length(sig_gene_list_pcadapt), "\n")

# Create all_gene_ids from your GO mapping (background genes)
all_gene_ids_pcadapt <- unique(gene_id_go_mapping$GeneName)
all_gene_ids_pcadapt <- all_gene_ids_pcadapt[!is.na(all_gene_ids_pcadapt) & all_gene_ids_pcadapt != ""]

cat("Background genes from GO mapping:", length(all_gene_ids_pcadapt), "\n")
cat("Significant genes in background:", sum(sig_gene_list_pcadapt %in% all_gene_ids_pcadapt), "\n")

# Create a binary vector where 1 = significant gene, 0 = not significant
gene_vector_pcadapt <- numeric(length(all_gene_ids_pcadapt))
names(gene_vector_pcadapt) <- all_gene_ids_pcadapt

# Populate the gene vector (EXACT REPLICATION)
for (i in 1:length(all_gene_ids_pcadapt)) {
  gene_id <- all_gene_ids_pcadapt[i]
  # Check if this gene is in our significant gene list
  gene_vector_pcadapt[i] <- ifelse(gene_id %in% sig_gene_list_pcadapt, 1, 0)
}

cat("Genes marked as significant in vector:", sum(gene_vector_pcadapt), "\n")

# Create gene2go dataframe for goseq (EXACT REPLICATION)
gene2go_pcadapt <- data.frame(
  gene = gene_id_go_mapping$GeneName,
  category = gene_id_go_mapping$GOterms
) %>%
  dplyr::distinct()

cat("Gene2GO mapping rows:", nrow(gene2go_pcadapt), "\n")

# Create bias values with sufficient uniqueness for goseq (EXACT REPLICATION)
# Using rank of gene IDs plus small random noise for bias
set.seed(123) # For reproducibility
bias_pcadapt <- rank(as.numeric(factor(names(gene_vector_pcadapt)))) + 
  runif(length(gene_vector_pcadapt), min=0, max=0.1)
names(bias_pcadapt) <- names(gene_vector_pcadapt)

# Create probability weighting function with permissive parameters (EXACT REPLICATION)
pwf_pcadapt <- nullp(gene_vector_pcadapt, bias.data = bias_pcadapt, plot.fit = FALSE)

# Run GO enrichment analysis with higher p-value threshold (EXACT REPLICATION)
GO_results_pcadapt <- goseq(pwf_pcadapt, gene2cat = gene2go_pcadapt, use_genes_without_cat = TRUE)

# Add FDR correction and process results (EXACT REPLICATION)
GO_results_all_pcadapt <- GO_results_pcadapt %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))

# Save all results
write.csv(GO_results_all_pcadapt, "PCAdapt_USA_GO_enrichment_all_results.csv", row.names = FALSE)

# Filter at a higher threshold to capture more results (EXACT REPLICATION)
GO_results_filtered_pcadapt <- GO_results_all_pcadapt %>%
  filter(FDR < 0.05) %>%  # More lenient threshold
  arrange(over_represented_pvalue)

cat("PCAdapt: Number of terms with FDR < 0.05:", nrow(GO_results_filtered_pcadapt), "\n")

if (nrow(GO_results_filtered_pcadapt) > 0) {
  # Generate visualization of GO enrichment results (EXACT REPLICATION)
  # Select the top significant GO terms
  top_GO_pcadapt <- GO_results_filtered_pcadapt %>%
    arrange(over_represented_pvalue) %>%
    head(30)
  
  # Create bar plot of enriched GO terms (EXACT REPLICATION)
  enrichment_plot_pcadapt <- ggplot(top_GO_pcadapt, 
                                    aes(x = reorder(term, -over_represented_pvalue), 
                                        y = -log10(over_represented_pvalue), fill = ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
    labs(title = "PCAdapt Over-represented GO Terms",
         x = NULL,
         y = "-log10(p-value)",
         fill = "Ontology") +
    theme_minimal()
  
  # Save the GO term enrichment plot
  ggsave("PCAdapt_USA_Selection_GO_top_terms.pdf", enrichment_plot_pcadapt, width = 12, height = 8)
  
  # Create a bubble plot with enrichment information (EXACT REPLICATION)
  top_GO_pcadapt <- top_GO_pcadapt %>%
    mutate(enrichment_ratio = numDEInCat / 
             (numInCat * sum(gene_vector_pcadapt) / length(gene_vector_pcadapt)))
  
  bubble_plot_pcadapt <- ggplot(top_GO_pcadapt, 
                                aes(x = enrichment_ratio, 
                                    y = reorder(term, enrichment_ratio),
                                    size = numDEInCat,
                                    color = -log10(over_represented_pvalue))) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(title = "PCAdapt GO Term Enrichment in Selection Candidates",
         x = "Enrichment Ratio (>1 = enriched)",
         y = "GO Term",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_minimal()
  
  # Save the bubble plot
  ggsave("PCAdapt_USA_Selection_GO_enrichment_bubble_plot.pdf", bubble_plot_pcadapt, width = 12, height = 8)
}

#===============================================================================
# STEP 2: EXACT SAME APPROACH FOR OVERLAP ANALYSIS
#===============================================================================

# Extract significant gene list from overlap results
sig_gene_list_overlap <- merged_overlap %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

cat("Overlap significant genes:", length(sig_gene_list_overlap), "\n")

# Use same background genes
all_gene_ids_overlap <- all_gene_ids_pcadapt  # Same background

cat("Overlap significant genes in background:", sum(sig_gene_list_overlap %in% all_gene_ids_overlap), "\n")

# Create a binary vector where 1 = significant gene, 0 = not significant
gene_vector_overlap <- numeric(length(all_gene_ids_overlap))
names(gene_vector_overlap) <- all_gene_ids_overlap

# Populate the gene vector (EXACT REPLICATION)
for (i in 1:length(all_gene_ids_overlap)) {
  gene_id <- all_gene_ids_overlap[i]
  # Check if this gene is in our significant gene list
  gene_vector_overlap[i] <- ifelse(gene_id %in% sig_gene_list_overlap, 1, 0)
}

cat("Overlap genes marked as significant in vector:", sum(gene_vector_overlap), "\n")

# Create gene2go dataframe for goseq (SAME AS BEFORE)
gene2go_overlap <- gene2go_pcadapt  # Same mapping

# Create bias values with sufficient uniqueness for goseq (EXACT REPLICATION)
set.seed(456) # Different seed for independence
bias_overlap <- rank(as.numeric(factor(names(gene_vector_overlap)))) + 
  runif(length(gene_vector_overlap), min=0, max=0.1)
names(bias_overlap) <- names(gene_vector_overlap)

# Create probability weighting function with permissive parameters (EXACT REPLICATION)
pwf_overlap <- nullp(gene_vector_overlap, bias.data = bias_overlap, plot.fit = FALSE)

# Run GO enrichment analysis (EXACT REPLICATION)
GO_results_overlap <- goseq(pwf_overlap, gene2cat = gene2go_overlap, use_genes_without_cat = TRUE)

# Add FDR correction and process results (EXACT REPLICATION)
GO_results_all_overlap <- GO_results_overlap %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))

# Save all results
write.csv(GO_results_all_overlap, "Overlap_Selection_GO_enrichment_all_results.csv", row.names = FALSE)

# Filter at a higher threshold to capture more results (EXACT REPLICATION)
GO_results_filtered_overlap <- GO_results_all_overlap %>%
  filter(FDR < 0.05) %>%  # More lenient threshold
  arrange(over_represented_pvalue)

cat("Overlap: Number of terms with FDR < 0.05:", nrow(GO_results_filtered_overlap), "\n")

if (nrow(GO_results_filtered_overlap) > 0) {
  # Generate visualization of GO enrichment results (EXACT REPLICATION)
  top_GO_overlap <- GO_results_filtered_overlap %>%
    arrange(over_represented_pvalue) %>%
    head(31)
  
  # Create bar plot of enriched GO terms (EXACT REPLICATION)
  enrichment_plot_overlap <- ggplot(top_GO_overlap, 
                                    aes(x = reorder(term, -over_represented_pvalue), 
                                        y = -log10(over_represented_pvalue), fill = ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
    labs(title = "Overlap Over-represented GO Terms",
         x = NULL,
         y = "-log10(p-value)",
         fill = "Ontology") +
    theme_minimal()
  
  # Save the GO term enrichment plot
  ggsave("Overlap_Selection_GO_top_terms.pdf", enrichment_plot_overlap, width = 12, height = 8)
  
  # Create a bubble plot with enrichment information (EXACT REPLICATION)
  top_GO_overlap <- top_GO_overlap %>%
    mutate(enrichment_ratio = numDEInCat / 
             (numInCat * sum(gene_vector_overlap) / length(gene_vector_overlap)))
  
  bubble_plot_overlap <- ggplot(top_GO_overlap, 
                                aes(x = enrichment_ratio, 
                                    y = reorder(term, enrichment_ratio),
                                    size = numDEInCat,
                                    color = -log10(over_represented_pvalue))) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(title = "Overlap GO Term Enrichment in Selection Candidates",
         x = "Enrichment Ratio (>1 = enriched)",
         y = "GO Term",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_minimal()
  
  # Save the bubble plot
  ggsave("Overlap_Selection_GO_enrichment_bubble_plot.pdf", bubble_plot_overlap, width = 12, height = 8)
}

# 0 significant GO term for overlap


#===============================================================================
# STEP 3: pRDA-ONLY GO ANALYSIS (Environmental Selection Pressures)
#===============================================================================
# Load pRDA candidate SNPs
snp_data_prda <- read.csv("./USA_new/prda_SLF_annotation.csv")

# CRITICAL FIX: Extract pRDA-ONLY genes (exclude overlap with PCAdapt)
# This matches your previous successful approach

# First, get all pRDA genes
all_prda_genes <- snp_data_prda %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

# Get PCAdapt genes for exclusion
all_pcadapt_genes <- final_results %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

# Extract pRDA-ONLY genes (exclude those also in PCAdapt)
sig_gene_list_prda <- setdiff(all_prda_genes, all_pcadapt_genes)

cat("Total pRDA genes:", length(all_prda_genes), "\n")
cat("PCAdapt genes to exclude:", length(all_pcadapt_genes), "\n")
cat("pRDA-ONLY genes (after removing PCAdapt overlap):", length(sig_gene_list_prda), "\n")

cat("pRDA-ONLY significant genes:", length(sig_gene_list_prda), "\n")

# Use same background genes as other analyses
all_gene_ids_prda <- all_gene_ids_pcadapt  # Same background for consistency

cat("pRDA-ONLY significant genes in background:", sum(sig_gene_list_prda %in% all_gene_ids_prda), "\n")

# Create a binary vector where 1 = significant gene, 0 = not significant
gene_vector_prda <- numeric(length(all_gene_ids_prda))
names(gene_vector_prda) <- all_gene_ids_prda

# Populate the gene vector (EXACT REPLICATION of your working method)
for (i in 1:length(all_gene_ids_prda)) {
  gene_id <- all_gene_ids_prda[i]
  # Check if this gene is in our significant gene list
  gene_vector_prda[i] <- ifelse(gene_id %in% sig_gene_list_prda, 1, 0)
}

cat("pRDA-ONLY genes marked as significant in vector:", sum(gene_vector_prda), "\n")

# Create gene2go dataframe for goseq (SAME AS BEFORE)
gene2go_prda <- gene2go_pcadapt  # Same mapping

# Create bias values with sufficient uniqueness for goseq (EXACT REPLICATION)
set.seed(789) # Different seed for independence
bias_prda <- rank(as.numeric(factor(names(gene_vector_prda)))) + 
  runif(length(gene_vector_prda), min=0, max=0.1)
names(bias_prda) <- names(gene_vector_prda)

# Create probability weighting function with permissive parameters (EXACT REPLICATION)
pwf_prda <- nullp(gene_vector_prda, bias.data = bias_prda, plot.fit = FALSE)

# Run GO enrichment analysis (EXACT REPLICATION)
GO_results_prda <- goseq(pwf_prda, gene2cat = gene2go_prda, use_genes_without_cat = TRUE)

# Add FDR correction and process results (EXACT REPLICATION)
GO_results_all_prda <- GO_results_prda %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))

# Save all results
write.csv(GO_results_all_prda, "pRDA_ONLY_USA_GO_enrichment_all_results.csv", row.names = FALSE)

# Filter at a higher threshold to capture more results (EXACT REPLICATION)
GO_results_filtered_prda <- GO_results_all_prda %>%
  filter(FDR < 0.05) %>%  # More lenient threshold
  arrange(over_represented_pvalue)

cat("pRDA-ONLY: Number of terms with FDR < 0.05:", nrow(GO_results_filtered_prda), "\n")

if (nrow(GO_results_filtered_prda) > 0) {
  # Generate visualization of GO enrichment results (EXACT REPLICATION)
  top_GO_prda <- GO_results_filtered_prda %>%
    arrange(over_represented_pvalue) %>%
    head(30)
  
  # Create bar plot of enriched GO terms (EXACT REPLICATION)
  enrichment_plot_prda <- ggplot(top_GO_prda, 
                                 aes(x = reorder(term, -over_represented_pvalue), 
                                     y = -log10(over_represented_pvalue), fill = ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("MF" = "#4daf4a", "BP" = "#3690c0", "CC" = "#c861a5")) + 
    labs(title = "pRDA-ONLY Over-represented GO Terms (Environmental Selection)",
         x = NULL,
         y = "-log10(p-value)",
         fill = "Ontology") +
    theme_minimal()
  
  # Save the GO term enrichment plot
  ggsave("pRDA_ONLY_USA_Selection_GO_top_terms.pdf", enrichment_plot_prda, width = 12, height = 8)
  
  # Create a bubble plot with enrichment information (EXACT REPLICATION)
  top_GO_prda <- top_GO_prda %>%
    mutate(enrichment_ratio = numDEInCat / 
             (numInCat * sum(gene_vector_prda) / length(gene_vector_prda)))
  
  bubble_plot_prda <- ggplot(top_GO_prda, 
                             aes(x = enrichment_ratio, 
                                 y = reorder(term, enrichment_ratio),
                                 size = numDEInCat,
                                 color = -log10(over_represented_pvalue))) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(title = "pRDA-ONLY GO Term Enrichment (Environmental Selection)",
         x = "Enrichment Ratio (>1 = enriched)",
         y = "GO Term",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_minimal()
  
  # Save the bubble plot
  ggsave("pRDA_ONLY_USA_Selection_GO_enrichment_bubble_plot.pdf", bubble_plot_prda, width = 12, height = 8)
}





























#===============================================================================
# STREAMLINED GO ANALYSIS USING YOUR EXACT WORKING APPROACH
# Uses: merged_overlap, final_results, snp_data_prda + GeneID + GeneName background
#===============================================================================

library(dplyr)
library(ggplot2)
library(goseq)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")

# Load your datasets (these are confirmed good)
final_results <- read.csv("./USA_new/Final_overlap_K2-K6_final.csv")
merged_overlap <- read.csv("./USA_new/Overlap_SNPs_pcadapt_prda.csv") 
snp_data_prda <- read.csv("./USA_new/prda_SLF_annotation.csv")

# Load gene-GO mapping
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection")
gene_id_go_mapping <- read.csv("GeneID_GeneName_GOterms.csv")

#===============================================================================
# YOUR EXACT GO ANALYSIS FUNCTION (No modifications)
#===============================================================================

perform_go_analysis <- function(gene_list, file_prefix, set_name) {
  
  # Print the number of genes with annotations
  cat("Number of genes with annotations in", set_name, ":", length(gene_list), "\n")
  
  # Extract all genes that have GO annotations as background (YOUR EXACT APPROACH)
  all_genes <- gene_id_go_mapping %>%
    select(GeneID, GeneName) %>%
    distinct()
  
  # Get all unique gene identifiers as background (GeneID + GeneName)
  all_gene_ids <- unique(c(all_genes$GeneID, all_genes$GeneName))
  
  cat("Background genes (GeneID + GeneName):", length(all_gene_ids), "\n")
  cat("Input genes found in background:", sum(gene_list %in% all_gene_ids), "\n")
  
  # Create a binary vector where 1 = gene from set, 0 = not in set
  gene_vector <- numeric(length(all_gene_ids))
  names(gene_vector) <- all_gene_ids
  
  # Populate the gene vector
  for (i in 1:length(all_gene_ids)) {
    gene_id <- all_gene_ids[i]
    gene_vector[i] <- ifelse(gene_id %in% gene_list, 1, 0)
  }
  
  cat("Genes marked as significant in vector:", sum(gene_vector), "\n")
  
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
  
  # Run GO enrichment analysis
  GO_results <- goseq(pwf, gene2cat = gene2go, use_genes_without_cat = TRUE)
  
  # Add FDR correction and process results
  GO_results_all <- GO_results %>%
    as.data.frame() %>%
    mutate(FDR = p.adjust(over_represented_pvalue, method = "BH"))
  
  GO_results_all <- GO_results_all %>%
    filter(!is.na(term) & term != "NA")
  
  # Save all results
  write.csv(GO_results_all, paste0(file_prefix, "_GO_enrichment_all_results.csv"), row.names = FALSE)
  
  # Filter significant results
  GO_results_filtered <- GO_results_all %>%
    filter(FDR < 0.05) %>%
    arrange(over_represented_pvalue)
  
  # Print the number of significant terms
  cat("Number of GO terms with FDR < 0.05 in", set_name, ":", nrow(GO_results_filtered), "\n")
  
  # Define top_GO for plotting
  top_GO <- if(nrow(GO_results_filtered) > 0) {
    GO_results_filtered %>%
      head(min(15, nrow(GO_results_filtered)))
  } else {
    GO_results_all %>%
      arrange(over_represented_pvalue) %>%
      head(30)
  }
  
  # Create plots if we have data
  if(nrow(top_GO) > 0) {
    # Bar plot
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
    
    ggsave(paste0(file_prefix, "_GO_top_terms.pdf"), enrichment_plot, width = 12, height = 8)
    
    # Bubble plot
    top_GO <- top_GO %>%
      mutate(enrichment_ratio = numDEInCat / (numInCat * sum(gene_vector) / length(gene_vector)))
    
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
  
  # Save detailed table
  top_GO_table <- top_GO %>%
    select(term, ontology, numDEInCat, numInCat, over_represented_pvalue, FDR, enrichment_ratio) %>%
    arrange(over_represented_pvalue)
  
  write.csv(top_GO_table, paste0(file_prefix, "_top_GO_terms_details.csv"), row.names = FALSE)
  
  return(GO_results_all %>% mutate(is_significant = FDR < 0.05))
}

#===============================================================================
# CREATE GENE LISTS FROM YOUR DATASETS
#===============================================================================

# Extract genes from overlap (both methods)
overlap_genes <- merged_overlap %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

# Extract genes from PCAdapt
pcadapt_all_genes <- final_results %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

# Extract genes from pRDA  
prda_all_genes <- snp_data_prda %>%
  filter(!is.na(Gene_Name), Gene_Name != "") %>%
  pull(Gene_Name) %>%
  unique()

# Create method-specific gene lists (YOUR EXACT APPROACH)
pcadapt_only_genes <- setdiff(pcadapt_all_genes, overlap_genes)
prda_only_genes <- setdiff(prda_all_genes, overlap_genes)

cat("Gene counts:\n")
cat("- Overlap genes:", length(overlap_genes), "\n")
cat("- PCAdapt-only genes:", length(pcadapt_only_genes), "\n")
cat("- pRDA-only genes:", length(prda_only_genes), "\n")

#===============================================================================
# RUN GO ANALYSIS FOR EACH GENE SET
#===============================================================================

# Change to output directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection")

cat("\n=== RUNNING GO ANALYSES ===\n")

# Run analyses
overlap_go <- perform_go_analysis(overlap_genes, "overlap", "Overlap SNPs")
pcadapt_only_go <- perform_go_analysis(pcadapt_only_genes, "pcadapt_only", "PCAdapt-only SNPs")  
prda_only_go <- perform_go_analysis(prda_only_genes, "prda_only", "pRDA-only SNPs")

#===============================================================================
# RESULTS SUMMARY
#===============================================================================
overlap_sig <- if(!is.null(overlap_go)) sum(overlap_go$is_significant) else 0
pcadapt_sig <- if(!is.null(pcadapt_only_go)) sum(pcadapt_only_go$is_significant) else 0
prda_sig <- if(!is.null(prda_only_go)) sum(prda_only_go$is_significant) else 0

cat("Significant GO terms (FDR < 0.05):\n")
cat("- Overlap:", overlap_sig, "\n")
cat("- PCAdapt-only:", pcadapt_sig, "\n") 
cat("- pRDA-only:", prda_sig, "\n")

# Save gene lists for reference
write.csv(data.frame(Gene = overlap_genes), "overlap_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = pcadapt_only_genes), "pcadapt_only_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = prda_only_genes), "prda_only_genes.csv", row.names = FALSE)

# Create summary
summary_df <- data.frame(
  Category = c("Overlap", "PCAdapt-only", "pRDA-only"),
  Gene_Count = c(length(overlap_genes), length(pcadapt_only_genes), length(prda_only_genes)),
  Significant_GO_Terms = c(overlap_sig, pcadapt_sig, prda_sig)
)

write.csv(summary_df, "GO_analysis_summary.csv", row.names = FALSE)
print(summary_df)


#===============================================================================
# CREATE PUBLICATION-QUALITY DOTPLOT FOR SLF INVASION GO TERMS
# Replicating your successful visualization style
#===============================================================================

library(ggplot2)
library(dplyr)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection")

# Load the GO results that gave you the successful results
overlap_go <- read.csv("overlap_GO_enrichment_all_results.csv")
pcadapt_only_go <- read.csv("pcadapt_only_GO_enrichment_all_results.csv")
prda_only_go <- read.csv("prda_only_GO_enrichment_all_results.csv")

#===============================================================================
# PREPARE DATA FOR DOTPLOT VISUALIZATION
#===============================================================================

# Create comprehensive dataset with all significant GO terms from all categories
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

cat("Significant GO terms found:\n")
cat("- Overlap:", sum(all_sig_terms$Category == "Overlap"), "\n")
cat("- PCAdapt-only:", sum(all_sig_terms$Category == "PCAdapt-only"), "\n") 
cat("- pRDA-only:", sum(all_sig_terms$Category == "pRDA-only"), "\n")
cat("- Total unique terms:", length(unique(all_sig_terms$term)), "\n")

if (nrow(all_sig_terms) == 0) {
  cat("No significant terms found with FDR < 0.05. Check your input files.\n")
  stop("No data to plot")
}

# Create a copy of the terms ordered by p-value for y-axis ordering (most significant at top)
ordered_terms <- all_sig_terms %>%
  arrange(over_represented_pvalue) %>%  # Smallest p-values first = most significant
  pull(term) %>%
  unique() %>%
  rev()  # This reverses the order so the most significant appears at the top

# Make term a factor to preserve the p-value ordering across all categories
all_sig_terms$term <- factor(all_sig_terms$term, levels = ordered_terms)

# Set Category as a factor with the desired order (Overlap first, then PCAdapt, then pRDA)
all_sig_terms$Category <- factor(all_sig_terms$Category, 
                                 levels = c("Overlap", "PCAdapt-only", "pRDA-only"))

# Save the data for reference
write.csv(all_sig_terms, "all_sig_terms_SLF_invasion.csv", row.names = FALSE)

#===============================================================================
# CREATE THE DOTPLOT VISUALIZATION
#===============================================================================

# Calculate plot height based on number of terms
plot_height <- max(8, nrow(all_sig_terms) * 0.25)

# Create the main dotplot
combined_plot <- ggplot(all_sig_terms,
                        aes(x = Category, y = term,
                            size = numDEInCat,
                            color = -log10(over_represented_pvalue))) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
  scale_size_continuous(range = c(2, 8), name = "Gene Count") +
  labs(title = "Functional Signatures of SLF Invasion: GO Term Enrichment",
       subtitle = paste( "PCAdapt-only:", sum(all_sig_terms$Category == "PCAdapt-only"),
                        "terms | pRDA-only:", sum(all_sig_terms$Category == "pRDA-only"), "terms"),
       x = NULL,
       y = "GO Terms (ordered by significance)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Display the plot
print(combined_plot)

# Save high-quality versions
ggsave("SLF_Invasion_GO_terms_dotplot.pdf", combined_plot, 
       width = 12, height = plot_height, dpi = 300)
ggsave("SLF_Invasion_GO_terms_dotplot.png", combined_plot, 
       width = 12, height = plot_height, dpi = 400)






####### To just do K=11
library(pcadapt)
library(qvalue)
library(dplyr)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")
# Load data
bed_file <- "./USA_final_annotated.bed"
annot_data <- annot_combined  # assuming this is already loaded
pop_data <- read.csv("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/USA_SLF_MS_Info.csv")

# Read .bed file
genetic_data <- read.pcadapt(bed_file, type = "bed")

# Run pcadapt with K=11
pcadapt_k11 <- pcadapt(input = genetic_data, K = 11, method = "mahalanobis", ploidy = 2)

# Calculate q-values
qv_k11 <- qvalue(pcadapt_k11$pvalues)
qval_adj_k11 <- qv_k11$qvalues

# Get SNPs passing FDR threshold
fdr_threshold <- 0.01
signif_indices_k11 <- which(qval_adj_k11 < fdr_threshold)

cat("✅ Number of SNPs with FDR < 0.01 at K=11:", length(signif_indices_k11), "\n")

# Extract info for significant SNPs
sig_snps_k11 <- annot_data[signif_indices_k11, ]
sig_snps_k11$p_value <- pcadapt_k11$pvalues[signif_indices_k11]
sig_snps_k11$q_value <- qval_adj_k11[signif_indices_k11]

# Save
write.csv(sig_snps_k11, "./USA_new/pcadapt_K11_significant_snps.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)

# Load data
snp_data <- sig_snps_k11
snp_data$chr_num <- as.numeric(gsub("chr_", "", snp_data$CHROM))
snp_data <- snp_data %>% arrange(chr_num, POS)
snp_data$chr_color <- factor(snp_data$chr_num %% 2)

# Create cumulative positions
chr_pos <- snp_data %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS), chr_center = mean(POS), .groups = "drop") %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

chr_pos$center_cum <- 0
chr_cum_pos <- 0
snp_data$pos_cum <- 0

for (i in 1:nrow(chr_pos)) {
  chr_snps <- snp_data$CHROM == chr_pos$CHROM[i]
  snp_data$pos_cum[chr_snps] <- snp_data$POS[chr_snps] + chr_cum_pos
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5e6
}

# Manhattan Plot
p <- ggplot(snp_data, aes(x = pos_cum, y = -log10(p_value), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#3690c0")) +
  scale_x_continuous(breaks = chr_pos$center_cum, labels = chr_pos$CHROM) +
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
    title = "Manhattan Plot: pcadapt K = 11",
    subtitle = paste0("FDR < 0.01: ", nrow(snp_data), " significant SNPs")
  )

# Save
ggsave("./USA_new/SLF_USA_Manhattan_K11.pdf", p, width = 14, height = 8, dpi = 300)
ggsave("./USA_new/SLF_USA_Manhattan_K11.png", p, width = 14, height = 8, dpi = 300)
print(p)








setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")

# Filter for chromosomes of interest
chromosomes_of_interest <- c("chr_1", "chr_6", "chr_8")
filtered_df <- final_results %>%
  filter(CHROM %in% chromosomes_of_interest)

# Loop through each chromosome and plot
for (chrom in chromosomes_of_interest) {
  chrom_data <- filtered_df %>% filter(CHROM == chrom)
  
  p <- ggplot(chrom_data, aes(x = POS, y = -log10(FDR_p_value))) +
    geom_point(alpha = 0.7, color = "#2b8cbe") +
    labs(title = paste("Detailed Manhattan Plot -", chrom),
         x = "Genomic Position",
         y = expression(-log[10](italic(p)))) +
    theme_minimal()
  
  ggsave(paste0("USA_LD/Manhattan_", chrom, ".pdf"), plot = p, width = 10, height = 4)
}

# Count SNPs per chromosome
chr_counts <- snp_data_prda %>%
  group_by(CHROM) %>%
  summarise(n_snps = n()) %>%
  mutate(total_snps = sum(n_snps),
         proportion = round(n_snps / total_snps, 4)) %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

# Print the result
print(chr_counts)



# Filter for chr_1, chr_6, chr_8
top_snps_chr_168 <- final_results %>%
  filter(CHROM %in% c("chr_1", "chr_6", "chr_8"))

# Count annotation types
annotation_counts <- top_snps_chr_168 %>%
  mutate(is_intron_or_intergenic = Annotation %in% c("intron_variant", "intergenic_region")) %>%
  group_by(CHROM, Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(CHROM, desc(n))

# Total number and percentage of intron/intergenic
summary_table <- top_snps_chr_168 %>%
  mutate(annotation_type = ifelse(Annotation %in% c("intron_variant", "intergenic_region"), "intron/intergenic", "other")) %>%
  group_by(CHROM, annotation_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = annotation_type, values_from = count, values_fill = 0) %>%
  mutate(total = `intron/intergenic` + other,
         percent_intron_intergenic = round(`intron/intergenic` / total * 100, 2))

# View results
print(annotation_counts)
print(summary_table)







######## Explore Combined (USA+SH)
# ----------------- SETUP -----------------
library(pcadapt)
library(qvalue)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection")

# ----------------- LOAD DATA -----------------
bed_file <- "./Final_SLF.bed"
annot_data <- annot_combined  # should have CHROM, POS, etc.
pop_data <- read.csv("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/Final_SLF_MS_Info.csv")

# Read .bed
genetic_data <- read.pcadapt(bed_file, type = "bed")

# ----------------- SCREE PLOT TO CHOOSE K -----------------
scree_res <- pcadapt(input = genetic_data, K = 20, method = "mahalanobis", ploidy = 2)
scree_df <- data.frame(K = 1:20, eigenvalues = scree_res$singular.values^2)

ggplot(scree_df, aes(x = K, y = eigenvalues)) +
  geom_line(color = "darkblue") +
  geom_point(size = 2) +
  labs(title = "Scree Plot: Combined SH + USA", x = "K", y = "Eigenvalue") +
  theme_minimal(base_size = 14)
ggsave("./Combined_SLF_ScreePlot.pdf", width = 8, height = 6)

# ----------------- CHOOSE K (MANUALLY) -----------------
# Based on elbow of scree plot; suppose you choose K = 5
K_opt <- 5
pcadapt_res <- pcadapt(input = genetic_data, K = K_opt, method = "mahalanobis", ploidy = 2)

# ----------------- CALCULATE Q-VALUES -----------------
qv <- qvalue(pcadapt_res$pvalues)
qval_adj <- qv$qvalues

# Get significant SNPs
fdr_threshold <- 0.01
signif_indices <- which(qval_adj < fdr_threshold)
cat("✅ Significant SNPs (FDR < 0.01):", length(signif_indices), "\n")

# ----------------- EXTRACT SIGNIFICANT SNP INFO -----------------
sig_snps <- annot_data[signif_indices, ]
sig_snps$p_value <- pcadapt_res$pvalues[signif_indices]
sig_snps$q_value <- qval_adj[signif_indices]

# Save significant SNPs
write.csv(sig_snps, "./Combined_SLF/pcadapt_K5_significant_snps.csv", row.names = FALSE)

# ----------------- MANHATTAN PLOT -----------------
sig_snps$chr_num <- as.numeric(gsub("chr_", "", sig_snps$CHROM))
sig_snps <- sig_snps %>% arrange(chr_num, POS)
sig_snps$chr_color <- factor(sig_snps$chr_num %% 2)

# Create cumulative positions
chr_pos <- sig_snps %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS), chr_center = mean(POS), .groups = "drop") %>%
  arrange(as.numeric(gsub("chr_", "", CHROM)))

chr_pos$center_cum <- 0
chr_cum_pos <- 0
sig_snps$pos_cum <- 0

for (i in 1:nrow(chr_pos)) {
  chr_snps <- sig_snps$CHROM == chr_pos$CHROM[i]
  sig_snps$pos_cum[chr_snps] <- sig_snps$POS[chr_snps] + chr_cum_pos
  chr_pos$center_cum[i] <- chr_cum_pos + chr_pos$chr_center[i]
  chr_cum_pos <- chr_cum_pos + chr_pos$chr_len[i] + 5e6
}

# Plot
manhattan_plot <- ggplot(sig_snps, aes(x = pos_cum, y = -log10(p_value), color = chr_color)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3690c0", "1" = "#3690c0")) +
  scale_x_continuous(breaks = chr_pos$center_cum, labels = chr_pos$CHROM) +
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
    title = paste("Manhattan Plot: pcadapt K =", K_opt),
    subtitle = paste0("FDR < 0.01: ", nrow(sig_snps), " significant SNPs")
  )

# Save plot
ggsave("./Combined_SLF/SLF_Combined_Manhattan_K5.pdf", manhattan_plot, width = 14, height = 8, dpi = 300)
ggsave("./Combined_SLF/SLF_Combined_Manhattan_K5.png", manhattan_plot, width = 14, height = 8, dpi = 300)
print(manhattan_plot)








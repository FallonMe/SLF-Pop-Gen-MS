# Load required libraries
library(readr)   
library(data.table)
library(dplyr)
library(ggplot2)
library(vegan)

# Set working directory if needed
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/Overlap")


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

# Load RDA results
print("Loading RDA data...")
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Final_USA_RDA.RData")

# First, read in the file with the 184 shared SNPs
print("Reading shared SNPs file...")
shared_snps <- read.csv("184overlap_prda_pcadapt_full_info.csv")

# Extract just the SNP names from the shared SNPs file
shared_snp_ids <- shared_snps$SNP

print(paste("Number of shared SNPs:", length(shared_snp_ids)))

# Assign colors based on environmental predictor
sel_prda <- cand_prda$snp
env_prda <- cand_prda$predictor

# Define colors
env_prda[env_prda == "Temperature"]   <- "#1f78b4"  # Blue
env_prda[env_prda == "Precipitation"] <- "#c861a5"  # Purple
env_prda[env_prda == "Urbanization"]  <- "#33a02c"  # Green

# Initialize colors for all SNPs as gray (non-candidates)
col_pred_prda <- rep('#f1eef6', length(rownames(prda_rda$CCA$v)))  # Default gray
names(col_pred_prda) <- rownames(prda_rda$CCA$v)  # Assign row names

# Initialize shapes vector - default is circle (pch=21)
shape_pred_prda <- rep(21, length(rownames(prda_rda$CCA$v)))
names(shape_pred_prda) <- rownames(prda_rda$CCA$v)

# Process SNPs in batches to improve performance
print("Processing SNPs and assigning colors/shapes...")
batch_size <- 1000
total_snps <- length(sel_prda)
num_batches <- ceiling(total_snps / batch_size)

for (batch in 1:num_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, total_snps)
  
  print(paste("Processing batch", batch, "of", num_batches, 
              "(SNPs", start_idx, "to", end_idx, ")"))
  
  for (i in start_idx:end_idx) {
    foo <- match(sel_prda[i], names(col_pred_prda))  # Ensure proper matching
    if (!is.na(foo)) {
      col_pred_prda[foo] <- env_prda[i]  # Assign correct color
      
      # Check if this SNP is in the shared list
      snp_base <- gsub("_[ACGT]$", "", sel_prda[i])
      
      if (snp_base %in% shared_snp_ids) {
        # If it's a shared SNP, change to triangle shape (pch=24)
        shape_pred_prda[foo] <- 24
      }
    }
  }
}

# Alternative approach: Add black outline/border to the shared SNPs
print("Setting up border colors...")
col_border_prda <- rep("gray32", length(rownames(prda_rda$CCA$v)))  # Default gray border
names(col_border_prda) <- rownames(prda_rda$CCA$v)

# Process shared SNPs for borders
for (batch in 1:num_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, total_snps)
  
  for (i in start_idx:end_idx) {
    foo <- match(sel_prda[i], names(col_border_prda))
    if (!is.na(foo)) {
      # Check if this is a shared SNP
      snp_base <- gsub("_[ACGT]$", "", sel_prda[i])
      if (snp_base %in% shared_snp_ids) {
        col_border_prda[foo] <- "black"  # Black border for shared SNPs
      }
    }
  }
}

# Set non-candidates to transparent
empty <- col_pred_prda
empty[empty == "#f1eef6"] <- rgb(0,1,0, alpha=0)  # Transparent for non-candidates
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")

# Define legend colors
bg_prda <- c("#1f78b4", "#c861a5", "#33a02c")

# Calculate loadings once to avoid repeating in the plot loops
print("Calculating SNP loadings...")
loadings_1_2 <- scores(prda_rda, choices=c(1,2), display="species", scaling=1)
loadings_1_3 <- scores(prda_rda, choices=c(1,3), display="species", scaling=1)

# Get candidate SNP names for faster filtering
candidate_snps <- names(col_pred_prda)[col_pred_prda != '#f1eef6']
print(paste("Number of candidate SNPs to plot:", length(candidate_snps)))

# OPTION 1: Plot with different shapes for shared SNPs
print("Creating plot with different shapes (axes 1 & 2)...")
pdf("Partial_RDA_SNP_shared_shapes_1_2.pdf", width = 10, height = 8)
plot(prda_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1),
     xlab= 'RDA1 (37.0%)', ylab= 'RDA2 (32.3%)')

# Only plot candidate SNPs
for (snp_name in candidate_snps) {
  snp_coords <- loadings_1_2[snp_name,]
  points(snp_coords[1], snp_coords[2], 
         pch=shape_pred_prda[snp_name], 
         bg=col_pred_prda[snp_name], 
         col="gray32", 
         cex=1.5)
}

text(prda_rda, scaling=1, display="bp", col="#191970", cex=1.2)

# Add to legend for the shared SNPs
legend("bottomright", 
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"), 
       bty = "n", 
       col = c("gray32", "gray32", "gray32", "gray32"), 
       pch = c(21, 21, 21, 24), 
       cex = 1.5, 
       pt.bg = c(bg_prda, "gray"))
dev.off()
print("Plot 1 completed.")

# OPTION 2: Plot with black outlines for shared SNPs
print("Creating plot with black borders (axes 1 & 2)...")
pdf("Partial_RDA_SNP_shared_borders_1_2.pdf", width = 10, height = 8)
plot(prda_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1),
     xlab= 'RDA1 (37.0%)', ylab= 'RDA2 (32.3%)')

# Only plot candidate SNPs
for (snp_name in candidate_snps) {
  snp_coords <- loadings_1_2[snp_name,]
  
  # Determine if it's a shared SNP for border width
  snp_base <- gsub("_[ACGT]$", "", snp_name)
  lwd_val <- ifelse(snp_base %in% shared_snp_ids, 2, 1)
  
  points(snp_coords[1], snp_coords[2], 
         pch=21, 
         bg=col_pred_prda[snp_name], 
         col=col_border_prda[snp_name], 
         cex=1.5,
         lwd=lwd_val)
}

text(prda_rda, scaling=1, display="bp", col="#191970", cex=1.2)

# Add legend including shared SNPs with black border
legend("bottomright", 
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"), 
       bty = "n", 
       col = c("gray32", "gray32", "gray32", "black"), 
       pch = c(21, 21, 21, 21), 
       cex = 1.5, 
       pt.bg = c(bg_prda, "gray"),
       lwd = c(1, 1, 1, 2))
dev.off()
print("Plot 2 completed.")

# Repeat the same approach for axes 1 & 3
print("Creating plot with different shapes (axes 1 & 3)...")
pdf("Partial_RDA_SNP_shared_shapes_1_3.pdf", width = 10, height = 8)
plot(prda_rda, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1), 
     xlab= 'RDA1 (32.9%)', ylab= 'RDA3 (30.1%)', choices=c(1,3))

# Only plot candidate SNPs
for (snp_name in candidate_snps) {
  snp_coords <- loadings_1_3[snp_name,]
  points(snp_coords[1], snp_coords[2], 
         pch=shape_pred_prda[snp_name], 
         bg=col_pred_prda[snp_name], 
         col="gray32", 
         cex=1.5)
}

text(prda_rda, scaling=1, display="bp", col="#191970", cex=1.2, choices=c(1,3))

# Add to legend for the shared SNPs
legend("bottomright", 
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"), 
       bty = "n", 
       col = c("gray32", "gray32", "gray32", "gray32"), 
       pch = c(21, 21, 21, 24), 
       cex = 1.5, 
       pt.bg = c(bg_prda, "gray"))
dev.off()
print("Plot 3 completed.")


# Option to export the SNP loadings for further visualization with ggplot
print("Preparing SNP loadings for ggplot...")
snp_loadings <- scores(prda_rda, choices = c(1:3), display = "species", scaling = 1)
snp_loadings_df <- as.data.frame(snp_loadings)

# Add a column identifying shared SNPs
snp_loadings_df$shared <- "no"
rownames_base <- gsub("_[ACGT]$", "", rownames(snp_loadings_df))
snp_loadings_df$shared[rownames_base %in% shared_snp_ids] <- "yes"

# Add a column identifying environment association
snp_loadings_df$environment <- "non-candidate"

# ✅ Define candidate SNPs properly
candidate_snps <- sel_prda[sel_prda %in% rownames(snp_loadings_df)]

# Assign environment labels to candidate SNPs
for (snp_name in candidate_snps) {
  i <- which(sel_prda == snp_name)
  if (length(i) > 0) {
    snp_loadings_df[snp_name, "environment"] <- as.character(env_prda[i[1]])
  }
}


# Option to export the SNP loadings for further visualization with ggplot
print("Preparing SNP loadings for ggplot...")
snp_loadings <- scores(prda_rda, choices=c(1:3), display="species", scaling=1)
snp_loadings_df <- as.data.frame(snp_loadings)

# Add a column identifying shared SNPs
snp_loadings_df$shared <- "no"
rownames_base <- gsub("_[ACGT]$", "", rownames(snp_loadings_df))
snp_loadings_df$shared[rownames_base %in% shared_snp_ids] <- "yes"

# Add a column identifying environment association
snp_loadings_df$environment <- "non-candidate"
for (snp_name in candidate_snps) {
  i <- which(sel_prda == snp_name)
  if (length(i) > 0) {
    snp_loadings_df[snp_name,"environment"] <- as.character(env_prda[i[1]])
  }
}

# Save the loadings data for potential ggplot use
write.csv(snp_loadings_df, "snp_loadings_with_shared_info.csv")
print("SNP loadings saved to CSV.")



library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# Load your files
snp_df <- read_csv("snp_loadings_with_shared_info.csv")
shared_snps <- read_csv("184overlap_prda_pcadapt_full_info.csv")

# Extract just chr_pos_ref from your SNP names in snp_df
snp_df <- snp_df %>%
  mutate(
    base_id = str_replace(...1, "_[ACGT]$", ""),  # remove trailing _ALT
    env_label = case_when(
      environment == "#1f78b4" ~ "Temperature",
      environment == "#c861a5" ~ "Precipitation",
      environment == "#33a02c" ~ "Urbanization",
      TRUE ~ NA_character_
    )
  )

# Match based on base SNP IDs
shared_ids <- shared_snps$SNP  # these are chr_pos_ref
snp_df <- snp_df %>%
  mutate(is_shared = ifelse(base_id %in% shared_ids, "yes", "no"))

# Now filter candidates
cand_df <- snp_df %>% filter(!is.na(env_label))

# Check how many shared SNPs are now matched
cat("✅ Shared SNPs among candidates:", nrow(filter(cand_df, is_shared == "yes")), "\n")

library(ggplot2)
library(dplyr)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Selection/Final_selection/Overlap")
cand_df <- read.csv("RDA_plot.csv") 

# Add a shape flag to distinguish shared SNPs
cand_df <- cand_df %>%
  mutate(
    outline = ifelse(is_shared == "yes", "Shared", "Not Shared")
  )

# write.csv(cand_df, "RDA_plot.csv")

# Plot
p <- ggplot(cand_df, aes(x = RDA1, y = RDA3)) +
  # All SNPs, shape = 21 for filled points
  geom_point(aes(fill = env_label), color = "grey30", shape = 21, size = 2, alpha = 0.8, stroke = 0.3) +
  
  # Overlay shared SNPs with thicker black outline
  geom_point(data = filter(cand_df, outline == "Shared"),
             aes(fill = env_label),
             shape = 21,
             size = 2.5,
             color = "black",
             stroke = 1.2) +
  
  # Fill color by environmental driver
  scale_fill_manual(
    name = "Environmental Predictor",
    values = c(
      "Temperature" = "#1f78b4",
      "Precipitation" = "#c861a5",
      "Urbanization" = "#33a02c"
    )
  ) +
  
  # Legends: no separate shape legend
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "grey30", stroke = 0.3, size = 3),
      title = "Environmental Predictor"
    )
  ) +
  
  # Axis labels
  labs(
    x = "RDA1 (37.0%)",
    y = "RDA3 (30.1%)"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save
ggsave("RDA1_RDA3_final.pdf", plot = p, width = 10, height = 8)
ggsave("RDA1_RDA3_final.png", plot = p, width = 10, height = 8, dpi = 600)

print(p)

# Plot
p <- ggplot(cand_df, aes(x = RDA1, y = RDA2)) +
  # All SNPs, shape = 21 for filled points
  geom_point(aes(fill = env_label), color = "grey30", shape = 21, size = 2, alpha = 0.8, stroke = 0.3) +
  
  # Overlay shared SNPs with thicker black outline
  geom_point(data = filter(cand_df, outline == "Shared"),
             aes(fill = env_label),
             shape = 21,
             size = 2.5,
             color = "black",
             stroke = 1.2) +
  
  # Fill color by environmental driver
  scale_fill_manual(
    name = "Environmental Predictor",
    values = c(
      "Temperature" = "#1f78b4",
      "Precipitation" = "#c861a5",
      "Urbanization" = "#33a02c"
    )
  ) +
  
  # Legends: no separate shape legend
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "grey30", stroke = 0.3, size = 3),
      title = "Environmental Predictor"
    )
  ) +
  
  # Axis labels
  labs(
    x = "RDA1 (37.0%)",
    y = "RDA2 (32.3%)"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save
ggsave("RDA1_RDA2_final.pdf", plot = p, width = 10, height = 8)
ggsave("RDA1_RDA2_final.png", plot = p, width = 10, height = 8, dpi = 600)

print(p)



# Add shape flag if not already done
cand_df <- cand_df %>%
  mutate(outline = ifelse(is_shared == "yes", "Shared", "Not Shared"))

# Plot
p <- ggplot(cand_df, aes(x = RDA1, y = RDA2)) +
  # All SNPs, shape = 21 for filled points
  geom_point(aes(fill = env_label), color = "grey30", shape = 21, size = 2, alpha = 0.8, stroke = 0.3) +
  
  # Overlay shared SNPs with thicker black outline
  geom_point(data = filter(cand_df, outline == "Shared"),
             aes(fill = env_label),
             shape = 21,
             size = 2.5,
             color = "black",
             stroke = 1.2) +
  
  # Fill color by environmental driver
  scale_fill_manual(
    name = "Environmental Predictor",
    values = c(
      "Temperature" = "#1f78b4",
      "Precipitation" = "#c861a5",
      "Urbanization" = "#33a02c"
    )
  ) +
  
  # Legend settings
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "grey30", stroke = 0.3, size = 3),
      title = "Environmental Predictor"
    )
  ) +
  
  # Axis labels
  labs(
    x = "RDA1 (37.0%)",
    y = "RDA2 (32.3%)"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save
ggsave("RDA1_RDA2_final.pdf", plot = p, width = 10, height = 8)
ggsave("RDA1_RDA2_final.png", plot = p, width = 10, height = 8, dpi = 600)

# Display plot
print(p)


# Plot: RDA2 vs RDA3
p <- ggplot(cand_df, aes(x = RDA2, y = RDA3)) +
  # All SNPs, shape = 21 for filled points
  geom_point(aes(fill = env_label), color = "grey30", shape = 21, size = 2, alpha = 0.8, stroke = 0.3) +
  
  # Overlay shared SNPs with thicker black outline
  geom_point(data = filter(cand_df, outline == "Shared"),
             aes(fill = env_label),
             shape = 21,
             size = 2.5,
             color = "black",
             stroke = 1.2) +
  
  # Fill color by environmental driver
  scale_fill_manual(
    name = "Environmental Predictor",
    values = c(
      "Temperature" = "#1f78b4",
      "Precipitation" = "#c861a5",
      "Urbanization" = "#33a02c"
    )
  ) +
  
  # Legend settings
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "grey30", stroke = 0.3, size = 3),
      title = "Environmental Predictor"
    )
  ) +
  
  # Axis labels
  labs(
    x = "RDA2 (32.3%)",
    y = "RDA3 (30.1%)"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save
ggsave("RDA2_RDA3_final.pdf", plot = p, width = 10, height = 8)
ggsave("RDA2_RDA3_final.png", plot = p, width = 10, height = 8, dpi = 600)

# Display plot
print(p)


fdr_USA_Only <- merge(overlap_all_k, annot_combined, by = "SNP", all.x = TRUE)
write.csv(fdr_USA_Only, "USA_Only_pcadapt_SigFDR_results.csv", row.names = FALSE)


# Clean up annot_combined$SNP_ID just in case
annot_combined <- annot_combined %>%
  mutate(SNP_ID = as.character(SNP_ID))

# Join annotation to cand_df
cand_df_annotated <- cand_df %>%
  left_join(annot_combined, by = c("base_id" = "SNP"))

# Check result
glimpse(cand_df_annotated)



library(ggrepel)

# Filter overlapping, annotated, and named functional SNPs
highlight_snps <- cand_df_annotated %>%
  filter(
    is_shared == "yes",
    !is.na(Gene_Name),
    Gene_Name != "",
    !grepl("^LDel_", Gene_Name),
    !grepl("-", Gene_Name)
  ) %>%
  group_by(Gene_Name) %>%
  slice_min(order_by = abs(RDA2) + abs(RDA3), n = 1, with_ties = FALSE) %>%
  ungroup()


# Plot
p <- ggplot(cand_df, aes(x = RDA2, y = RDA3)) +
  # All SNPs, shape = 21 for filled points
  geom_point(aes(fill = env_label), color = "grey30", shape = 21, size = 3, alpha = 0.8, stroke = 0.3) +
  
  # Overlay shared SNPs with thicker black outline
  geom_point(data = filter(cand_df, outline == "Shared"),
             aes(fill = env_label),
             shape = 21,
             size = 3.5,
             color = "black",
             stroke = 1.2) +
  
  # Fill color by environmental driver
  scale_fill_manual(
    name = "Environmental Predictor",
    values = c(
      "Temperature" = "#1f78b4",
      "Precipitation" = "#c861a5",
      "Urbanization" = "#33a02c"
    )
  ) +
  
  # Legend settings
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "grey30", stroke = 0.3, size = 3),
      title = "Environmental Predictor"
    )
  ) +
  
  # Axis labels
  labs(
    x = "RDA2 (32.3%)",
    y = "RDA3 (30.1%)"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  )


# Add labels to plot
p <- p +
  geom_text_repel(
    data = highlight_snps,
    aes(x = RDA2, y = RDA3, label = Gene_Name),
    size = 3.0,
    fontface = "bold",
    color = "black",
    box.padding = 0.6,
    point.padding = 0.5,
    max.overlaps = 50,
    segment.size = 0.3,
    segment.color = "grey50"
  )

# Display plot
print(p)

# Save
ggsave("RDA2_RDA3_final_annotated.pdf", plot = p, width = 10, height = 8)
ggsave("RDA2_RDA3_final_annotated.png", plot = p, width = 10, height = 8, dpi = 600)


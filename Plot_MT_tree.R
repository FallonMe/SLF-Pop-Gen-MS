# Phylogenetic visualization for SLF mitochondrial haplotypes
# All trees rooted with Lycorma_meliae_outgroup and matched to sample metadata

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Set working directory and load tree and metadata
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/MT/Final/ALL510")
tree <- read.tree("WG_with_outgroup.contree")
tree <- root(tree, outgroup = "Lycorma_meliae_outgroup", resolve.root = TRUE)

metadata <- read.csv("Final_MT_meta.csv", header = TRUE) %>%
  mutate(across(c(SampleID, Haplotype, Sample_Code, Country, Lineage, City), as.character))

# Match tree tips to metadata
tree_data <- data.frame(label = tree$tip.label) %>%
  left_join(
    metadata %>%
      group_by(Haplotype) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(Haplotype, SampleID, Sample_Code, Country, Lineage, City),
    by = c("label" = "Haplotype")
  ) %>%
  mutate(Label = ifelse(!is.na(Sample_Code), Sample_Code, abbreviate(City, minlength = 4)))

stopifnot(!any(is.na(tree_data$Label)))

# Color schemes
lineage_colors <- c("EP"="#e6ab02", "SH"="#a65628", "SB"="#08589e",
                    "LP"="#c861a5", "GP"="#ff7f00", "YP"="#6a3d9a", "unknown"="gray60")
country_colors <- c("CN"="#6a3d9a", "USA"="#e41a1c", "SK"="#009E73", "JP"="#e6ab02",
                    "Lycorma_meliae_outgroup"="green")

# Rectangular ML tree
p_rect <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") %<+% tree_data +
  geom_tippoint(aes(color = Lineage), size = 1.0, alpha = 0.8) +
  scale_color_manual(name = "Lineage", values = lineage_colors, na.value = "gray80") +
  ggnewscale::new_scale_color() +
  geom_tiplab(aes(label = Label, color = Country), size = 1.8, hjust = -0.05) +
  scale_color_manual(name = "Country", values = country_colors, na.value = "gray50") +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(15, 15, 15, 15),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("ML Phylogeny of SLF (Sample Code Tip Labels)")

ggsave("SLF_ML_tree_branch_lengths_sampleCode.pdf", plot = p_rect, width = 20, height = 18, dpi = 300)

# Radial tree
p_radial <- ggtree(tree, layout = "radial", branch.length = "branch.length") %<+% tree_data +
  geom_tippoint(aes(color = Lineage), size = 1.3, alpha = 0.8) +
  scale_color_manual(name = "Lineage", values = lineage_colors, na.value = "gray60") +
  ggnewscale::new_scale_color() +
  geom_tiplab(aes(label = Label, color = Country), size = 1.8, hjust = -0.1, show.legend = FALSE) +
  scale_color_manual(name = "Country", values = country_colors, na.value = "gray60") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ggtitle("Radial Phylogeny of SLF (Sample Code Tip Labels)")

ggsave("SLF_radial_ML_tree_SampleID.pdf", plot = p_radial, width = 14, height = 14, dpi = 300)


# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pcadapt)
library(ape)
library(phytools)
library(ggtree)
library(ggnewscale)

setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/MT/Final/ALL510")
metadata <- read.csv("Final_MT_meta.csv", header = TRUE)

#with outgroup IQ tree
tree <- read.tree("WG_with_outgroup.contree")

tree <- root(tree, outgroup = "Lycorma_meliae_outgroup", resolve.root = TRUE)

# Clean up the metadata and prepare for joining
metadata_clean <- metadata %>%
  mutate(
    Haplotype_clean = gsub(" ", "_", Haplotype),
    SampleID = as.character(SampleID)
  )

# Extract base haplotype name from tip labels if needed
tree$tip.label <- gsub("_n[0-9]+$", "", tree$tip.label)

# Get one representative sample per haplotype for the tree
# (since tree tips are haplotypes, not individual samples)
haplo_metadata <- metadata_clean %>%
  group_by(Haplotype_clean) %>%
  slice_head(n = 1) %>%  # Take first sample for each haplotype
  ungroup()

# Create mapping between tree tips and metadata
tree_data <- data.frame(label = tree$tip.label) %>%
  left_join(haplo_metadata, by = c("label" = "Haplotype_clean"))

# Alternative approach: If you want each sample as a tip rather than haplotypes
# You might need to create a different tree or modify tip labels

# Create the circular phylogenetic tree
p <- ggtree(tree, layout = "circular", branch.length = "branch.length") %<+% tree_data +
  # Add tip points colored by lineage
  geom_tippoint(aes(color = Lineage), size = 1.5, alpha = 0.9) +
  
  # Add tip labels showing sample IDs 
  geom_tiplab(aes(label = City), 
              size = 1.8, 
              hjust = -0.1,
              offset = 0.005) +
  
  # Customize the appearance
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  # Custom color scale for your specific lineages
  scale_color_manual(
    name = "Lineage",
    values = c(
      "EP" = "#E8F4FD",      # East Plain (light blue)
      "SH" = "#00BFFF",      # Shanghai (cyan) 
      "SB" = "#FFB347",      # Sichuan Basin (orange)
      "LP" = "#90EE90",      # Loess Plateau (light green)
      "GP" = "#FFB6C1",      # Guizhou Plateau (light pink)
      "YP" = "#DDA0DD",      # Yunnan Plateau (plum)
      "unknown" = "grey70"   # Unknown samples
    ),
    na.value = "grey50",
    breaks = c("EP", "SH", "SB", "LP", "GP", "YP", "unknown")
  ) +
  
  # Add a title
  ggtitle("Spotted Lanternfly Mitochondrial Phylogeny by Lineage") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Create a clean circular phylogenetic tree with vivid colors
p_vivid <- ggtree(tree, layout = "circular", branch.length = "none") %<+% tree_data +
  # Add tip points colored by lineage - uniform size
  geom_tippoint(aes(color = Lineage), size = 1.5, alpha = 0.9) +
  
  # Add bigger tip labels showing sample IDs 
  geom_tiplab(aes(label = City), 
              size = 2.5, 
              hjust = -0.05,
              offset = 0.01) +
  
  # Minimal clean theme - no scales, backgrounds, etc.
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  
  # Vivid color palette as requested
  scale_color_manual(
    name = "Lineage",
    values = c(
      "EP" = "#e6ab02",      # Yellow-orange
      "SH" = "#a65628",      # Brown 
      "SB" = "#08589e",      # Blue
      "LP" = "#c861a5",      # Pink
      "GP" = "#ff7f00",      # Orange
      "YP" = "#6a3d9a",      # Purple
      "unknown" = "#4daf4a"  # Green
    ),
    na.value = "grey60"
  )

# Print the vivid plot
print(p_vivid)

# Alternative version: Color tip labels by country (like your code example)
# Country color palette
area_colors <- c(
  "JP" = "#e6ab02",
  "SK" = "#009E73", 
  "USA" = "#e41a1c",
  "CN" = "#6a3d9a"
)

# Create version with tip labels colored by country
p_country_labels <- ggtree(tree, layout = "circular", branch.length = "none") %<+% tree_data +
  # Points colored by lineage
  geom_tippoint(aes(color = Lineage), size = 1.5, alpha = 0.9) +
  
  # Tip labels colored by country
  geom_tiplab(aes(label = City, color = Country), 
              size = 2.8, 
              hjust = -0.1,
              offset = 0.02) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # First scale for lineage (points)
  scale_color_manual(
    name = "Lineage",
    values = c(
      "EP" = "#e6ab02", "SH" = "#a65628", "SB" = "#08589e",
      "LP" = "#c861a5", "GP" = "#ff7f00", "YP" = "#6a3d9a",
      "unknown" = "gray50"
    ),
    guide = guide_legend(order = 1)
  ) +
  
  # Second scale for country (labels) - this requires ggnewscale
  ggnewscale::new_scale_color() +
  geom_tiplab(aes(label = City, color = Country), 
              size = 2.8, 
              hjust = -0.1,
              offset = 0.02) +
  scale_color_manual(
    name = "Country",
    values = area_colors,
    guide = guide_legend(order = 2)
  )

print(p_country_labels)

# Save plots
ggsave("SLF_hap_circular_phylogeny_city.pdf", plot = p_vivid, 
       width = 12, height = 12, units = "in", dpi = 300)

ggsave("SLF_hap_circular_phylogeny_city_labels.pdf", plot = p_country_labels, 
       width = 12, height = 12, units = "in", dpi = 300)


##### Tip Label by sample ID
# Make sure SampleID and Lineage are clean

library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)
library(ggnewscale)

# Set working directory and read tree and metadata
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/MT/Final/ALL510")
tree <- read.tree("WG_with_outgroup.contree")
tree <- root(tree, outgroup = "Lycorma_meliae_outgroup", resolve.root = TRUE)

# Load metadata
metadata <- read.csv("Final_MT_meta.csv", header = TRUE) %>%
  mutate(across(c(SampleID, Haplotype, Lineage, Country, City), as.character))

# Step 1: Extract tree tip labels (they are haplotype names like "CXLVII")
tips <- tree$tip.label

# Step 2: For each tip, match a single representative sample
haplo_to_sample <- metadata %>%
  filter(Haplotype %in% tips) %>%
  group_by(Haplotype) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Haplotype, SampleID, Country, Lineage, City)

# Step 3: Build tree_data with label = tip (haplotype) and sample info
tree_data <- data.frame(label = tips) %>%
  left_join(haplo_to_sample, by = c("label" = "Haplotype"))

# Sanity check
stopifnot(!any(is.na(tree_data$SampleID)))

# Color palettes
lineage_colors <- c(
  "EP" = "#e6ab02", "SH" = "#a65628", "SB" = "#08589e",
  "LP" = "#c861a5", "GP" = "#ff7f00", "YP" = "#6a3d9a"
)
country_colors <- c(
  "JP" = "#e6ab02", "SK" = "#009E73", 
  "USA" = "#e41a1c", "CN" = "#6a3d9a",
  "Lycorma_meliae_outgroup" = "green"
)

# Plot: Show representative SampleID instead of haplotype name
p <- ggtree(tree, layout = "circular", branch.length = "none") %<+% tree_data +
  geom_tippoint(aes(color = Lineage), size = 1.5, alpha = 0.9) +
  scale_color_manual(name = "Lineage", values = lineage_colors, na.value = "grey60") +
  ggnewscale::new_scale_color() +
  geom_tiplab(aes(label = SampleID, color = Country), 
              size = 2.5, hjust = -0.1, offset = 0.02) +
  scale_color_manual(name = "Country", values = country_colors, na.value = "grey60") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save output
ggsave("SLF_MT_phylogeny_tiplab_SampleID.pdf", plot = p, width = 12, height = 12, dpi = 300)

# Show in viewer
print(p)


## 3.3 SNP-based phylogenetic tree
#Tree built part was done on hpc scripts

# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pcadapt)
library(ape)
library(phytools)

## set seed
set.seed(9)

# Set working directory
setwd("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Pop Structure")

# Load population data
popdata_USA <- read.csv("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/USA_SLF_MS_Info.csv", header = TRUE)
popdata_USA$Location <- as.factor(popdata_USA$Location)

area_colors <- c(
  "Bronx" = "#252525", 
  "Brooklyn" = "#e6ab02",
  "Camden" = "#a65628",
  "Connecticut" = "#e41a1c",
  "Governor Island" = "#08589e",
  "Manhattan" = "#c861a5",
  "Staten Island" = "#3690c0", 
  "West Orange" = "#009E73"
)

# Read the tree file
treefile_USA <- read.tree("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Pop Structure/USA_Final.SNPTree.contree")

# Map each tip label to its area based on popdata
tip_areas_USA <- popdata_USA$Location[match(treefile_USA$tip.label, popdata_USA$SLF_ID)]

if(any(is.na(tip_areas_USA))) {
  warning("Some tree tip labels did not match the SLF_ID in your metadata!")
}

# Replace tip labels with area names
treefile_USA$tip.label <- as.character(tip_areas_USA)

# Match the areas to the colors defined in location_colors
tip_colors <- area_colors[as.character(tip_areas_USA)]

# Midpoint root the tree
#Midpoint rooting is a common and acceptable approach if you lack an obvious outgroup. It will help in visualizing the tree in a more balanced manner
midpoint_tree_USA <- midpoint.root(treefile_USA)


# -----------------------------
# Assuming all libraries are loaded and the previous code has been run,
# and that 'midpoint_tree_USA' and 'tip_colors' are defined.
# -----------------------------

# Define the nodes I want to collapse, basically Manhattan samples for a better visulization
nodes <- c(145, 144, 150, 164, 165, 170, 171, 173, 100, 191, 189, 183, 110, 118, 120, 124, 133)

# Create a working copy of the midpoint-rooted tree for collapsing clades
collapsed_tree <- midpoint_tree_USA

# Paint each subtree from the given nodes so that when plotted,
# those clades are assigned a state ("1") that will be rendered as transparent.
for(i in seq_along(nodes)) {
  collapsed_tree <- paintSubTree(collapsed_tree, node = nodes[i], state = "1", anc.state = "0")
}

# Open a PDF device to save the plot
pdf("Collapsed_Arc_midpoint_tree_plot_with_uncollapsed_tip_labels.pdf", width = 12, height = 10)

# Plot the tree using an arc layout.
# Use a custom state-to-color mapping:
#   State 0 (unpainted) → "black"
#   State 1 (painted)   → "transparent"
# We disable tip labels in the plot (ftype = "off") so we can add them later selectively.
plot(collapsed_tree, setNames(c("black", "transparent"), 0:1), type = "arc",
     arc_height = 1, fsize = 0.6, ftype = "off")

# Retrieve the plotting coordinates (for nodes and tips) from the last plot
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Loop over each specified node to overlay a polygon (using the convex hull)
# that visually “collapses” the clade descending from that node.
for(i in seq_along(nodes)) {
  current_node <- nodes[i]
  
  # Get all descendant nodes of the current node; keep only tip indices.
  tips <- getDescendants(collapsed_tree, current_node)
  tips <- tips[tips <= Ntip(collapsed_tree)]
  
  # Combine the current node with its descendant tip indices.
  dd <- c(current_node, tips)
  
  # Calculate the convex hull of these points: the indices (with respect to dd)
  # that define the polygonal boundary.
  hull_indices <- chull(pp$xx[dd], pp$yy[dd])
  
  # Overlay the polygon for the clade with a black border and fill it with the Manhattan color.
  polygon(pp$xx[dd][hull_indices], pp$yy[dd][hull_indices],
          col = "#c861a5", border = "black", lwd = 1)
}

# --- Now add tip labels for the uncollapsed samples ---

# Identify all tips that fall within any collapsed clade.
collapsed_tips <- integer(0)
for(n in nodes) {
  tips <- getDescendants(collapsed_tree, n)
  tips <- tips[tips <= Ntip(collapsed_tree)]
  collapsed_tips <- union(collapsed_tips, tips)
}

# The uncollapsed tips are those tip indices not found in any collapsed clade.
uncollapsed_tips <- setdiff(1:Ntip(collapsed_tree), collapsed_tips)

# Add text labels for these uncollapsed tips.
# The tip labels (which we previously set to the area names) and their colors are retrieved from the tree.
text(x = pp$xx[uncollapsed_tips], y = pp$yy[uncollapsed_tips],
     labels = collapsed_tree$tip.label[uncollapsed_tips],
     col = tip_colors[uncollapsed_tips], cex = 0.6, pos = 4)

# Close the PDF device
dev.off()

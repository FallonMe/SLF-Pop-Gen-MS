## SNP-based phylogenetic tree
Tree built part was done on hpc scripts, plot on USA samples first

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
treefile_USA <- read.tree("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Pop Structure/New_LDPrune_pcangsd.tree")

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

# Plot the midpoint rooted tree
pdf("SNP_midpoint_rooted_tree_plot.pdf", width = 12, height = 10)
plot.phylo(midpoint_tree_USA, type = "fan", tip.color = tip_colors, rotate.tree = 90)
dev.off()

pdf("ARC_midpoint_rooted_tree_plot.pdf", width = 12, height = 10)
plotTree(midpoint_tree_USA , arc_height = 1, fsize = 0.6, type = "arc", arc_height = 1, fsize = 0.6)
nodelabels(bg = "white", cex = 0.6)
dev.off()

# -----------------------------
# Assuming all libraries are loaded and the previous code has been run,
# and that 'midpoint_tree_USA' and 'tip_colors' are defined.
# -----------------------------

# Define the nodes you want to collapse
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


```

Do it for all samples
```{r}
# Load required packages
library(phytools)
library(RColorBrewer)

# Define a helper function to add transparency (if not already defined)
make.transparent <- function(col, alpha=0.4) {
  rgb_val <- col2rgb(col)/255
  rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha=alpha)
}

# -----------------------------
# 1. Read and root the tree
# -----------------------------
tree <- read.tree("/Users/fallonmong/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/Pop Structure/Final.SNPTree.contree")
tree <- midpoint.root(tree)

# First plot the tree with node labels so you can verify node numbers
plotTree(tree, type="arc", arc_height=1, fsize=0.6)
nodelabels(bg="white", cex=0.6)

# -----------------------------
# 2. Open a new plotting device for the clade highlighting
#    (Suppress tip labels to get clean coordinate retrieval)
# -----------------------------
plotTree(tree, type="arc", arc_height=1, ftype="off", lwd=1)

# -----------------------------
# 3. Define the cladebox functions
# -----------------------------
cladebox <- function(node, col="#0000FF40", part=0.5, ...) {
  # Determine the number of points along the arc (default 30)
  if(hasArg(npoints)) {
    npoints <- list(...)$npoints
  } else {
    npoints <- 30
  }
  # Retrieve plotting coordinates from the last plot
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  # Create a minimal phylo object from pp to use getDescendants()
  phy <- list(edge=pp$edge, Nnode=pp$Nnode, tip.label=1:pp$Ntip)
  class(phy) <- "phylo"
  
  # Get the descendant tips for the given node
  dd <- getDescendants(phy, node)
  dd <- dd[dd <= Ntip(phy)]
  
  # Compute the polar angles for all descendant tip positions
  theta <- atan(pp$yy[dd] / pp$xx[dd])
  theta[pp$xx[dd] < 0] <- pi + theta[pp$xx[dd] < 0]
  theta[pp$xx[dd] > 0 & pp$yy[dd] < 0] <- 2*pi + theta[pp$xx[dd] > 0 & pp$yy[dd] < 0]
  
  # Expand the angle range a little to form the box
  max_theta <- max(theta) + 0.25 * 2*pi/Ntip(phy) * part
  min_theta <- min(theta) - 0.25 * 2*pi/Ntip(phy) * part
  
  # Compute a buffer based on the distance from the origin
  max_h <- max(sqrt(pp$xx^2 + pp$yy^2))
  h1 <- sqrt(pp$xx[node]^2 + pp$yy[node]^2) - 0.02 * max_h
  h2 <- max(sqrt(pp$xx[dd]^2 + pp$yy[dd]^2)) + 0.02 * max_h
  
  # Compute a sequence of angles to form the arc boundary
  tt <- seq(min_theta, max_theta, length.out = npoints)
  out_arc.x <- h2 * cos(tt)
  out_arc.y <- h2 * sin(tt)
  in_arc.x <- h1 * cos(tt)
  in_arc.y <- h1 * sin(tt)
  
  # Combine the outer and inner arcs to form a closed polygon
  xx <- c(out_arc.x, in_arc.x[npoints:1])
  yy <- c(out_arc.y, in_arc.y[npoints:1])
  
  # Draw the polygon (no border by default)
  polygon(xx, yy, lwd=5, col=col, border=FALSE)
}

# Function for multiple nodes: unions the descendant tips of all specified nodes.
cladeboxMulti <- function(nodes, col="#0000FF40", part=0.5, npoints=30) {
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  phy <- list(edge=pp$edge, Nnode=pp$Nnode, tip.label=1:pp$Ntip)
  class(phy) <- "phylo"
  
  dd_all <- integer(0)
  for (n in nodes) {
    dd <- getDescendants(phy, n)
    dd <- dd[dd <= Ntip(phy)]
    dd_all <- union(dd_all, dd)
  }
  
  # Use the first node as a representative for coordinate buffering
  rep_node <- nodes[1]
  
  theta <- atan(pp$yy[dd_all] / pp$xx[dd_all])
  theta[pp$xx[dd_all] < 0] <- pi + theta[pp$xx[dd_all] < 0]
  theta[pp$xx[dd_all] > 0 & pp$yy[dd_all] < 0] <- 2*pi + theta[pp$xx[dd_all] > 0 & pp$yy[dd_all] < 0]
  
  max_theta <- max(theta) + 0.25 * 2*pi/Ntip(phy) * part
  min_theta <- min(theta) - 0.25 * 2*pi/Ntip(phy) * part
  max_h <- max(sqrt(pp$xx^2 + pp$yy^2))
  h1 <- sqrt(pp$xx[rep_node]^2 + pp$yy[rep_node]^2) - 0.02 * max_h
  h2 <- max(sqrt(pp$xx[dd_all]^2 + pp$yy[dd_all]^2)) + 0.02 * max_h
  
  tt <- seq(min_theta, max_theta, length.out = npoints)
  out_arc.x <- h2 * cos(tt)
  out_arc.y <- h2 * sin(tt)
  in_arc.x <- h1 * cos(tt)
  in_arc.y <- h1 * sin(tt)
  xx <- c(out_arc.x, in_arc.x[npoints:1])
  yy <- c(out_arc.y, in_arc.y[npoints:1])
  
  polygon(xx, yy, lwd=5, col=col, border=FALSE)
}

# -----------------------------
# 4. Highlight your groups
# -----------------------------
# Define colors for the groups (you can adjust as desired)
usa_col <- make.transparent("#377eb8", 0.5)         # USA samples: node 137
shang_urban_col <- make.transparent("#ff7f00", 0.5)      # Shanghai Urban: node 121
shang_rural1_col <- make.transparent("#6a3d9a", 0.5)      # Shanghai Rural 1: node 134
shang_rural2_col <- make.transparent("#4daf4a", 0.5)     # Shanghai Rural 2: nodes 133, 135, 136

# Group 1: USA samples
cladebox(node = 137, col = usa_col)
# Group 2: Shanghai Urban
cladebox(node = 121, col = shang_urban_col)
# Group 3: Shanghai Rural 1
cladebox(node = 134, col = shang_rural1_col)
# Group 4: Shanghai Rural 2 (combine nodes 138, 135, 136)
cladeboxMulti(nodes = c(133, 135, 136), col = shang_rural2_col)

# -----------------------------
# 5. Overlay the tree outline for clarity
# -----------------------------
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=4, color="white")
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=2, color="black")

# -----------------------------
# 6. Save the final plot to PDF and PNG
# -----------------------------
# Save as PDF
pdf("FinalTreePlot.pdf", width=10, height=8)
# Re-run the plotting commands in the PDF device:
plotTree(tree, type="arc", arc_height=1, ftype="off", lwd=1)
cladebox(node = 137, col = usa_col)
cladebox(node = 121, col = shang_urban_col)
cladebox(node = 134, col = shang_rural1_col)
cladeboxMulti(nodes = c(133, 135, 136), col = shang_rural2_col)
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=4, color="white")
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=2, color="black")
dev.off()

# Save as a high-quality PNG
png("FinalTreePlot.png", width=2000, height=1600, res=300)
plotTree(tree, type="arc", arc_height=1, ftype="off", lwd=1)
cladebox(node = 137, col = usa_col)
cladebox(node = 121, col = shang_urban_col)
cladebox(node = 134, col = shang_rural1_col)
cladeboxMulti(nodes = c(133, 135, 136), col = shang_rural2_col)
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=4, color="white")
plotTree(tree, type="arc", arc_height=1, ftype="off", add=TRUE, lwd=2, color="black")
dev.off()


#OK. That’s not bad. Now let’s make one for Twitter.

pdf("Twitter_TreePlot.pdf", width = 10, height = 8)

# Set up a Twitter-style plot (white foreground, black background)
par(fg = "white", bg = "black")

# Plot the tree in "arc" style with a transparent branch color to capture coordinates.
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", lwd = 1, color = "transparent")

# Get a 4-color palette from Dark2
cols <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")

# Highlight the groups using your cladebox functions:
cladebox(node = 137, col = cols[1])                      # USA samples
cladebox(node = 121, col = cols[2])                      # Shanghai Urban
cladebox(node = 134, col = cols[3])                      # Shanghai Rural 1
cladeboxMulti(nodes = c(133, 135, 136), col = cols[4])   # Shanghai Rural 2

# Finally, overlay the tree outline for clarity
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", add = TRUE, lwd = 5, color = "white")
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", add = TRUE, lwd = 3, color = "black")

dev.off()

# ----- Save Twitter-style plot to high-resolution PNG -----
png("Twitter_TreePlot.png", width = 2000, height = 1600, res = 300)

# Set up a Twitter-style plot (white foreground, black background)
par(fg = "white", bg = "black")

# Plot the tree in "arc" style with a transparent branch color to capture coordinates.
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", lwd = 1, color = "transparent")

# Get the same 4-color palette from Dark2
cols <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")

# Highlight the groups using your cladebox functions:
cladebox(node = 137, col = cols[1])                      # USA samples
cladebox(node = 121, col = cols[2])                      # Shanghai Urban
cladebox(node = 134, col = cols[3])                      # Shanghai Rural 1
cladeboxMulti(nodes = c(133, 135, 136), col = cols[4])   # Shanghai Rural 2

# Finally, overlay the tree outline for clarity
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", add = TRUE, lwd = 5, color = "white")
plotTree(tree, type = "arc", arc_height = 0.5, ftype = "off", add = TRUE, lwd = 3, color = "black")

dev.off()

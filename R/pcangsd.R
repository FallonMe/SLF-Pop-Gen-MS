setwd("/scratch/work/winchell/Lanternfly/WGS/ANGSD/New_PCAngsd")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ape)
library(phytools)

area_colors <- c(
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

# Load population data
popdata <- read.csv("/scratch/work/winchell/Lanternfly/WGS/ANGSD/ANGSD_Final_SLF_MS_Info.csv", header = TRUE)  # Population metadata

popbed <- popdata[,c("SLF_ID","Location")]
poplist.names <- popbed[,2]

# Load covariance matrix
cov <- as.matrix(read.table("New_LDPrune_pcangsd.cov", header = FALSE))

# Perform PCA
e <- eigen(cov)

# Extract eigenvectors and eigenvalues
eigenvectors <- as.data.frame(e$vectors)
colnames(eigenvectors) <- paste0("PC", 1:ncol(eigenvectors))
eigenvectors <- cbind(popbed, eigenvectors)

# Proportion of variance explained
explained_variance <- e$values / sum(e$values) * 100
explained_variance

write.csv(explained_variance, "explained_variance_PCAngsd.csv")


# Subset for Shanghai (samples 1–20)
cov_SH <- cov[1:20, 1:20]
popbed_SH <- popbed[1:20, ]

# Subset for USA (samples 21–118)
cov_USA <- cov[21:118, 21:118]
popbed_USA <- popbed[21:118, ]

# Shanghai PCA
e_SH <- eigen(cov_SH)
evec_SH <- as.data.frame(e_SH$vectors)
colnames(evec_SH) <- paste0("PC", 1:ncol(evec_SH))
evec_SH <- cbind(popbed_SH, evec_SH)
var_SH <- e_SH$values / sum(e_SH$values) * 100

# USA PCA
e_USA <- eigen(cov_USA)
evec_USA <- as.data.frame(e_USA$vectors)
colnames(evec_USA) <- paste0("PC", 1:ncol(evec_USA))
evec_USA <- cbind(popbed_USA, evec_USA)
var_USA <- e_USA$values / sum(e_USA$values) * 100

library(ggplot2)

pca_SH_ellipse <- ggplot(evec_SH, aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Location), type = "norm", linetype = "solid", size = 0.7) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd - Shanghai",
    x = paste0("PC1 (", round(var_SH[1], 2), "%)"),
    y = paste0("PC2 (", round(var_SH[2], 2), "%)"),
    color = "Location"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCAngsd_SH_PC12_ellipse.pdf", pca_SH_ellipse, width = 6, height = 5)
ggsave("PCAngsd_SH_PC12.pdf", pca_SH_ellipse, width = 6, height = 6)


pca_USA_ellipse <- ggplot(evec_USA, aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Location), type = "norm", linetype = "solid", size = 0.7) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd - USA",
    x = paste0("PC1 (", round(var_USA[1], 2), "%)"),
    y = paste0("PC2 (", round(var_USA[2], 2), "%)"),
    color = "Location"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCAngsd_USA_PC12_ellipse.pdf", pca_USA_ellipse, width = 6, height = 5)
ggsave("PCAngsd_USA_PC12.pdf", pca_USA_ellipse, width = 6, height = 6)

# Simple PCA plot
PCAngsd_pca12 <- ggplot(eigenvectors, aes(x = PC1, y = PC2, color = Location, label = SLF_ID)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd",
    x = paste0("PC1 (", round(explained_variance[1], 2), "% variance)"),
    y = paste0("PC2 (", round(explained_variance[2], 2), "% variance)"),
    color = "Location"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12)
  )

ggsave("new_PCAngsd_PCA12.pdf", plot = PCAngsd_pca12 )


# Simple PCA plot
PCAngsd_pca13 <- ggplot(eigenvectors, aes(x = PC1, y = PC3, color = Location, label = SLF_ID)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd",
    x = paste0("PC1 (", round(explained_variance[1], 2), "% variance)"),
    y = paste0("PC3 (", round(explained_variance[3], 2), "% variance)"),
    color = "Location"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12)
  )

ggsave("New_PCAngsd_PCA13.pdf", plot = PCAngsd_pca13 )


# Simple PCA plot
PCAngsd_pca23 <- ggplot(eigenvectors, aes(x = PC2, y = PC3, color = Location, label = SLF_ID)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd",
    x = paste0("PC2 (", round(explained_variance[2], 2), "% variance)"),
    y = paste0("PC3 (", round(explained_variance[3], 2), "% variance)"),
    color = "Location"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12)
  )

ggsave("PCAngsd_PCA23.pdf", plot = PCAngsd_pca23 )



# PCAngsd construct neighbour-joining tree of samples from estimated covariance matrix estimated based on indivdual allele frequencies.
library(ape)
library(tidyverse)
library(phytools)

tree_file <- "/scratch/work/winchell/Lanternfly/WGS/ANGSD/New_PCAngsd/New_LDPrune_pcangsd.tree"
treefile <- read.tree(tree_file)

# plot(tree, main = "Neighbor-Joining Tree")

# Check for unmatched tip labels
unmatched_labels <- setdiff(treefile$tip.label, popdata$SLF_ID)
if (length(unmatched_labels) > 0) {
  stop(paste("Error: The following tip labels could not be matched to popdata$SLF_ID:", 
             paste(unmatched_labels, collapse = ", ")))
}

# Map each tip label to its area based on popdata
tip_areas <- popdata$Location[match(treefile$tip.label, popdata$SLF_ID)]

# Check for unmatched areas
if (any(is.na(tip_areas))) {
  stop("Error: Some tree tip labels could not be matched to popdata$Location. Check your data.")
}

# Replace tip labels with area names
treefile$tip.label <- as.character(tip_areas)

# Match the areas to the colors defined in location_colors
missing_areas <- setdiff(unique(tip_areas), names(area_colors))
if (length(missing_areas) > 0) {
  stop(paste("Error: Missing color definitions for the following areas:", paste(missing_areas, collapse = ", ")))
}
tip_colors <- area_colors[as.character(tip_areas)]

# Midpoint root the tree
midpoint_tree <- midpoint.root(treefile)

# Plot the midpoint rooted tree
pdf("New_PCAngsd_midpoint_rooted_tree_plot.pdf", width = 10, height = 8)
plot.phylo(midpoint_tree, type = "f", tip.color = tip_colors)
dev.off()

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
tree <- midpoint.root(treefile)

# Save the first plot with node labels
pdf("tree_with_node_labels.pdf", width = 8, height = 6)
plotTree(tree, type = "arc", arc_height = 1, fsize = 0.6, tip.color = tip_colors)
nodelabels(bg = "white", cex = 0.5)
dev.off()

# Save the second plot for clade highlighting (suppressing tip labels)
pdf("tree_clade_highlighting.pdf", width = 8, height = 6)
plotTree(tree, type = "arc", arc_height = 1, ftype = "off", lwd = 1)
dev.off()

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

# Admix, PCAngsd suggest k=4


# Subset for Shanghai (samples 1–20)
cov_SH <- cov[1:20, 1:20]
popbed_SH <- popbed[1:20, ]

# Subset for USA (samples 21–118)
cov_USA <- cov[21:118, 21:118]
popbed_USA <- popbed[21:118, ]

# Shanghai PCA
e_SH <- eigen(cov_SH)
evec_SH <- as.data.frame(e_SH$vectors)
colnames(evec_SH) <- paste0("PC", 1:ncol(evec_SH))
evec_SH <- cbind(popbed_SH, evec_SH)
var_SH <- e_SH$values / sum(e_SH$values) * 100

# USA PCA
e_USA <- eigen(cov_USA)
evec_USA <- as.data.frame(e_USA$vectors)
colnames(evec_USA) <- paste0("PC", 1:ncol(evec_USA))
evec_USA <- cbind(popbed_USA, evec_USA)
var_USA <- e_USA$values / sum(e_USA$values) * 100

library(ggplot2)

pca_SH_ellipse <- ggplot(evec_SH, aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Location), type = "norm", linetype = "solid", size = 0.7) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd - Shanghai",
    x = paste0("PC1 (", round(var_SH[1], 2), "%)"),
    y = paste0("PC2 (", round(var_SH[2], 2), "%)"),
    color = "Location"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCAngsd_SH_PC12_ellipse.pdf", pca_SH_ellipse, width = 6, height = 5)
ggsave("PCAngsd_SH_PC12.pdf", pca_SH_ellipse, width = 6, height = 6)


pca_USA_ellipse <- ggplot(evec_USA, aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Location), type = "norm", linetype = "solid", size = 0.7) +
  scale_color_manual(values = area_colors) +
  labs(
    title = "PCAngsd - USA",
    x = paste0("PC1 (", round(var_USA[1], 2), "%)"),
    y = paste0("PC2 (", round(var_USA[2], 2), "%)"),
    color = "Location"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCAngsd_USA_PC12_ellipse.pdf", pca_USA_ellipse, width = 6, height = 5)
ggsave("PCAngsd_USA_PC12.pdf", pca_USA_ellipse, width = 6, height = 6)


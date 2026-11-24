# Load required libraries
library(adegenet)
library(vcfR)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)

# Set seed for reproducibility
set.seed(0)
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/DAPC.RData")

# Read the VCF file
vcf <- read.vcfR("/scratch/work/winchell/Lanternfly/WGS/13_Filtered_VCF/Shanghai.vcf")
vcf_all <- read.vcfR("/scratch/work/winchell/Lanternfly/WGS/13_Filtered_VCF/Shanghai.vcf")

popdata_SH <- read.csv("/scratch/work/winchell/Lanternfly/WGS/13_Filtered_VCF/VCF_Stats/SH_SLF_MS_Info.csv", header=T) 
#link populations

popgen_sort <- popdata_SH[order(popdata_SH$SLF_ID),c("SLF_ID","Location")] 
#ordered and reduced to just columns we need

# Convert VCF to genlight object
slf.gen = vcfR2genlight(vcf)

#assigns the 'area' column of the sorted dataframe to the 'pop' slot of the genlight object.
pop(slf.gen) <- popgen_sort[,2] #assign populations to gen file
slf.gen@pop #verify population attribute

toRemove <- is.na(glMean(slf.gen, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
slf.gen.new <- slf.gen[, !toRemove]
glPca(slf.gen.new) # this should work

#find cluster
clusters_result <- find.clusters(slf.gen.new, parallel = require("parallel"), n.cores=10)

# cross-validation for the DAPC 
optx_slf_SH <- xvalDapc(tab(slf.gen.new, NA.method = "mean"), pop(slf.gen.new), parallel="multicore", ncpus=8L)  

#look for leveling off window
optx_slf_SH[2:7]

dapc.SH <- dapc(slf.gen.new,  n.pca=2, n.da=2, var.contrib=T, scale=F,parallel = require("parallel"), n.cores=8)

# Define colors for each population
location_colors <- c(
  "Jinqiao" = "#ff7f00",
  "Sheshan" = "#6a3d9a",
  "Tianma" = "#4daf4a"
)

#visualize DAPC
# Set up a PDF device
pdf("SH_DAPC.pdf")
# Plot the DAPC scatter plot
scatter(dapc.SH, clabel=F, legend=T, col=location_colors, cex=1.5, mstree=F) #retain 2 df 
# Turn off the device
dev.off() # Figure 3A

# Save to a high-quality, editable vector format (change file name as needed)
pdf("Final_DAPC_SH.pdf", width=7, height=7)  # 180 mm ≈ 7.1 inches

# Plot DAPC with improved visual settings
scatter(
  dapc.SH,
  clabel = 0,                      # Turn off default labels (use legend instead)
  legend = TRUE,
  col = location_colors,          # Preset vector of colors
  cex = 2,                       # Point size (consider cex=1.2–1.5 for print)
  mstree = FALSE,
  scree.da = FALSE,               # Optional: hide screeplot
  #cstar = 0,                      # No lines connecting clusters
  solid = 0.7                     # Transparency for better overlap visualization
)

# Finalize the figure
dev.off()


# Load required libraries
library(ggplot2)
library(reshape2)

# Extract population information
pops <- pop(slf.gen.new)  # Access population information for your updated data

# Create ordered individual names
orderedindivid <- rep(NA, nInd(slf.gen.new))  # Use nInd() to dynamically get the number of individuals
for (i in 1:nInd(slf.gen.new)) {
  orderedindivid[i] <- paste(pops[i], i, sep = "")
}

# Reorder individuals within each population
rownames(dapc.SH$posterior) <- orderedindivid

# Convert posterior probabilities to a dataframe
dapc.results <- as.data.frame(dapc.SH$posterior)
dapc.results$Original_Pop <- pops  # Add original population information
dapc.results$indNames <- rownames(dapc.results)  # Add individual names
dapc.results <- reshape2::melt(dapc.results)  # Convert to long format for ggplot
colnames(dapc.results) <- c("Original_Pop", "Sample", "Assigned_Pop", "Posterior_membership_probability")

# Define colors for populations
plotcolors1 <- c(
  "Jinqiao" = "#ff7f00",
  "Sheshan" = "#6a3d9a",
  "Tianma" = "#4daf4a"
)

# Create the posterior probability plot
posterior_probability_plot <- ggplot(dapc.results, aes(x = Sample, y = Posterior_membership_probability, fill = Assigned_Pop)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = plotcolors1) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  ylab("Posterior Membership Probability") +
  xlab("Individuals")

# Save the plot as PNG and PDF
ggsave(filename = "SH_Posterior_Probability_Plot.png",
       plot = posterior_probability_plot,
       width = 14, height = 5, dpi = 300)

ggsave(filename = "SH_Posterior_Probability_Plot.pdf",
       plot = posterior_probability_plot,
       width = 14, height = 5)


# Define colors for each population
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

# Save to a high-quality, editable vector format (change file name as needed)
pdf("Wide_Final_USA_DAPC.pdf", width=8, height=6)  # 180 mm ≈ 7.1 inches

# Plot DAPC with improved visual settings
scatter(
  dapc.USA,
  clabel = 0,                      # Turn off default labels (use legend instead)
  legend = TRUE,
  col = location_colors,          # Preset vector of colors
  cex = 2,                       # Point size (consider cex=1.2–1.5 for print)
  mstree = FALSE,
  scree.da = FALSE,               # Optional: hide screeplot
  #cstar = 0,                      # No lines connecting clusters
  solid = 0.7                     # Transparency for better overlap visualization
)

# Finalize the figure
dev.off()

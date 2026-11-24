# SLF Mitochondrial Haplotype Analysis and Visualization

# Load necessary libraries
library(ape)
library(pegas)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vcfR)
library(factoextra)
library(RColorBrewer)

# Set working directory
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/MT/Final")

# Load aligned mitochondrial sequences
dna <- read.dna("WG_511_mafftstrict_aligned.fasta", format = "fasta")

# Generate haplotypes
haps <- haplotype(dna)

# Save haplotype FASTA
rownames(haps) <- paste0("Hap", seq_len(nrow(haps)))
write.dna(as.matrix(haps), "Final_WG_haplotypes.fasta", format = "fasta", nbcol = -1, colsep = "")

# Map samples to haplotypes
hap_index <- attr(haps, "index")
sample_names <- rownames(dna)
hap_map <- data.frame(
  SampleID = unlist(lapply(hap_index, function(i) sample_names[i])),
  Haplotype = rep(rownames(haps), lengths(hap_index))
)
write.csv(hap_map, "Final_haplotype_membership.csv", row.names = FALSE)

# Filter for shared haplotypes (≥2 samples)
hap_counts <- table(hap_map$Haplotype)
shared_haps <- names(hap_counts[hap_counts >= 2])
shared_samples <- hap_map %>% filter(Haplotype %in% shared_haps) %>% arrange(Haplotype)
write.csv(shared_samples, "shared_haplotypes_samples.csv", row.names = FALSE)

# Create and plot haplotype network
net <- haploNet(haps)
plot(net, size = attr(net, "freq"), scale.ratio = 0.5, labels = TRUE, cex = 0.7)
setwd("~/Documents/Documents - Fallon MacBook Pro/SLF Pop Gen 2025/From_VCF/LD")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(nls2)

# Load PopLDdecay output files for Shanghai and USA populations
sh_ld <- read.table("SH.stat", header=FALSE)
usa_ld <- read.table("USA.stat", header=FALSE)

# Add proper column names
colnames(sh_ld) <- c("dist_bin", "r2_mean", "V3", "V4", "V5", "n_pairs")
colnames(usa_ld) <- c("dist_bin", "r2_mean", "V3", "V4", "V5", "n_pairs")

# Add population identifier
sh_ld$population <- "Shanghai"
usa_ld$population <- "USA"

# Combine data
all_ld <- rbind(sh_ld, usa_ld)

# Convert distance bin to actual distance in kb
# The distance bin (V1) appears to be in kb already, but let's confirm this
# by examining the range
range(sh_ld$dist_bin)
range(usa_ld$dist_bin)

# Visualize LD decay curves
p1 <- ggplot(all_ld, aes(x=dist_bin, y=r2_mean, color=population)) +
  geom_line() +
  labs(x="Distance (kb)", y="Mean r²", 
       title="LD Decay Curves for Shanghai and USA Populations") +
  theme_bw()

print(p1)

# Fit nonlinear model to estimate recombination rate
# Function to fit LD decay model
estimate_recomb_rate <- function(ld_data, pop_name) {
  # Convert distance bin to base pairs for the model
  ld_data$dist_bp <- ld_data$dist_bin * 1000
  
  # Fit Hill function model for LD decay
  # This is based on the expected theoretical decay of LD with distance
  # The model is r² ≈ 1/(1 + 4Nec*distance)
  # We'll use a slightly more complex version that accounts for initial r² and asymptotic r²
  
  # Define a grid of starting values
  # Refined grid with biologically meaningful increments from 0.1 to 30 cM/Mb
  start_grid <- expand.grid(
    d_0 = c(0.8, 0.9, 1.0),  # Initial r² at distance 0
    b = c(0.1, 0.2, 0.3),    # Asymptotic r²
    c = c(0.000001, 0.000005, 0.00001, 0.00003, 0.00005, 0.0001, 0.0002, 0.0003)  # Corresponding to 0.1, 0.5, 1, 3, 5, 10, 20, 30 cM/Mb
  )
  # Find best model using grid search
  model <- nls2(r2_mean ~ b + (d_0 - b)*exp(-c*dist_bp),
                data=ld_data,
                start=start_grid,
                algorithm="brute-force")
  
  # Extract parameters
  params <- coef(model)
  c_param <- params["c"]
  
  # Convert c to cM/Mb (100 cM = 1 Morgan)
  # This is a simplification - c is related to 4Nec, but we're extracting
  # just the c component as an approximation
  cM_per_Mb <- c_param * 100000  # Convert to cM/Mb
  
  # Return results
  return(list(
    population = pop_name,
    model = model,
    c_param = c_param,
    cM_per_Mb = cM_per_Mb
  ))
}

# Run for each population
sh_results <- estimate_recomb_rate(sh_ld, "Shanghai")
usa_results <- estimate_recomb_rate(usa_ld, "USA")

# Extract and compare recombination rates
results_df <- data.frame(
  Population = c("Shanghai", "USA"),
  c_param = c(sh_results$c_param, usa_results$c_param),
  cM_per_Mb = c(sh_results$cM_per_Mb, usa_results$cM_per_Mb)
)

print(results_df)

# Create predictions from fitted models for plotting
sh_pred <- data.frame(
  dist_bin = seq(min(sh_ld$dist_bin), max(sh_ld$dist_bin), length.out=100),
  population = "Shanghai"
)
sh_pred$dist_bp <- sh_pred$dist_bin * 1000
sh_pred$r2_pred <- predict(sh_results$model, newdata=sh_pred)

usa_pred <- data.frame(
  dist_bin = seq(min(usa_ld$dist_bin), max(usa_ld$dist_bin), length.out=100),
  population = "USA"
)
usa_pred$dist_bp <- usa_pred$dist_bin * 1000
usa_pred$r2_pred <- predict(usa_results$model, newdata=usa_pred)

# Combine predictions
all_pred <- rbind(sh_pred, usa_pred)

# Enhanced publication-quality plot with customized theme - font-safe version
p_publication <- ggplot() +
  geom_point(data=all_ld, aes(x=dist_bin, y=r2_mean, color=population), 
             alpha=0.4, size=1) +
  geom_line(data=all_pred, aes(x=dist_bin, y=r2_pred, color=population), 
            size=1) +
  scale_color_manual(values=c("Shanghai"="#c861a5", "USA"="#3690c0"),
                     name="Population") +
  labs(x="Distance (kb)", 
       y=expression(paste("Linkage Disequilibrium (", r^2, ")")), 
       title="LD Decay in Native and Invasive SLF Populations") +
  theme_minimal() +
  theme(
    # Remove font family specification to avoid errors
    text = element_text(size=12),  # Removed 'family="Arial"'
    plot.title = element_text(size=14, face="bold"),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill="white", color=NA),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", fill=NA)
  ) +
  # Add recombination rate annotations with the correct variable names
  annotate("text", x=max(all_ld$dist_bin)*0.7, y=max(all_ld$r2_mean)*0.3,
           label=paste("Shanghai:", round(sh_results$cM_per_Mb, 2), "cM/Mb"),
           color="#c861a5", hjust=0, size=4) +
  annotate("text", x=max(all_ld$dist_bin)*0.7, y=max(all_ld$r2_mean)*0.2,
           label=paste("USA:", round(usa_results$cM_per_Mb, 2), "cM/Mb"),
           color="#3690c0", hjust=0, size=4)

# Save the publication-quality plot
pdf("SLF_LD_decay_publication.pdf", width=8, height=6, useDingbats=FALSE)
print(p_publication)
dev.off()

# Also save as high-resolution PNG for easy viewing
png("SLF_LD_decay_publication.png", width=8, height=6, units="in", res=300)
print(p_publication)
dev.off()
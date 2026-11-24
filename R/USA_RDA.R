library(readr)   
library(data.table)
library(dplyr)
library(ggplot2)
library(vegan) # Used to run RDA
library(factoextra)
library(corrplot)
library(patchwork)
library(lme4)
library(psych)

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


#GEA analysis using all environmental data including WorldClim, data
#Population structure (K3) estimated using snmf, tessR and PCA
#Geographic distances calculated as MEM value from the coordinates. # MEM 8, 44, 54, 60, 74 are significant.

#####Moran's eigenvector maps (MEMs)######
#MEMs are used to model the spatial auto correlations in RDA
library(adespatial)
library(geosphere)
library(spdep)

setwd("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/")
#Read the Enviornment data
env_usa <- read_csv("USA_Final_SLF_MS_Info.csv")
str(env_usa)

#bio1-11 are temperature related env variables, bio12-19 are precipitation related env variables, others such as ALAN(artificial light at night), AIS(artificial impervious surface), TCC(tree canopy cover) are urbanization related variables. 
# so here we separate this three category
# Manually selecting only temperature-related columns (bio_1 to bio_11)
temp_vars <- env_usa %>%
  select(SLF_ID, Country, City, Location, detailed_area,lat, lon, 
         bio_1, bio_2, bio_3, bio_4, bio_5, bio_6, bio_7, bio_8, bio_9, bio_10, bio_11)

# Manually selecting only precipitation-related columns (bio_12 to bio_19)
precip_vars <- env_usa %>%
  select(SLF_ID, Country, City, Location, detailed_area, lat, lon, 
         bio_12, bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19)

# Manually selecting urbanization-related variables
urban_vars <- env_usa %>%
  select(SLF_ID, Country, City, Location,detailed_area, lat, lon, ALAN, AIS, TCC, Pop_Den)

# Save to CSV
write_csv(temp_vars, "USA_Temperature_Variables.csv")
write_csv(precip_vars, "USA_Precipitation_Variables.csv")
write_csv(urban_vars, "USA_Urbanization_Variables.csv")

#Explore the Env PCA

# Load required libraries
library(ggplot2)
library(dplyr)
library(factoextra)
library(tidyr)
library(readr)
library(psych)
library(reshape)

# Define colors for each location
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

# PCA for Temperature Variables
tempPCA_USA <- prcomp(temp_vars[,8:18], scale=T)
precipPCA_USA <- prcomp(precip_vars[,8:15], scale=T)
urbanPCA_USA <- prcomp(urban_vars[,8:11], scale=T)

# Summarize PCAs
summary(tempPCA_USA) 
summary(precipPCA_USA) 
summary(urbanPCA_USA) 

# Visualize Temperature PCA and save the result
tempPCA_USA$rotation[,1:4]

pdf("USA_fviz_eig_Temp_PCA.pdf", width = 8, height = 8)
fviz_eig(tempPCA_USA, addlabels = TRUE, ylim = c(0, 100))
dev.off()

# Visualize Precipitation PCA and save the result
precipPCA_USA$rotation[,1:4]

pdf("USA_fviz_eig_precip_PCA.pdf", width = 8, height = 8)
fviz_eig(precipPCA_USA, addlabels = TRUE, ylim = c(0, 100))
dev.off()

# Visualize Urbanization PCA and save the result
urbanPCA_USA$rotation[,1:3]

pdf("USA_fviz_eig_urban_PCA.pdf", width = 8, height = 8)
fviz_eig(urbanPCA_USA, addlabels = TRUE, ylim = c(0, 100))
dev.off()

# Create correlation plot for Temperature PCA
M_temp <- tempPCA_USA$rotation[,1:4]
melted_cormat_temp <- reshape::melt(M_temp)

level.order_temp <- rev(rownames(M_temp))

corplot_temp <- ggplot(melted_cormat_temp, aes(x=X2, y=factor(X1, level = level.order_temp),
                                               fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(min(melted_cormat_temp$value),max(melted_cormat_temp$value)), space = "Lab",
                       name="Value") + theme_minimal()+ xlab(label="") + ylab(label="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank())+
  coord_fixed(ratio=.25) + geom_text(aes(X2, X1, label = round(value,2)),
                                     color = "black", size = 3) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.hjust = 0.5)) 

ggsave('USA_corplot_temp.pdf', plot = corplot_temp, device = "pdf")

# Convert to dataframe format for easy reading
correlation_df_temp <- as.data.frame(M_temp)

# Save to CSV file
write.csv(correlation_df_temp, "USA_Temperature_PCA_Correlations.csv", row.names = TRUE)

# Correlation plot for Precipitation PCA
M_precip <- precipPCA_USA$rotation[, 1:4]
M_precip_correlation_df <- as.data.frame(M_precip)

# Save to CSV file
write.csv(M_precip_correlation_df, "USA_Precip_PCA_Correlations.csv", row.names = TRUE)

melted_cormat_precip <- reshape::melt(M_precip)

level.order_precip <- rev(rownames(M_precip))

corplot_precip <- ggplot(melted_cormat_precip, aes(x = X2, y = factor(X1, levels = level.order_precip),
                                                   fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(min(melted_cormat_precip$value), max(melted_cormat_precip$value)), space = "Lab",
                       name = "Value") + theme_minimal() + xlab(label = "") + ylab(label = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank()) +
  coord_fixed(ratio = .25) + geom_text(aes(X2, X1, label = round(value, 2)),
                                       color = "black", size = 3) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.hjust = 0.5)) 

ggsave('USA_corplot_precip.pdf', plot = corplot_precip, device = "pdf")

# Correlation plot for Urbanization PCA
M_urban <- urbanPCA_USA$rotation[, 1:3]  # Urbanization has fewer variables

M_urban_correlation_df <- as.data.frame(M_urban)

# Save to CSV file
write.csv(M_urban_correlation_df, "USA_Urban_PCA_Correlations.csv", row.names = TRUE)

melted_cormat_urban <- reshape::melt(M_urban)

level.order_urban <- rev(rownames(M_urban))

corplot_urban <- ggplot(melted_cormat_urban, aes(x = X2, y = factor(X1, levels = level.order_urban),
                                                 fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(min(melted_cormat_urban$value), max(melted_cormat_urban$value)), space = "Lab",
                       name = "Value") + theme_minimal() + xlab(label = "") + ylab(label = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank()) +
  coord_fixed(ratio = .25) + geom_text(aes(X2, X1, label = round(value, 2)),
                                       color = "black", size = 3) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.hjust = 0.5)) 

ggsave('USA_corplot_urban.pdf', plot = corplot_urban, device = "pdf")

# Get a correlation plot for all env factors
pdf("USA_correlation_plot_forALL.pdf", width=20, height=8)
corPlot(env_usa[,12:33], scale = TRUE)
dev.off()

# Do PCA on all env factors
envPCA_USA <- prcomp(env_usa[,12:33], scale=T) 
summary(envPCA_USA)
eigenvals(envPCA_USA)  # eigenvalues

# Visualize PCA and save the result
envPCA_USA$rotation[,1:4]

pdf("USA_biplot_envPCA.pdf", width = 20, height = 20)
biplot(envPCA_USA)
dev.off() 

pdf("USA_fviz_eig_envPCA.pdf", width = 7, height = 7)
fviz_eig(envPCA_USA, addlabels = TRUE, ylim = c(0, 100))
dev.off()

# Keep first 4 PCs and save the correlation plot
M_env <- envPCA_USA$rotation[,1:4]
melted_cormat_env <- reshape::melt(M_env)

level.order_env <- rev(rownames(M_env))

corplotv2 <- ggplot(melted_cormat_env, aes(x=X2, y=factor(X1, level = level.order_env),
                                           fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(min(melted_cormat_env$value),max(melted_cormat_env$value)), space = "Lab",
                       name="Value") + theme_minimal()+ xlab(label="") + ylab(label="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank())+
  coord_fixed(ratio=.25) + geom_text(aes(X2, X1, label = round(value,2)),
                                     color = "black", size = 3) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.hjust = 0.5)) 

ggsave('USA_corplot_env.pdf', plot = corplotv2, device = "pdf")
ggsave('USA_corplot_env.png', plot = corplotv2, device = "png", dpi = 300)

# Create dataframes with PCA scores
tempPCA_df_USA <- data.frame(temp_vars, tempPCA_USA$x) 
precipPCA_df_USA <- data.frame(precip_vars, precipPCA_USA$x) 
urbanPCA_df_USA <- data.frame(urban_vars, urbanPCA_USA$x) 

write_csv(tempPCA_df_USA, "USA_tempPCA_df.csv")
write_csv(precipPCA_df_USA, "USA_precipPCA_df.csv")
write_csv(urbanPCA_df_USA, "USA_urbanPCA_df.csv")

# Create scatter plots for Temperature PCA by location
tempPCA_area <- ggplot(tempPCA_df_USA, aes(x = PC1, y = PC2, color = Location)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(tempPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(tempPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = location_colors) +
  theme(legend.position = "right")

ggsave("USA_tempPCA_area_plot.pdf", tempPCA_area, width = 10, height = 8, units = "in")
ggsave("USA_tempPCA_area_plot.png", tempPCA_area, width = 10, height = 8, dpi = 300)

# Precipitation PCA plot
precipPCA_area <- ggplot(precipPCA_df_USA, aes(x = PC1, y = PC2, color = Location)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(precipPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(precipPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = location_colors) +
  theme(legend.position = "right")

ggsave("USA_precipPCA_area_plot.pdf", precipPCA_area, width = 10, height = 8, units = "in")
ggsave("USA_precipPCA_area_plot.png", precipPCA_area, width = 10, height = 8, dpi = 300)

# Urbanization PCA plot
urbanPCA_area <- ggplot(urbanPCA_df_USA, aes(x = PC1, y = PC2, color = Location)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(urbanPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(urbanPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = location_colors) +
  theme(legend.position = "right")

ggsave("USA_urbanPCA_area_plot.pdf", urbanPCA_area, width = 10, height = 8, units = "in")
ggsave("USA_urbanPCA_area_plot.png", urbanPCA_area, width = 10, height = 8, dpi = 300)

# Manhattan-specific plots
# Subset to plot only Manhattan areas
manhattan_temp_df <- subset(tempPCA_df_USA, detailed_area %in% c("Lower Manhattan", "Upper Manhattan", "Midtown Manhattan"))

# Define colors for the Manhattan plot
plotcolors <- c("Lower Manhattan" = "#c861a5", "Upper Manhattan" = "#009E73", "Midtown Manhattan" = "#3690c0")

# Temperature PCA plot for Manhattan areas
tempPCA_manhattan <- ggplot(manhattan_temp_df, aes(x = PC1, y = PC2, color = detailed_area)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(tempPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(tempPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = plotcolors) +
  theme(legend.position = "right")

ggsave("USA_tempPCA_manhattan_plot.pdf", tempPCA_manhattan, width = 10, height = 8, units = "in")
ggsave("USA_tempPCA_manhattan_plot.png", tempPCA_manhattan, width = 10, height = 8, dpi = 300)

# Subset Manhattan areas for precipitation and urbanization
manhattan_precip_df <- subset(precipPCA_df_USA, detailed_area %in% c("Lower Manhattan", "Upper Manhattan", "Midtown Manhattan"))
manhattan_urban_df <- subset(urbanPCA_df_USA, detailed_area %in% c("Lower Manhattan", "Upper Manhattan", "Midtown Manhattan"))

# Precipitation PCA plot for Manhattan areas
precipPCA_manhattan <- ggplot(manhattan_precip_df, aes(x = PC1, y = PC2, color = detailed_area)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(precipPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(precipPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = plotcolors) +
  theme(legend.position = "right")

ggsave("USA_precipPCA_manhattan_plot.pdf", precipPCA_manhattan, width = 10, height = 8, units = "in")
ggsave("USA_precipPCA_manhattan_plot.png", precipPCA_manhattan, width = 10, height = 8, dpi = 300)

# Urbanization PCA plot for Manhattan areas
urbanPCA_manhattan <- ggplot(manhattan_urban_df, aes(x = PC1, y = PC2, color = detailed_area)) +
  geom_jitter(size = 3, alpha = 0.8) + 
  theme_bw() +  # White background theme
  labs(x = paste0("PC1 (", round(summary(urbanPCA_USA)$importance[2,1]*100, 1), "%)"), 
       y = paste0("PC2 (", round(summary(urbanPCA_USA)$importance[2,2]*100, 1), "%)")) + 
  scale_color_manual(values = plotcolors) +
  theme(legend.position = "right")

ggsave("USA_urbanPCA_manhattan_plot.pdf", urbanPCA_manhattan, width = 10, height = 8, units = "in")
ggsave("USA_urbanPCA_manhattan_plot.png", urbanPCA_manhattan, width = 10, height = 8, dpi = 300)

# Create final CSV with PC1 values for selection analysis
# Extract PC1 scores from each PCA object
temp_PC1 <- tempPCA_USA$x[, 1]
precip_PC1 <- precipPCA_USA$x[, 1]
urban_PC1 <- urbanPCA_USA$x[, 1]

# Create a data frame with the PC1 scores and a unique identifier
pc_scores <- data.frame(
  SLF_ID = temp_vars$SLF_ID,
  temp_PC1 = temp_PC1,
  precip_PC1 = precip_PC1,
  urban_PC1 = urban_PC1
)

# Extract the key environmental metadata columns from the original env data
env_metadata <- env_usa %>% 
  select(SLF_ID, Country, City, Location, detailed_area, date_col, lat, lon, sex, color, habitat)

# Merge the metadata with the PC1 scores using SLF_ID as the key
env_with_PC1 <- left_join(env_metadata, pc_scores, by = "SLF_ID")

# Save the final merged data to a CSV file for later selection analyses
write_csv(env_with_PC1, "USA_Environmental_PC1s.csv")

load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")
## load geographic coordinates
geo = env_usa[,1:8]
head(geo)

geo <- unique(geo)

#Step 2: Compute Mean Geographic Coordinates per Deme
grp.coord <- cbind(
  tapply(geo$lon, INDEX = geo$lat, mean), 
  tapply(geo$lat, INDEX = geo$lat, mean)
)

#Step 3: Create the Spatial Network (Gabriel Graph)
#This defines neighborhood connectivity using a Gabriel graph.
nb <- graph2nb(gabrielneigh(grp.coord, nnmult = 100), sym=TRUE)
plot(nb, coords = grp.coord)

#Step 4: Compute the Spatial Weights Matrix
listW   <- nb2listw(nb, style="W")   # Convert neighborhood graph to weights list

#Step 5: Calculate Distance-Based Spatial Weights
disttri <- nbdists(nb, grp.coord, longlat = TRUE)  # Compute pairwise distances
fdist   <- lapply(disttri, function(x) x^(-1))  # Apply inverse distance weighting
listW   <- nb2listw(nb, glist=fdist, style="W")  # Apply revised weights

#Step 6: Compute Moran's Eigenvector Maps (MEMs)
MEM <- scores.listw(listW, MEM.autocor = "all") #Extracts spatial eigenvectors capturing spatial autocorrelation.
rownames(MEM) <- rownames(grp.coord)

#Step 7: Assign MEMs to Samples
mem.full <- as.data.frame(apply(MEM, 2, function(x) { x[match(geo$lat, rownames(MEM))] }))
rownames(mem.full) <- geo$SLF_ID


# Read the geno_data, converted using plink
slf_gen <- read.table("USA_final_annotated_raw.raw", header=T)
dim(slf_gen)
str(slf_gen)

slf.gen.matrix <- as.matrix(slf_gen)

gen <- slf.gen.matrix[, -c(2:6)] #extract SNP columns
colnames(gen)[1] <- 'SLF_ID'

sum(is.na(gen))

str(gen)

# Exclude the first column (SLF_ID) and convert the rest to numeric
gen_numeric <- apply(gen[, -1], 2, as.numeric)
str(gen_numeric)

#Both LFMM and RDA require complete data frames (i.e., no missing genetic data)
gen.imp <- apply(gen_numeric, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

sum(is.na(gen.imp)) # No NAs
# should be 0

#Step 8: Forward Selection of Significant MEMs
#gen.imp is the SLF genotype data 
fwsel.mem <- forward.sel(Y = as.matrix(gen.imp), X = as.matrix(mem.full), alpha = 0.05) # significant if p-value < 0.05
head(mem.full)

# extract the significant MEMs (p < 0.05)
sel.mem <- mem.full[,colnames(mem.full) %in% fwsel.mem$variables]
write.csv(sel.mem, "USA_selected_mem.csv", row.names = TRUE)

load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")

# "MEM43","MEM53","MEM59","MEM64" are significant.



############ Run RDA for USA samples only, no clear pop structure, so no need to do the pop structure corretion ####
# But beacuse we have a lot of samples are geographically close, so it's nice to do the spatial autocorrection #
# So let's run a simple RDA without controlling space, and one RDA with space correction (we will use the the significant MEMs we created above) #
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")

env_pcs_df_area <- read.csv("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_Environmental_PC1s.csv")

mem_selected <- read.csv("USA_selected_mem.csv")   # "MEM43","MEM53","MEM59","MEM64" are significant, which we saved above

# Extract relevant MEM variables
mem_selected <- mem_selected[, c("MEM43", "MEM53", "MEM59", "MEM64")]
rownames(mem_selected) <- env_pcs_df_area$SLF_ID  # Ensure matching IDs
#mem_selected <- mem_selected[ -c(1) ] 

# Extract genotype data (Remove metadata columns)
slf.gen.matrix <- as.matrix(slf_gen[, -c(1:6)])  
rownames(slf.gen.matrix) <- slf_gen$FID  # Set row names for alignment

# Ensure all datasets are aligned
env_pcs_df_area <- env_pcs_df_area[order(env_pcs_df_area$SLF_ID),]
slf.gen.matrix <- slf.gen.matrix[order(rownames(slf.gen.matrix)),]
mem_selected <- mem_selected[order(rownames(mem_selected)),]

# Merge environmental predictors (First PCs of temp, precipitation, urbanization)
env_pred <- env_pcs_df_area[, c("temp_PC1", "precip_PC1", "urban_PC1")]
rownames(env_pred) <- env_pcs_df_area$SLF_ID

# Merge all covariates into a single data frame
all_covariates <- cbind(env_pred, mem_selected)
rownames(all_covariates) <- env_pcs_df_area$SLF_ID

# Ensure no missing values
sum(is.na(all_covariates))  # Should return 0
write.csv(all_covariates, "USA_all_covariates.csv", row.names = TRUE)

save.image("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")

########## Run Below as script on hpc #############
## simple RDA without controlling space is done before on hpc###
## simple RDA###
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")

## simple RDA###
library(vegan)

# **Step 1: Run Simple RDA (Without Controlling for Population & Space)**
rda.full <- rda(gen.imp ~ temp_PC1 + precip_PC1 + urban_PC1, 
                data = all_covariates, scale = TRUE)

saveRDS(rda.full, "rda_full_model.rds")

# Overall significance of the RDA
signif.full <- anova.cca(rda.full, permutations = 5000, parallel=getOption("mc.cores"))
signif.full

# Compute Adjusted R²
rda_full_r2 <- RsquareAdj(rda.full)
write.csv(data.frame(R2 = rda_full_r2$r.squared, Adj_R2 = rda_full_r2$adj.r.squared), 
          "rda_full_r2.csv", row.names = FALSE)

# Compute Permutation Tests
rda_full_anova <- anova.cca(rda.full, permutations = 5000)
write.csv(rda_full_anova, "rda_full_anova.csv", row.names = FALSE)

rda_full_anova_axis <- anova.cca(rda.full, by = "axis", permutations = 5000)
write.csv(rda_full_anova_axis, "rda_full_anova_axis.csv", row.names = FALSE)

rda_full_anova_terms <- anova.cca(rda.full, by = "terms", permutations = 5000)
write.csv(rda_full_anova_terms, "rda_full_anova_terms.csv", row.names = FALSE)

# **Step 2: Run Partial RDA (Controlling for Population & Space)**
prda_rda <- rda(gen.imp ~ temp_PC1 + precip_PC1 + urban_PC1 + 
                  Condition(MEM43 + MEM53 + MEM59 + MEM64), 
                data = all_covariates, scale = TRUE)
saveRDS(prda_rda, "prda_rda_model.rds")


# Compute Adjusted R²
prda_rda_r2 <- RsquareAdj(prda_rda)
write.csv(data.frame(R2 = prda_rda_r2$r.squared, Adj_R2 = prda_rda_r2$adj.r.squared), 
          "prda_rda_r2.csv", row.names = FALSE)

signif.prda_rda <- anova.cca(prda_rda, permutations = 5000, parallel=getOption("mc.cores"))
signif.prda_rda

# Compute Permutation Tests
prda_anova <- anova.cca(prda_rda, permutations = 5000)
write.csv(prda_anova, "prda_anova.csv", row.names = FALSE)

rda_env_anova_axis <- anova.cca(prda_rda, by = "axis", permutations = 5000)
write.csv(rda_env_anova_axis, "rda_env_anova_axis.csv", row.names = FALSE)

rda_env_anova_terms <- anova.cca(prda_rda, by = "terms", permutations = 5000)
write.csv(rda_env_anova_terms, "rda_env_anova_terms.csv", row.names = FALSE)

# Save Eigenvalues
eigen_full <- summary(eigenvals(rda.full, model = "constrained"))
write.csv(eigen_full, "rda_full_eigenvalues.csv", row.names = FALSE)

eigen_prda <- summary(eigenvals(prda_rda, model = "constrained"))
write.csv(eigen_prda, "prda_rda_eigenvalues.csv", row.names = FALSE)

# Save R Workspace
save.image("USA_RDA.RData")

############Above is done on hpc cluster#######
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/New_USA_RDA.RData")

summary(eigenvals(prda_rda, model = "constrained")) 

# Check for multicollinearity
#For a given predictor (p), multicollinearity can assessed by computing a score called the variance inflation factor (or VIF), which measures how much the variance of a regression coefficient is inflated due to multicollinearity in the model.
#The smallest possible value of VIF is one (absence of multicollinearity). 
#As a rule of thumb, a VIF value that exceeds 5 or 10 indicates a problematic amount of collinearity (James et al. 2014).
vif.cca(rda.full)
#temp_PC1 precip_PC1  urban_PC1 
#1.939648   1.295147   2.077893 

vif.cca(prda_rda)
#     MEM43      MEM53      MEM59      MEM64   temp_PC1 precip_PC1  urban_PC1 
#1.005957   1.040193   1.039653   1.050951   2.080855   1.188964   2.006161 

#We’ll start with simple triplots from vegan. 
#Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. 
#This scales the SNP and individual scores by the square root of the eigenvalues. See Borcard et al. (2011) or the vegan help for more information on scaling in RDA plots.

pdf("RDA_Scaling.pdf")
plot(rda.full, scaling=3)  # default is axes 1 and 2
dev.off()

pdf("RDA_Scaling_3.pdf")
plot(rda.full, choices = c(1, 3), scaling=3)  # axes 1 and 3
dev.off()

#Here, the SNPs are in red (in the center of each plot), and the individuals are the black circles. The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.
env_pcs_df_area$Location <- factor(env_pcs_df_area$Location )

# Use the Location factor (renamed here as eco) for plotting
eco <- env_pcs_df_area$Location

# Create a vector of colors for each sample by matching the factor levels.
# Convert 'eco' to character to ensure proper indexing.
bg <- location_colors[as.character(eco)]

### Plot for axes 1 & 2 ###
pdf("RDA_axes_1_and_2.pdf", width = 10, height = 8)
plot(rda.full, type = "n", scaling = 3)
# Plot SNPs (species scores)
points(rda.full, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
# Plot sites (SLF) with background colors determined by 'bg'
points(rda.full, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg)
# Add predictor biplot labels
text(rda.full, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
# Add legend using the defined colors for each level of 'eco'
legend("bottomleft", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1,
       pt.bg = location_colors[levels(eco)])
dev.off()

### Plot for axes 1 & 3 ###
pdf("RDA_axes_1_and_3.pdf", width = 10, height = 8)
plot(rda.full, type = "n", scaling = 3, choices = c(1, 3))
points(rda.full, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(rda.full, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg)
text(rda.full, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomleft", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1,
       pt.bg = location_colors[levels(eco)])
dev.off()


#Identify candidate SNPs involved in local adaptation
#We’ll use the loadings of the SNPs in the ordination space to determine which SNPs are candidates for local adaptation. The SNP loadings are stored as species in the RDA object. We’ll extract the SNP loadings from the three significant constrained axes.

load.rda <- scores(rda.full, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

#If we look at histograms of the loadings on each RDA axis, we can see their (relatively normal) distributions. SNPs loading at the center of the distribution are not showing a relationship with the environmental predictors; those loading in the tails are, and are more likely to be under selection as a function of those predictors (or some other predictor correlated with them).
pdf("Loadings_on_RDA1.pdf", width = 7, height = 7)
hist(load.rda[,1], main="Loadings on RDA1")
dev.off()

pdf("Loadings_on_RDA2.pdf", width = 7, height = 7)
hist(load.rda[,2], main="Loadings on RDA2")
dev.off()

pdf("Loadings_on_RDA3.pdf", width = 7, height = 7)
hist(load.rda[,3], main="Loadings on RDA3")
dev.off()



#I’ve written a simple function to identify SNPs that load in the tails of these distributions. 
#We’ll start with a 3 standard deviation cutoff (two-tailed p-value = 0.0027). As with all cutoffs, this can be modified to reflect the goals of the analysis and our tolerance for true positives vs. false positives. For example, if you needed to be very conservative and only identify those loci under very strong selection (i.e., minimize false positive rates), you could increase the number of standard deviations to 3.5 (two-tailed p-value = 0.0005). This would also increase the false negative rate. 
#If you were less concerned with false positives, and more concerned with identifying as many potential candidate loci as possible (including those that may be under weaker selection), you might choose a 2.5 standard deviation cutoff (two-tailed p-value = 0.012).

# Loadings for candidate SNPs
load.rda <- scores(rda.full, choices = 1:3, display = "species")

# Function for outlier selection
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

# We’ll start with a 3 standard deviation cutoff (two-tailed p-value = 0.0027)
z <- 3

cand1 <- outliers(load.rda[,1], z) #1128
cand2 <- outliers(load.rda[,2], z) #582
cand3 <- outliers(load.rda[,3], z) #1133

ncand <- length(cand1) + length(cand2) + length(cand3) 
ncand #2843

# Organize results into dataframe
cand1 <- cbind.data.frame(rep(1, times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# Correlations with environmental PCs
foo <- matrix(nrow=(ncand), ncol=3)
colnames(foo) <- c("Temperature","Precipitation","Urbanization")

pred <- env_pcs_df_area[,12:14] 

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand, foo)

length(cand$snp[duplicated(cand$snp)])  #27 duplicated SNPs

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) # 4 duplicates on axis 2
table(foo[foo[,1]==3,2]) # 23 duplicates on axis 3

# Remove duplicate SNPs
cand <- cand[!duplicated(cand$snp),] #2816 in the end

write.csv(cand, "Full_RDA_Candidate_SNPs.csv", row.names = F)

#Next, we’ll see which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor) 
#Precipitation   Temperature  Urbanization 
#796          1372           648

write.csv(cand, "final.cand_FullRDA.3PC.csv", row.names = F)


#Plot the SNPs
#Let’s look at RDA plots again, but this time focus in on the SNPs in the ordination space. We’ll color code the SNPs based on the predictor variable that they are most strongly correlated with. 

sel <- cand$snp
env <- cand$predictor
env[env == "Temperature"]   <- "#1f78b4"
env[env == "Precipitation"] <- "#c861a5"
env[env == "Urbanization"]  <- "#33a02c"

# Initialize colors for all SNPs as gray (non-candidates)
col_pred <- rep('#f1eef6', length(rownames(rda.full$CCA$v)))  # Default gray
names(col_pred) <- rownames(rda.full$CCA$v)  # Assign row names

# Assign colors to candidate SNPs
for (i in 1:length(sel)) {
  foo <- match(sel[i], names(col_pred))  # Ensure proper matching
  if (!is.na(foo)) {
    col_pred[foo] <- env[i]  # Assign correct color
  }
}

# Set non-candidates to transparent
empty <- col_pred
empty[empty == "#f1eef6"] <- rgb(0,1,0, alpha=0)  # Transparent for non-candidates
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")

# Define legend colors
bg <- c("#1f78b4", "#c861a5", "#33a02c")

# axes 1 & 2
pdf("SLF_RDA_axes_1_and_2.pdf", width = 10, height = 8)
plot(rda.full, type="n", scaling=1,xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (36.88%)', ylab= 'RDA2 (32.81%)')
points(rda.full, display="species", pch=21, cex=1.5, col="gray32", bg=col_pred, scaling=1)
points(rda.full, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=1)
text(rda.full, scaling=1, display="bp",col="#191970", cex=1.2)
legend("bottomright", legend = c("Temperature", "Precipitation", "Urbanization"),
       bty = "n", col = "gray32", pch = 21, cex = 2, pt.bg = bg)
dev.off()

# axes 1 & 3
pdf("SLF_RDA_axes_1_and_3.pdf", width = 10, height = 8)
plot(rda.full, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1),xlab= 'RDA1 (36.88%)', ylab= 'RDA3 (30.31%)', choices=c(1,3))
points(rda.full, display="species", pch=21, cex=1.5, col="gray32", bg=col_pred, scaling=1,choices=c(1,3))
points(rda.full, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=1)
text(rda.full, scaling=1, display="bp", col="#191970", cex=1.2, choices=c(1,3))
legend("bottomright", legend = c("Temperature", "Precipitation", "Urbanization"),
       bty = "n", col = "gray32", pch = 21, cex = 1.5, pt.bg = bg)
dev.off()

cand <- cand[order(-abs(cand$correlation)), ]
cand$ID <- substr(cand$snp, 1, nchar(cand$snp) - 2)

# Remove the allele suffix (_C, _T, _G, _A) from `cand_prda$snp`
cand$ID <- gsub("_[ACGT]$", "", cand$snp)  # Remove last character if it's an allele

#load the SNP annotation file
annot_combined <- read.csv("Final_SLF_Annotation_Info.csv")  
annot_combined$SNP <- gsub("_[ACGT]$", "", annot_combined$SNP) 
cand$ID <- gsub("_[ACGT]$", "", cand$ID)

# Check unique formats before merging
head(cand$ID)
head(annot_combined$SNP)

# Merge annotation with candidate SNP results
SLF_annotation_FullRDA <- merge(annot_combined, cand, by.x = "SNP", by.y = "ID", all.x = FALSE)  # Use FALSE to get only matches

# Remove rows where `predictor` is missing
SLF_annotation_FullRDA <- SLF_annotation_FullRDA[!is.na(SLF_annotation_FullRDA$predictor), ]

# Save updated annotated candidate SNPs
write.csv(SLF_annotation_FullRDA, "SLF_annotation_FullRDA.csv", row.names = FALSE)

# Check the frequency table for predictors
table(SLF_annotation_FullRDA$predictor)
#Precipitation   Temperature  Urbanization 
#796          1373           648


#https://popgen.nescent.org/2018-03-27_RDA_GEA.html


#In summary, Redundancy Analysis (RDA) is a highly effective genotype-environment association method, which can identify unique haplotypes associated with the multivariate environment. In a simulation study, RDA showed a superior combination of low false positive and high true positive rates across weak, moderate, and strong multilocus selection. These results were robust across the levels of population structure, demographic histories, sampling designs, and sample sizes tested (Forester et al., 2018).



###Part2b: partial RDA

## Load data
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/USA_RDA.RData")

# **Step 2: Run Partial RDA (Controlling for Population & Space)**

summary(eigenvals(prda_rda, model = "constrained")) 

library(vegan)
# Generate simple triplots for Partial RDA
pdf("Partial_RDA_Scaling.pdf")
plot(prda_rda, scaling=3)  # default axes 1 and 2
dev.off()

pdf("Partial_RDA_Scaling_3.pdf")
plot(prda_rda, choices = c(1, 3), scaling=3)  # axes 1 and 3
dev.off()

# Assign locations for color mapping
env_pcs_df_area$Location <- factor(env_pcs_df_area$Location)
eco <- env_pcs_df_area$Location

# Create a vector of colors for each sample by matching factor levels
bg_prda <- location_colors[as.character(eco)]

# **Plot Partial RDA - Axes 1 & 2**
pdf("Partial_RDA_axes_1_and_2.pdf", width = 10, height = 8)
plot(prda_rda, type = "n", scaling = 3)
points(prda_rda, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(prda_rda, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg_prda)
text(prda_rda, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomleft", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1,
       pt.bg = location_colors[levels(eco)])
dev.off()

# **Plot Partial RDA - Axes 1 & 3**
pdf("Partial_RDA_axes_1_and_3.pdf", width = 10, height = 8)
plot(prda_rda, type = "n", scaling = 3, choices = c(1, 3))
points(prda_rda, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(prda_rda, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg_prda)
text(prda_rda, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomleft", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1,
       pt.bg = location_colors[levels(eco)])
dev.off()

# Step 2: Identify Candidate SNPs in Partial RDA
# Extract SNP loadings **only from RDA1**
# We’ll start with a 3 standard deviation cutoff (two-tailed p-value = 0.0027)
z <- 3

# Function to identify outlier SNPs (3 standard deviations)
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

# Extract SNP loadings from first three constrained axes
load_prda <- scores(prda_rda, choices=c(1:3), display="species")


# **Plot Loadings Histograms**
pdf("Loadings_on_Partial_RDA1.pdf", width = 7, height = 7)
hist(load_prda[,1], main="Loadings on Partial RDA1")
dev.off()

pdf("Loadings_on_Partial_RDA2.pdf", width = 7, height = 7)
hist(load_prda[,2], main="Loadings on Partial RDA2")
dev.off()

pdf("Loadings_on_Partial_RDA3.pdf", width = 7, height = 7)
hist(load_prda[,3], main="Loadings on Partial RDA3")
dev.off()

# Apply 3 SD cutoff
cand1_prda <- outliers(load_prda[,1], z) # RDA1 1254 SNPs
cand2_prda <- outliers(load_prda[,2], z) # RDA2 502 SNPs
cand3_prda <- outliers(load_prda[,3], z) # RDA3 945 SNPs

# Count total number of candidate SNPs
ncand_prda <- length(cand1_prda) + length(cand2_prda) + length(cand3_prda)
ncand_prda # Number of candidate SNPs, 2701 SNPs

# Organize results into dataframe
cand1_prda <- cbind.data.frame(rep(1, times=length(cand1_prda)), names(cand1_prda), unname(cand1_prda))
cand2_prda <- cbind.data.frame(rep(2, times=length(cand2_prda)), names(cand2_prda), unname(cand2_prda))
cand3_prda <- cbind.data.frame(rep(3, times=length(cand3_prda)), names(cand3_prda), unname(cand3_prda))

colnames(cand1_prda) <- colnames(cand2_prda) <- colnames(cand3_prda) <- c("axis","snp","loading")

cand_prda <- rbind(cand1_prda, cand2_prda, cand3_prda)
cand_prda$snp <- as.character(cand_prda$snp)

# Correlations with environmental PCs
foo_prda <- matrix(nrow=(ncand_prda), ncol=3)
colnames(foo_prda) <- c("Temperature","Precipitation","Urbanization")

for (i in 1:length(cand_prda$snp)) {
  nam <- cand_prda[i,2]
  snp.gen <- gen.imp[,nam]
  foo_prda[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand_prda <- cbind.data.frame(cand_prda, foo_prda)

length(cand_prda$snp[duplicated(cand_prda$snp)])  #26 duplicates

foo <- cbind(cand_prda$axis, duplicated(cand_prda$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #17 duplicates on axis 2
table(foo[foo[,1]==3,2]) #9 duplicates on axis 3

# Remove duplicate SNPs
cand_prda <- cand_prda[!duplicated(cand_prda$snp),] 

write.csv(cand_prda, "Partial_RDA_Candidate_SNPs.csv", row.names = F)

#Next, we’ll see which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand_prda$snp)) {
  bar <- cand_prda[i,]
  cand_prda[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand_prda[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(cand_prda)[7] <- "predictor"
colnames(cand_prda)[8] <- "correlation"

table(cand_prda$predictor) 

write.csv(cand_prda, "final.cand_prda.3PC.csv", row.names = F)


# Step 3: Visualize Candidate SNPs in Partial RDA
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

# Assign colors to candidate SNPs
for (i in 1:length(sel_prda)) {
  foo <- match(sel_prda[i], names(col_pred_prda))  # Ensure proper matching
  if (!is.na(foo)) {
    col_pred_prda[foo] <- env_prda[i]  # Assign correct color
  }
}

# Set non-candidates to transparent
empty <- col_pred_prda
empty[empty == "#f1eef6"] <- rgb(0,1,0, alpha=0)  # Transparent for non-candidates
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")

# Define legend colors
bg_prda <- c("#1f78b4", "#c861a5", "#33a02c")

# **Plot SNPs for Partial RDA (Axes 1 & 2)**
pdf("Partial_RDA_SNP_axes_1_and_2.pdf", width = 10, height = 8)
plot(prda_rda, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1),xlab= 'RDA1 (37.0%)', ylab= 'RDA2 (32.3%)')
points(prda_rda, display="species", pch=21, cex=1.5, col="gray32", bg=col_pred_prda, scaling=1)
text(prda_rda, scaling=1, display="bp", col="#191970", cex=1.2)
legend("bottomright", legend = c("Temperature", "Precipitation", "Urbanization"), bty = "n", col = "gray32", pch = 21, cex = 1.5, pt.bg = bg_prda)
dev.off()

# **Plot SNPs for Partial RDA (Axes 1 & 3)**
pdf("Partial_RDA_SNP_axes_1_and_3.pdf", width = 10, height = 8)
plot(prda_rda, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1), xlab= 'RDA1 (32.9%)', ylab= 'RDA3 (30.1%)', choices=c(1,3))
points(prda_rda, display="species", pch=21, cex=1.5, col="gray32", bg=col_pred_prda, scaling=1, choices=c(1,3))
text(prda_rda, scaling=1, display="bp", col="#191970", cex=1.2, choices=c(1,3))
legend("bottomright", legend = c("Temperature", "Precipitation", "Urbanization"),bty = "n", col = "gray32", pch = 21, cex = 1.5, pt.bg = bg_prda)
dev.off()

cand_prda <- cand_prda [order(-abs(cand_prda $correlation)), ]
cand_prda$ID <- substr(cand_prda $snp, 1, nchar(cand_prda $snp) - 2)

annot_combined <- read.csv('Final_SLF_Annotation_Info.csv')

# Ensure SNP IDs match between `cand_prda` and `annot_combined`

# Remove the allele suffix (_C, _T, _G, _A) from `cand_prda$snp`
cand_prda$ID <- gsub("_[ACGT]$", "", cand_prda$ID)  # Remove last character if it's an allele
annot_combined$SNP <- gsub("_[ACGT]$", "", annot_combined$SNP) 

# Check unique formats before merging
head(cand_prda$ID)
head(annot_combined$SNP)

# Merge annotation with candidate SNP results
prda_SLF_annotation_RDA <- merge(annot_combined, cand_prda, by.x = "SNP", by.y = "ID", all.x = FALSE)  # Use FALSE to get only matches

# Remove rows where `predictor` is missing
prda_SLF_annotation_RDA <- prda_SLF_annotation_RDA[!is.na(prda_SLF_annotation_RDA$predictor), ]

# Save updated annotated candidate SNPs
write.csv(prda_SLF_annotation_RDA, "prda_SLF_annotation.csv", row.names = FALSE)

# Check the frequency table for predictors
table(prda_SLF_annotation_RDA$predictor)
#Precipitation   Temperature  Urbanization 
#848          1262           566 

#https://popgen.nescent.org/2018-03-27_RDA_GEA.html


#Part 3: Find the overlap for RDA and pcadpt
# Load required libraries
library(VennDiagram)
library(ggplot2)
library(dplyr)

# Read in data
srda_outliers <- read.csv("SLF_annotation_FullRDA.csv")$SNP  # Simple RDA outliers, 2817
prda_outliers <- read.csv("prda_SLF_annotation.csv")$SNP  # Partial RDA outliers, 2675
pcadapt_outliers <- read.csv("Pcadpt_fdr_USA_Only_SigFDR.csv")$SNP  # PCAdapt outliers, 7349

# Identify overlaps between methods
srda_prda <- intersect(srda_outliers, prda_outliers)
srda_pcadapt <- intersect(srda_outliers, pcadapt_outliers)
prda_pcadapt <- intersect(prda_outliers, pcadapt_outliers)
all_overlap <- Reduce(intersect, list(srda_outliers, prda_outliers, pcadapt_outliers))

# Check overlap counts
length(srda_prda)      # Overlap between simple RDA and partial RDA
length(srda_pcadapt)   # Overlap between simple RDA and PCAdapt
length(prda_pcadapt)   # Overlap between partial RDA and PCAdapt
length(all_overlap)    # SNPs found in all three methods

# Create a named list for overlaps
overlap_list <- list(
  SRDA_PRDA = srda_prda,
  SRDA_PCAdapt = srda_pcadapt,
  PRDA_PCAdapt = prda_pcadapt,
  All_Methods = all_overlap
)

# Convert to long format dataframe
shared_snp_df <- bind_rows(lapply(names(overlap_list), function(method) {
  data.frame(
    Method = method,
    SNP = overlap_list[[method]],
    stringsAsFactors = FALSE
  )
}))

# Save overlap data
write.csv(shared_snp_df, "shared_outliers.csv", row.names = FALSE)

# Save all unique outliers from all methods
all_outliers <- unique(c(srda_outliers, prda_outliers, pcadapt_outliers))
write.csv(data.frame(SNP = all_outliers), "all_unique_outliers.csv", row.names = FALSE)

# Create a temporary directory for the Venn diagram
dir.create("venn_temp", showWarnings = FALSE)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(
    sRDA = srda_outliers,
    pRDA = prda_outliers,
    PCAdapt = pcadapt_outliers
  ),
  category.names = c("sRDA", "pRDA", "PCAdapt"),
  filename = "venn_temp/outliers_venn_diagram.png",
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

# Create a summary table of the overlaps
overlap_summary <- data.frame(
  Comparison = c("sRDA", "pRDA", "PCAdapt", 
                 "sRDA & pRDA", "sRDA & PCAdapt", "pRDA & PCAdapt", 
                 "All methods"),
  SNP_Count = c(length(srda_outliers), length(prda_outliers), length(pcadapt_outliers),
                length(srda_prda), length(srda_pcadapt), length(prda_pcadapt),
                length(all_overlap))
)

# Print summary table
print(overlap_summary)

# Save summary table
write.csv(overlap_summary, "overlap_summary.csv", row.names = FALSE)

# Bonus: Create an UpSet plot (alternative to Venn diagram for better visualization)
# If you have the UpSetR package installed
if(requireNamespace("UpSetR", quietly = TRUE)) {
  library(UpSetR)
  
  # Create a list for the UpSet plot
  input_list <- list(
    SimpleRDA = srda_outliers,
    PartialRDA = prda_outliers,
    PCAdapt = pcadapt_outliers
  )
  
  # Create the UpSet plot
  upset_plot <- upset(fromList(input_list), 
                      nsets = 3, 
                      nintersects = NA, 
                      order.by = "freq",
                      keep.order = TRUE)
  
  # Save the plot
  pdf("upset_plot_outliers.pdf", width = 8, height = 6)
  print(upset_plot)
  dev.off()
}

# Load the list of 141 shared SNPs found in all three methods
all_overlap_snps <- read.csv("shared_outliers.csv") %>%
  filter(Method == "All_Methods") %>%
  pull(SNP)

# Assuming you have annotation information in these files
srda_annotation <- read.csv("SLF_annotation_FullRDA.csv")
prda_annotation <- read.csv("prda_SLF_annotation.csv")
pcadapt_annotation <- read.csv("Pcadpt_fdr_USA_Only_SigFDR.csv")

shared_snp_annotation <- srda_annotation %>%
  filter(SNP %in% all_overlap_snps) %>%
  select(SNP, Gene_Name) %>%  # Adjust column names as needed
  distinct()

write.csv(shared_snp_annotation, "Slection_shared_snp_annotation.csv", row.names = F)

# Merge annotation with candidate SNP results
ALL_overlap_annotation_SNP <- merge(annot_combined, shared_snp_annotation, by = "SNP", all.x = FALSE)  # Use FALSE to get only matches

# Save updated annotated candidate SNPs
write.csv(ALL_overlap_annotation_SNP, "ALL_overlap_annotation_SNP.csv", row.names = FALSE)
save.image("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Final_USA_RDA.RData")






# Load necessary libraries
library(dplyr)

# Load the list of 141 shared SNPs found in all three methods
all_overlap_snps <- read.csv("shared_outliers.csv") %>%
  filter(Method == "All_Methods") %>%
  pull(SNP)

# Load your annotation files with predictor information
srda_annotation <- read.csv("SLF_annotation_FullRDA.csv")
prda_annotation <- read.csv("prda_SLF_annotation.csv")

# Extract SNP and predictor information from both RDA analyses
srda_predictors <- srda_annotation %>%
  filter(SNP %in% all_overlap_snps) %>%
  select(SNP, Gene_Name, predictor, correlation) %>%
  rename(srda_predictor = predictor, 
         srda_correlation = correlation)

prda_predictors <- prda_annotation %>%
  filter(SNP %in% all_overlap_snps) %>%
  select(SNP, predictor, correlation) %>%
  rename(prda_predictor = predictor,
         prda_correlation = correlation)

# Merge the predictor information from both analyses
shared_snps_with_predictors <- srda_predictors %>%
  left_join(prda_predictors, by = "SNP")

# Load full annotation information
annot_combined <- read.csv("Final_SLF_Annotation_Info.csv")

# Merge with complete annotation information
ALL_overlap_annotation_SNP <- annot_combined %>%
  right_join(shared_snps_with_predictors, by = c("SNP", "Gene_Name"))

# Create summary counts of predictors
srda_predictor_counts <- shared_snps_with_predictors %>%
  count(srda_predictor) %>%
  rename(Predictor = srda_predictor, 
         Simple_RDA_Count = n)

prda_predictor_counts <- shared_snps_with_predictors %>%
  count(prda_predictor) %>%
  rename(Predictor = prda_predictor, 
         Partial_RDA_Count = n)

# Combine counts into one table
predictor_summary <- full_join(srda_predictor_counts, prda_predictor_counts, by = "Predictor")

# Check for consistency between methods
consistency_check <- shared_snps_with_predictors %>%
  mutate(consistent = srda_predictor == prda_predictor) %>%
  count(consistent)

consistent_predictors <- shared_snps_with_predictors %>%
  filter(srda_predictor == prda_predictor) %>%
  count(srda_predictor) %>%
  rename(Predictor = srda_predictor, 
         Consistent_Count = n)

# Add this to the summary table
predictor_summary <- predictor_summary %>%
  left_join(consistent_predictors, by = "Predictor")

# Save all results
write.csv(shared_snps_with_predictors, "shared_141_SNPs_with_predictors.csv", row.names = FALSE)
write.csv(ALL_overlap_annotation_SNP, "ALL_overlap_annotation_SNP.csv", row.names = FALSE)
write.csv(predictor_summary, "shared_141_SNPs_predictor_summary.csv", row.names = FALSE)

# Print summary information
cat("Number of shared SNPs:", length(all_overlap_snps), "\n")
cat("\nPredictor counts in Simple RDA:\n")
print(srda_predictor_counts)
cat("\nPredictor counts in Partial RDA:\n")
print(prda_predictor_counts)
cat("\nConsistency between methods:\n")
print(consistency_check)
cat("\nDetailed predictor summary:\n")
print(predictor_summary)

# Create a stacked bar plot to visualize predictor distribution
library(ggplot2)
library(tidyr)

# Reshape data for plotting
plot_data <- shared_snps_with_predictors %>%
  select(SNP, srda_predictor, prda_predictor) %>%
  pivot_longer(cols = c(srda_predictor, prda_predictor),
               names_to = "method",
               values_to = "predictor") %>%
  mutate(method = ifelse(method == "srda_predictor", "Simple RDA", "Partial RDA"))

# Create the plot
predictor_plot <- ggplot(plot_data, aes(x = method, fill = predictor)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Temperature" = "#1f78b4", 
                               "Precipitation" = "#c861a5",
                               "Urbanization" = "#33a02c")) +
  labs(title = "Environmental Predictors Associated with 141 Shared Selection Candidates",
       x = NULL,
       y = "Number of SNPs",
       fill = "Environmental\nPredictor") +
  theme_minimal()

# Save the plot
ggsave("shared_141_SNPs_predictor_distribution.pdf", predictor_plot, width = 8, height = 6)

# Create a Sankey/alluvial plot to show changes in predictor assignment between methods
if(requireNamespace("ggalluvial", quietly = TRUE)) {
  library(ggalluvial)
  
  # Count combinations of predictors between methods
  predictor_flows <- shared_snps_with_predictors %>%
    count(srda_predictor, prda_predictor) %>%
    rename(Simple_RDA = srda_predictor,
           Partial_RDA = prda_predictor,
           Count = n)
  
  # Create alluvial plot
  alluvial_plot <- ggplot(predictor_flows,
                          aes(axis1 = Simple_RDA, axis2 = Partial_RDA,
                              y = Count)) +
    geom_alluvium(aes(fill = Simple_RDA), alpha = 0.8) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Simple RDA", "Partial RDA"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = c("Temperature" = "#1f78b4", 
                                 "Precipitation" = "#c861a5",
                                 "Urbanization" = "#33a02c")) +
    labs(title = "Shifts in Environmental Predictor Assignment Between Methods",
         subtitle = "For 141 SNPs identified by all three selection methods",
         y = "Number of SNPs",
         fill = "Simple RDA\nPredictor") +
    theme_minimal()
  
  # Save the plot
  ggsave("shared_141_SNPs_predictor_flow.pdf", alluvial_plot, width = 9, height = 7)
}

# Save workspace
save.image("shared_141_SNPs_analysis.RData")


#Part 4: Find the overlap for pRDA and pcadpt


#Part 5: Plot the highlight 185 SNPs 
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Final_USA_RDA.RData")
load("/scratch/work/winchell/Lanternfly/WGS/16_NoDownsample_VCF_Analysis/Final_annotated_VCF/Separated/Highlighted_Shared_SNPs_Results.RData")

# Load shared SNPs (assumes column name "SNP")
shared_snps <- read.csv("185overlap_prda_pcadapt_full_info.csv")
shared_snp_ids_raw <- shared_snps$SNP
shared_snp_ids <- gsub("_[ACGT]+$", "", shared_snp_ids_raw)

# Clean rownames of RDA SNPs
snp_ids_raw <- rownames(prda_rda$CCA$v)
snp_ids_clean <- gsub("_[ACGT]+$", "", snp_ids_raw)
rownames(prda_rda$CCA$v) <- snp_ids_clean

# Remove both allele suffixes (_T_G or _A_A) — leave just chr and position
clean_snp_ids <- function(x) sub("_[ACGT]+(_[ACGT]+)?$", "", x)

# Apply to candidate SNPs and to the rownames of RDA loadings
cand_snp_ids <- clean_snp_ids(cand_prda$snp)
rownames(prda_rda$CCA$v) <- clean_snp_ids(rownames(prda_rda$CCA$v))
cand_env <- cand_prda$predictor

# Sanity check
stopifnot(all(cand_snp_ids %in% rownames(prda_rda$CCA$v)))
setdiff(cand_snp_ids, rownames(prda_rda$CCA$v)) %>% head()

# Map predictors to colors
color_map <- setNames(c("#1f78b4", "#c861a5", "#33a02c"),
                      c("Temperature", "Precipitation", "Urbanization"))
cand_colors <- color_map[as.character(cand_env)]

# Initialize color and shape mappings
snp_color <- rep("#f1eef6", length(rownames(prda_rda$CCA$v)))
names(snp_color) <- rownames(prda_rda$CCA$v)

snp_shape <- rep(21, length(rownames(prda_rda$CCA$v)))  # circle default
names(snp_shape) <- rownames(prda_rda$CCA$v)

snp_border <- rep("gray32", length(rownames(prda_rda$CCA$v)))
names(snp_border) <- rownames(prda_rda$CCA$v)

# Apply colors and shapes to candidate SNPs
for (i in seq_along(cand_snp_ids)) {
  snp <- cand_snp_ids[i]
  if (snp %in% names(snp_color)) {
    snp_color[snp] <- cand_colors[i]
    if (snp %in% shared_snp_ids) {
      snp_shape[snp] <- 22  # square for shared
      snp_border[snp] <- "black"
    }
  }
}

# Diagnostic check
cat("Number of candidate SNPs:", length(cand_snp_ids), "\n")
cat("Matched with RDA SNPs:", sum(cand_snp_ids %in% rownames(prda_rda$CCA$v)), "\n")
cat("Shared SNPs found:", sum(shared_snp_ids %in% rownames(prda_rda$CCA$v)), "\n")

# Get coordinates
snp_coords <- scores(prda_rda, choices = c(1, 2), display = "species", scaling = 1)

# Plot (axes 1 & 2)
pdf("Partial_RDA_SNP_Shapes_1_2_red.pdf", width = 10, height = 8)
#pdf("Partial_RDA_SNP_Shapes_1_2.pdf", width = 10, height = 8)
plot(prda_rda, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "RDA1 (37.0%)", ylab = "RDA2 (32.3%)")

for (i in seq_len(nrow(snp_coords))) {
  snp <- rownames(snp_coords)[i]
  if (snp_color[snp] != "#f1eef6") {
    points(snp_coords[i, 1], snp_coords[i, 2],
           pch = snp_shape[snp],
           bg = snp_color[snp],
           col = snp_border[snp],
           cex = 2)
  }
}
#text(prda_rda, scaling = 1, display = "bp", col = "#191970", cex = 1.2)

legend("bottomright",
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"),
       pt.bg = c("#1f78b4", "#c861a5", "#33a02c", "gray"),
       col = c("gray32", "gray32", "gray32", "red"),
       pch = c(21, 21, 21, 22),
       bty = "n",
       cex = 2)
dev.off()

# Get coordinates
snp_coords <- scores(prda_rda, choices = c(1, 3), display = "species", scaling = 1)

pdf("Partial_RDA_SNP_Shapes_1_3_red.pdf", width = 10, height = 8)
plot(prda_rda, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "RDA1 (37.0%)", ylab = "RDA3 (30.1%)",choices=c(1,3))

for (i in seq_len(nrow(snp_coords))) {
  snp <- rownames(snp_coords)[i]
  
  if (snp_color[snp] != "#f1eef6") {
    is_shared <- snp %in% shared_snp_ids
    snp_pch <- ifelse(is_shared, 21)  # 22 = square, 21 = circle
    snp_col <- ifelse(is_shared, "red", "gray32")  # border color
    
    points(snp_coords[i, 1], snp_coords[i, 2],
           pch = snp_pch,
           bg = snp_color[snp],
           col = snp_col,
           cex = 2,
           lwd = 2)
  }
}
legend("bottomright",
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"),
       pt.bg = c("#1f78b4", "#c861a5", "#33a02c", "gray"),
       col = c("gray32", "gray32", "gray32", "red"),
       pch = c(21, 21, 21, 21),
       bty = "n",
       pt.cex = 2,
       pt.lwd = 2,
       cex = 1.8)
dev.off()


# Axes 1 & 2
pdf("Partial_RDA_SNP_redBorderShared_1_3.pdf", width = 10, height = 8)
snp_coords <- scores(prda_rda, choices = c(1, 3), display = "species", scaling = 1)
plot(prda_rda, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "RDA1 (37.0%)",ylab = "RDA3 (30.1%)",choices=c(1,3))

for (i in seq_len(nrow(snp_coords))) {
  snp <- rownames(snp_coords)[i]
  
  if (snp_color[snp] != "#f1eef6") {
    is_shared <- snp %in% shared_snp_ids
    border_col <- ifelse(is_shared, "red", "gray32")
    
    points(snp_coords[i, 1], snp_coords[i, 2],
           pch = 21,
           bg = snp_color[snp],
           col = border_col,
           cex = 2,
           lwd = 2)
  }
}

legend("bottomright",
       legend = c("Temperature", "Precipitation", "Urbanization", "Shared SNPs"),
       pt.bg = c("#1f78b4", "#c861a5", "#33a02c", "red"),
       col = c("gray32", "gray32", "gray32", "black"),
       pch = 21,
       pt.cex = 2,
       pt.lwd = 2,
       bty = "n",
       cex = 1.8)

dev.off()

# Export SNP loadings
snp_loadings <- scores(prda_rda, choices = 1:3, display = "species", scaling = 1)
snp_df <- as.data.frame(snp_loadings)
snp_df$shared <- ifelse(rownames(snp_df) %in% shared_snp_ids, "yes", "no")
snp_df$environment <- ifelse(rownames(snp_df) %in% cand_snp_ids,
                             cand_env[match(rownames(snp_df), cand_snp_ids)],
                             "non-candidate")
write.csv(snp_df, "snp_loadings_cleaned.csv", row.names = TRUE)


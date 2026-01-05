library(vegan)
library(apcluster)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)


#set location of your dataset and resulting data - change this for your own location
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

###################
###   Setup     ###
###################

figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"

# Load MAG abundance table
metaT <- read_excel('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_abundance_table.xlsx')

# Drop Condensed Tannin samples
ct_samples <- c('STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127', 'GTDB')

metaT <- metaT %>% select(-any_of(ct_samples))

# Drop rows where sum of all numeric columns is 0
numeric_cols <- names(metaT)[-1]
metaT <- metaT %>% 
  filter(rowSums(select(., all_of(numeric_cols))) > 0)

# Prepare dataframe
df <- as.data.frame(metaT)
rownames(df) <- df[[1]]  # First column as row names (MAG IDs)
df <- df[, -1]  # Remove first column

metadata <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")
metadata <- metadata %>%
  mutate(treatment = case_when(
    treatment == "unamended" ~ "U",
    treatment == "CT" ~ "T",
    treatment == "catechin" ~ "C",
    TRUE ~ treatment
  )) %>% filter(treatment != 'T')



treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')

rownames(df) <- df$vOTU
df$vOTU <- NULL
colnames(df) <- treatment_timepoint

# make unique column names for df to investigate separately
colnames(df) <- make.unique(names(df))

#log transformation
#note this is ln
t_df <- t(df) %>% as_tibble(rownames = NA)
# log_df <- log(t_df + 1)
log_df <- decostand(t_df, method = "hellinger")

#bray curtis
bray_dist <- vegdist(log_df, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)


#some data will need to be a matrix or dataframe later
bray_matrix_df <- as.matrix(bray_dist)
# write.table(bray_matrix_df, file = "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU_aggregated_bray_curtis.csv")

# pdf(file = paste(figDir, "SX_df_bray.pdf", sep=''), width = 10, height = 10);
# heatmap(bray_matrix_df, main = "Bray Curtis")
# dev.off()

# Plot boxplots of bray-curtis distance between samples from different treatments but same day to compare to MAGs
# Are viruses more sensitive to changes in the environment than the microbes?

# Convert to data frame for easier manipulation
bray_df <- as.data.frame(as.table(bray_matrix_df))
names(bray_df) <- c("Sample1", "Sample2", "Distance")

# Extract metadata from sample names
extract_metadata <- function(sample) {
  parts <- strsplit(as.character(sample), "_")[[1]]
  treatment <- substr(parts[1], 1, 1)
  day <- parts[2]
  # Remove any replicate information from day
  day <- gsub("\\..*", "", day)
  return(data.frame(Treatment = treatment, Day = day))
}

# Apply to both sample columns
meta1 <- do.call(rbind, lapply(bray_df$Sample1, extract_metadata))
meta2 <- do.call(rbind, lapply(bray_df$Sample2, extract_metadata))

# Add to main data frame
bray_df$Treatment1 <- meta1$Treatment
bray_df$Day1 <- meta1$Day
bray_df$Treatment2 <- meta2$Treatment
bray_df$Day2 <- meta2$Day

# Filter for U vs C comparisons on the same day, excluding day0
# Now we'll specifically select only U vs C (not C vs U) to avoid duplicates
bray_df_filtered <- bray_df[
  bray_df$Day1 == bray_df$Day2 & 
    bray_df$Day1 != "day0" & 
    bray_df$Treatment1 == "U" & 
    bray_df$Treatment2 == "C", ]

# Create a factor for days (in correct order)
bray_df_filtered$Day <- factor(bray_df_filtered$Day1, 
                               levels = c("day7", "day14", "day21", "day35"))

# Create the plot with both boxplots and points
ggplot(bray_df_filtered, aes(x = Day, y = Distance)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +  # Hide default outliers
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.6, color = "steelblue") +
  labs(title = "Bray-Curtis distances between unamended and catechin samples",
       y = "Bray-Curtis Distance",
       x = "Day") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title




# Determine Optimal number of clusters #

# Elbow method
fviz_nbclust(bray_matrix_df, kmeans, method = "wss") +
  #   geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(bray_matrix_df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)

#calculate gap statistic based on number of clusters
#gap_stat <- clusGap(bray_matrix_df,
#                    FUN = kmeans,
#                    nstart = 25,
#                    K.max = 20,
#                    B = 500)

#plot number of clusters vs. gap statistic
#fviz_gap_stat(gap_stat)



#affinity propagation
df_ap_clust <- apcluster(negDistMat(r=2), bray_matrix_df, details = TRUE, seed=42)
plot(df_ap_clust)
df_ap_clust

#r here is a power used to scale the data and allow for more easily recognizable separations
# pdf(file = paste(figDir, "SX_df_ap_cluster.pdf", sep=''), width = 10, height = 10);
# heatmap(df_ap_clust)
# dev.off()
# df_ap_clust

# df_ap_clustK <- apclusterK(negDistMat(r=2), bray_matrix_df, K=4, details = TRUE, seed=42)
# df_ap_clustK

km <- kmeans(bray_matrix_df, centers = 5, nstart = 25)
km

df_types <- data.frame()

types <- data.frame()
#add each cluster to tibble
df_clusters <- 6
df_types <- data.frame()  # Initialize an empty data frame

for (i in 1:df_clusters) {
  df_temp_types <- data.frame("site" = names(df_ap_clust[[i]]),
                              "type" = paste0("Cluster ", i))
  
  df_types <- bind_rows(df_types, df_temp_types)
}


df_type_no_change <- df_types
sites_df <- df_types$site
df_types <- as.data.frame(df_types)
df_types <- as.data.frame(df_types[, -1])
rownames(df_types) <- sites_df
names(df_types)[1] <- "type"
fac_df <-factor(df_types$type)
cols_df <- c("red", "blue", "orange", "green", "violet", "gray", "cyan")

# Plot ordination
# set seed
set.seed(42)

# ordination
NMDS_df <- metaMDS(log_df, distance = "bray", k=2, trymax = 200, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_df

# the stress value is printed here
# noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
# k is the value of dimensions you want the data to be constrained to for a reported stress value
# we must reach convergence

# pdf(file = paste(figDir, "SX_NMDS_stress_df.pdf", sep=''), width = 12, height = 9);
# stress_df <- stressplot(NMDS_df)
# plot(stress_df)
# dev.off()
# dev.off()
# stress plot is just a relationship between ordination distance and dissimilarity
# goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

## plot ordination ##

# wrangle metadata
days <- sapply(strsplit(rownames(log_df), "day"), function(x) x[2])
days <- sapply(strsplit(days, '\\.'), function(x) x[1])
treatments <- sapply(strsplit(rownames(log_df), "_"), function(x) x[1])

df_types2 = as.data.frame(cbind(days, treatments))
colnames(df_types2) <- c('day', 'treatment')
rownames(df_types2) <- rownames(log_df)

#load the type data
fac_df <-factor(df_types2$treatment)
pch_fac_df <- factor(df_types2$day, levels=c('0', '7', '14', '21', '35'))
pch_df <- c(4, 15, 17, 18, 19)
cols_df <- c("#FC9B2D", "#7ACBC3") # C U

## With stat_ellipse for each affinity propagation cluster ##

# Ensure NMDS results are in a data frame
NMDS_points <- as.data.frame(scores(NMDS_df, display = "sites"))

# Add metadata (matching by row names)
NMDS_points$treatment <- df_types2$treatment
NMDS_points$day <- df_types2$day
matched_clusters <- df_types[match(rownames(NMDS_points), rownames(df_types)), ]
NMDS_points$cluster <- matched_clusters  # Add clusters from df_types

# Define shape and color mappings
pch_df <- c(4, 15, 17, 18, 19)
names(pch_df) <- c('0', '7', '14', '21', '35')

cols_df <- c("#FC9B2D", "#7ACBC3") # Catechin, Control, Unamended
names(cols_df) <- levels(factor(df_types2$treatment))

# Compute convex hulls
hull_data <- NMDS_points %>%
  group_by(cluster) %>%
  slice(chull(NMDS1, NMDS2))

# Run PCA
PCA_df <- prcomp(log_df, center = FALSE, scale. = FALSE)
PCA_scores <- as.data.frame(PCA_df$x[, 1:2])  # Get first 2 components
colnames(PCA_scores) <- c("PC1", "PC2")
explained_var <- summary(PCA_df)$importance[2, ] * 100
x_label <- paste0("PC1 (", round(explained_var[1], 1), "%)")
y_label <- paste0("PC2 (", round(explained_var[2], 1), "%)")

# Combine with metadata (assuming you have a metadata dataframe called metadata_df)
PCA_points <- cbind(PCA_scores, df_types2)  # metadata_df must have 'treatment', 'day', and 'cluster'

# Calculate convex hulls for each cluster
matched_clusters_PCA <- df_types[match(rownames(PCA_points), rownames(df_types)), ]
PCA_points$cluster <- matched_clusters_PCA  # Add clusters from df_types

hull_data <- PCA_points %>%
  group_by(cluster) %>%
  slice(chull(PC1, PC2))

pdf(file = paste0(figDir, "S3-B_PCA_treatment_day_polygon_daylabel_MAG.pdf"), width = 12, height = 12)
ggplot(PCA_points, aes(x = PC1, y = PC2, color = treatment, shape = factor(day))) +
  geom_polygon(data = hull_data, aes(x = PC1, y = PC2, group = cluster),
               fill = NA, color = "black", linetype = "solid", linewidth = 1) +  # Add hulls
  geom_point(size = 15, stroke = 2, alpha = 0.9) +  # Increase stroke for hollow shapes
  scale_color_manual(values = cols_df) +
  theme_minimal() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        legend.position = "none",  # Remove all legends
        panel.grid = element_blank(),  # Remove gridlines
        axis.line = element_line(color = "black", linewidth = 0.5)) +  # Add axis lines
  labs(x = x_label,
       y = y_label,
       title = "")

dev.off()


p_legend <- ggplot(PCA_points, aes(x = PC1, y = PC2, color = treatment, shape = factor(day))) +
  geom_point(size = 8, stroke = 2, alpha = 0.9) +
  scale_color_manual(values = cols_df) +
  labs(color = "Treatment", shape = "Day") +
  theme_minimal() +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        legend.key.size = unit(2, "lines"))

legend <- get_legend(p_legend)

pdf(file = paste0(figDir, "legend_treatment_day.pdf"), width = 4, height = 6)
grid::grid.draw(legend)
dev.off()

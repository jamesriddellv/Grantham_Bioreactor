library(vegan)
library(apcluster)
library(corrplot)
library(readxl)
library(dplyr)
library(ggplot2)
library(factoextra)
library(cluster)

#set location of your dataset and resulting data - change this for your own location
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

###################
###   Setup     ###
###################

figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"

# Define input files and their corresponding output names
input_configs <- list(
  list(
    file = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_vOTUs_wide.tsv',
    output_suffix = '1_read'
  ),
  list(
    file = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1_gene_per_10kb_vOTUs_wide.tsv',
    output_suffix = '1_gene_per_10kb'
  ),
  list(
    file = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_prop10_vOTUs_wide.tsv',
    output_suffix = 'prop10'
  )
)

# Load metadata (shared across all analyses)
metadata <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")
metadata <- metadata %>%
  mutate(treatment = case_when(
    treatment == "unamended" ~ "U",
    treatment == "CT" ~ "T",
    treatment == "catechin" ~ "C",
    TRUE ~ treatment
  )) %>% filter(treatment != 'T')

###################
###  Functions  ###
###################

# Function to perform the entire analysis pipeline
perform_pca_analysis <- function(input_file, output_suffix, metadata, figDir) {
  
  cat("\n=== Processing:", input_file, "===\n")
  
  # Read data
  df <- read.csv(input_file, sep='\t')
  
  # Prepare data
  treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')
  
  rownames(df) <- df$vOTU
  df$vOTU <- NULL
  colnames(df) <- treatment_timepoint
  
  # Make unique column names
  colnames(df) <- make.unique(names(df))
  
  # Transform data
  t_df <- t(df) %>% as_tibble(rownames = NA)
  log_df <- decostand(t_df, method = "hellinger")
  
  # Calculate Bray-Curtis distance
  bray_dist <- vegdist(log_df, method = "bray", binary = FALSE, diag = FALSE, 
                       upper = FALSE, na.rm = FALSE)
  bray_matrix_df <- as.matrix(bray_dist)
  
  # Affinity propagation clustering
  set.seed(42)
  df_ap_clust <- apcluster(negDistMat(r=2), bray_matrix_df, details = TRUE, seed=42)
  
  # Extract clusters
  df_clusters <- length(df_ap_clust)
  df_types <- data.frame()
  
  for (i in 1:df_clusters) {
    df_temp_types <- data.frame("site" = names(df_ap_clust[[i]]),
                                "type" = paste0("Cluster ", i))
    df_types <- bind_rows(df_types, df_temp_types)
  }
  
  # Format cluster data
  sites_df <- df_types$site
  df_types <- as.data.frame(df_types[, -1])
  rownames(df_types) <- sites_df
  names(df_types)[1] <- "type"
  
  # Wrangle metadata
  days <- sapply(strsplit(rownames(log_df), "day"), function(x) x[2])
  days <- sapply(strsplit(days, '\\.'), function(x) x[1])
  treatments <- sapply(strsplit(rownames(log_df), "_"), function(x) x[1])
  
  df_types2 <- as.data.frame(cbind(days, treatments))
  colnames(df_types2) <- c('day', 'treatment')
  rownames(df_types2) <- rownames(log_df)
  
  # Run PCA
  PCA_df <- prcomp(log_df, center = FALSE, scale. = FALSE)
  PCA_scores <- as.data.frame(PCA_df$x[, 1:2])
  colnames(PCA_scores) <- c("PC1", "PC2")
  explained_var <- summary(PCA_df)$importance[2, ] * 100
  x_label <- paste0("PC1 (", round(explained_var[1], 1), "%)")
  y_label <- paste0("PC2 (", round(explained_var[2], 1), "%)")
  
  # Combine with metadata
  PCA_points <- cbind(PCA_scores, df_types2)
  
  # Add clusters
  matched_clusters_PCA <- df_types[match(rownames(PCA_points), rownames(df_types)), ]
  PCA_points$cluster <- matched_clusters_PCA
  
  # Calculate convex hulls
  hull_data <- PCA_points %>%
    group_by(cluster) %>%
    slice(chull(PC1, PC2))
  
  # Define colors
  cols_df <- c("#FC9B2D", "#7ACBC3") # Catechin, Unamended
  names(cols_df) <- levels(factor(df_types2$treatment))
  
  # Generate plot
  output_file <- paste0(figDir, "S2-B_PCA_treatment_day_polygon_daylabel_", output_suffix, ".pdf")
  cat("Saving plot to:", output_file, "\n")
  
  pdf(file = output_file, width = 12, height = 12)
  p <- ggplot(PCA_points, aes(x = PC1, y = PC2, color = treatment, shape = factor(day))) +
    geom_polygon(data = hull_data, aes(x = PC1, y = PC2, group = cluster),
                 fill = NA, color = "black", linetype = "solid", linewidth = 1) +
    geom_point(size = 15, stroke = 2, alpha = 0.9) +
    scale_color_manual(values = cols_df) +
    theme_minimal() +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5)) +
    labs(x = x_label,
         y = y_label,
         title = "")
  print(p)
  dev.off()
  
  cat("Completed:", output_suffix, "\n")
  
  # Return results for optional further analysis
  return(list(
    PCA_scores = PCA_points,
    hull_data = hull_data,
    explained_var = explained_var
  ))
}

###################
###  Run Loop   ###
###################

# Loop through all configurations
results_list <- list()

for (config in input_configs) {
  results <- perform_pca_analysis(
    input_file = config$file,
    output_suffix = config$output_suffix,
    metadata = metadata,
    figDir = figDir
  )
  
  results_list[[config$output_suffix]] <- results
}

cat("\n=== All analyses complete! ===\n")
cat("Generated", length(input_configs), "PCA plots in:", figDir, "\n")
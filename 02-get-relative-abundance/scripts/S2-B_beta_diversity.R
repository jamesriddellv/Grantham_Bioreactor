library(vegan)
library(corrplot)
library(readxl)
library(dplyr)
library(ggplot2)

# Set location
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

###################
###   Setup     ###
###################

figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"

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

# Load metadata
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

perform_pca_analysis <- function(input_file, output_suffix, metadata, figDir) {
  
  cat("\n=== Processing:", input_file, "===\n")
  
  # Read data
  df <- read.csv(input_file, sep='\t')
  
  # Prepare data
  treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')
  rownames(df) <- df$vOTU
  df$vOTU <- NULL
  colnames(df) <- make.unique(treatment_timepoint)
  
  # Transform data
  t_df <- t(df) %>% as_tibble(rownames = NA)
  log_df <- decostand(t_df, method = "hellinger")
  
  # Wrangle metadata from rownames
  days <- sapply(strsplit(rownames(log_df), "day"), function(x) x[2])
  days <- sapply(strsplit(days, '\\.'), function(x) x[1])
  treats <- sapply(strsplit(rownames(log_df), "_"), function(x) x[1])
  
  meta_combined <- data.frame(
    day = factor(days, levels = c("0", "7", "14", "21", "35")),
    treatment = factor(treats),
    row.names = rownames(log_df)
  )
  
  # Run PCA
  PCA_df <- prcomp(log_df, center = FALSE, scale. = FALSE)
  PCA_scores <- as.data.frame(PCA_df$x[, 1:2])
  colnames(PCA_scores) <- c("PC1", "PC2")
  
  explained_var <- summary(PCA_df)$importance[2, ] * 100
  x_label <- paste0("PC1 (", round(explained_var[1], 1), "%)")
  y_label <- paste0("PC2 (", round(explained_var[2], 1), "%)")
  
  PCA_points <- cbind(PCA_scores, meta_combined)
  
  # Aesthetic Mappings
  cols_df <- c("C" = "#FC9B2D", "U" = "#7ACBC3")
  pch_df  <- c('0'=19, '7'=7, '14'=17, '21'=15, '35'=3)
  
  # Generate plot
  output_file <- paste0(figDir, "S2-B_PCA_", output_suffix, ".pdf")
  cat("Saving plot to:", output_file, "\n")
  
  pdf(file = output_file, width = 11, height = 10)
  p <- ggplot(PCA_points, aes(x = PC1, y = PC2, color = treatment, shape = day)) +
    geom_point(size = 12, stroke = 2.5, alpha = 0.85) +
    scale_color_manual(values = cols_df) +
    scale_shape_manual(values = pch_df) +
    theme_minimal() +
    theme(axis.text = element_text(size = 22, color = "black"),
          axis.title = element_text(size = 24, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          panel.grid.major = element_line(color = "grey92"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.8)) +
    labs(x = x_label, y = y_label, color = "Treatment", shape = "Day")
  
  print(p)
  dev.off()
  
  return(list(PCA_scores = PCA_points, explained_var = explained_var))
}

###################
###   Run Loop  ###
###################

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
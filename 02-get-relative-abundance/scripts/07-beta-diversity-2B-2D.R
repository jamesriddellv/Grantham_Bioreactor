library(vegan)

library(dplyr)

library(ggplot2)



# ── 0. Load and Prepare Data ──────────────────────────────────────────────────



setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"



# Load counts and metadata

df <- read.csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')

metadata <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")



# Clean metadata

metadata <- metadata %>%
  
  mutate(treatment = case_when(
    
    treatment == "unamended" ~ "U",
    
    treatment == "CT"        ~ "T",
    
    treatment == "catechin"  ~ "C",
    
    TRUE                     ~ treatment
    
  )) %>% 
  
  filter(treatment != 'T')



# Prepare rownames and column names for the count matrix

rownames(df) <- df$vOTU

df$vOTU <- NULL

treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')

colnames(df) <- make.unique(treatment_timepoint)



# Transpose for vegan compatibility

t_df <- as.data.frame(t(df))



# ── 1. Transformation and Distance ───────────────────────────────────────────



hel_df    <- decostand(t_df, method = "hellinger")

bray_dist <- vegdist(hel_df, method = "bray")



# ── 2. Define Formal Clusters ────────────────────────────────────────────────



# Create a clean sample info dataframe based on your specific cluster logic

sample_info <- data.frame(row_names = rownames(hel_df)) %>%
  
  mutate(
    
    treatment = sub("_.*", "", row_names),
    
    day = as.numeric(regmatches(row_names, regexpr("\\d+", row_names))),
    
    formal_cluster = case_when(
      
      (treatment == "C" & day == 7) | (treatment == "U" & day == 7) ~ "Cluster 1 (Day 7 pre-response)",
      
      treatment == "C" & day %in% c(14, 21, 35) ~ "Cluster 2 (C response)",
      
      treatment == "U" & day %in% c(14, 21, 35) ~ "Cluster 3 (U response)",
      
      
      
      TRUE ~ NA_character_
      
    )
    
  )

# Plot Figure S2 heatmap.
library(pheatmap)
library(RColorBrewer)

# 1. Prepare and Sort Metadata for Ordering
# We ensure 'day' is a factor or sorted numeric so the heatmap follows chronological order
plot_info <- sample_info %>%
  arrange(treatment, day) %>%
  # Convert day to factor so it's treated as a discrete group in the legend
  mutate(day = as.factor(day)) 

# 2. Convert Distance Object to Matrix and Reorder
# We must reorder the distance matrix to match our sorted metadata
bray_mat <- as.matrix(bray_dist)
bray_mat <- bray_mat[plot_info$row_names, plot_info$row_names]

write.csv(bray_mat, "/users/PAS1573/riddell26/data/S3_bray_matrix.csv", row.names = TRUE)

# 3. Create Annotation Dataframe for the Heatmap
# This tells pheatmap which samples belong to which treatment/day
ann_df <- plot_info %>%
  select(treatment, day) %>%
  as.data.frame()
rownames(ann_df) <- plot_info$row_names
# 4. Define Colors Dynamically (Crash-proof)
# Find all unique days in your actual data
unique_days <- levels(plot_info$day)

# Generate exactly enough blue shades for the number of unique days
day_colors <- colorRampPalette(brewer.pal(9, "Blues"))(length(unique_days))
# Name the colors so pheatmap knows which color goes to which day
names(day_colors) <- unique_days 

ann_colors = list(
  treatment = c(U = "#7ACBC3", C = "#FC9B2D"),
  day = day_colors
)

# 5. Plot the Heatmap
heatmap_plot <- pheatmap(
  bray_mat,
  # Update this line to include a middle color stop
  color = colorRampPalette(c("white", "gold", "firebrick3"))(100), 
  annotation_col = ann_df,       
  annotation_row = ann_df,       
  annotation_colors = ann_colors, 
  cluster_rows = FALSE,          
  cluster_cols = FALSE,          
  main = "Bray-Curtis Dissimilarity",
  display_numbers = FALSE,
  fontsize = 10
)

# Save the plot
ggsave(filename = paste0(figDir, "SX_vOTU_bray_heatmap.pdf"), plot = heatmap_plot, width = 10, height = 8)


# ── 3. Statistical Testing ────────────────────────────────────────────────────



# Filtering for the formal clusters only

test_idx <- !is.na(sample_info$formal_cluster)

bray_dist_filt <- as.dist(as.matrix(bray_dist)[test_idx, test_idx])

metadata_filt  <- sample_info[test_idx, ]



# A. PERMDISP (Variance Check)

disp_mod  <- betadisper(bray_dist_filt, metadata_filt$formal_cluster)

disp_test <- permutest(disp_mod, permutations = 999)



# B. Global PERMANOVA (Overall Difference)

global_perm <- adonis2(bray_dist_filt ~ formal_cluster, 
                       
                       data = metadata_filt, 
                       
                       permutations = 999)



# C. Pairwise PERMANOVA (1-on-1 Comparisons)

# We manually loop through the 3 clusters to compare them

clusters_to_test <- unique(metadata_filt$formal_cluster)

cluster_pairs <- combn(clusters_to_test, 2)

pairwise_results <- data.frame(Comparison = character(), P_Value = numeric(), R2 = numeric())



for(i in 1:ncol(cluster_pairs)) {
  
  pair <- cluster_pairs[,i]
  
  
  
  # Subset to only the two clusters in the current pair
  
  pair_idx <- metadata_filt$formal_cluster %in% pair
  
  pair_dist <- as.dist(as.matrix(bray_dist_filt)[pair_idx, pair_idx])
  
  pair_data <- metadata_filt[pair_idx, ]
  
  
  
  # Run PERMANOVA for this pair
  
  pair_res <- adonis2(pair_dist ~ formal_cluster, data = pair_data, permutations = 999)
  
  
  
  # Store result
  
  pairwise_results <- rbind(pairwise_results, data.frame(
    
    Comparison = paste(pair, collapse = " vs "),
    
    P_Value = pair_res$`Pr(>F)`[1],
    
    R2 = pair_res$R2[1]
    
  ))
  
}



# Apply p-value correction (Bonferroni) to handle multiple comparisons

pairwise_results$P_Adj <- p.adjust(pairwise_results$P_Value, method = "bonferroni")



print("--- Global Results ---")

print(global_perm)

print("--- Pairwise Results (with Bonferroni Correction) ---")

print(pairwise_results)



# --- Pairwise Non-Parametric Comparisons ---



# Extract distances and groups from your dispersion model

dist_to_centroid <- disp_mod$distances

groups <- disp_mod$group



# Perform Pairwise Wilcoxon tests (Non-parametric)

# This specifically compares medians of distances between 1&2, 2&3, and 1&3

pairwise_results <- pairwise.wilcox.test(dist_to_centroid, 
                                         
                                         groups, 
                                         
                                         p.adjust.method = "fdr")



print("--- Pairwise Dispersion Comparisons (Wilcoxon with FDR) ---")

print(pairwise_results)

# ── 3.E ANOSIM (C1 vs C3) ───────────────────────────────────────────────────

# 1. Identify the clusters we want to compare
target_clusters <- c("Cluster 1 (Day 7 pre-response)", "Cluster 3 (U response)")

# 2. Filter the distance matrix and metadata
anosim_idx   <- metadata_filt$formal_cluster %in% target_clusters
# Ensure the subset distance matrix is converted back to 'dist' class
anosim_dist  <- as.dist(as.matrix(bray_dist_filt)[anosim_idx, anosim_idx])
anosim_meta  <- metadata_filt[anosim_idx, ]

# 3. Run ANOSIM
# Use 999 permutations for a robust p-value
anosim_res <- anosim(anosim_dist, anosim_meta$formal_cluster, permutations = 999)

# 4. Print Results
print("--- ANOSIM Results: Pre-response vs. Unamended ---")
print(anosim_res)

# Extract specific values for your records
anosim_R <- anosim_res$statistic
anosim_p <- anosim_res$signif

cat("\nANOSIM R-statistic:", anosim_R)
cat("\nANOSIM p-value:", anosim_p, "\n")

# --- Visualization: Dispersion per Cluster ---

pdf(file = paste(figDir, "SX_cluster_dispersion.pdf", sep=''), width = 10, height = 10);

# Set up margins for better label visibility

par(mar = c(5, 5, 4, 2))



# 1. Create the boxplot with your specified colors

boxplot(dist_to_centroid ~ groups, 
        
        main = "Dispersion per Cluster", 
        
        ylab = "Distance to Centroid", 
        
        xlab = "Cluster",
        
        col = c("#FC9B2D", "#7ACBC3", "grey"), # Your requested colors
        
        border = "black",
        
        outline = FALSE, # Hide default outliers to avoid double-plotting with jitter
        
        las = 1,
        
        frame.plot = FALSE)



# 2. Add individual points (Jitter) to show the actual samples

# This makes the n=3 Day 0 vs. High-variation Day 35 very clear

stripchart(dist_to_centroid ~ groups, 
           
           vertical = TRUE, 
           
           add = TRUE, 
           
           pch = 21, 
           
           bg = "white", # Points with white fill to pop against colors
           
           col = "black", 
           
           cex = 1.1)

dev.off()

# ── 4. PCA Projection ────────────────────────────────────────────────────────



pca_res <- prcomp(hel_df, center = FALSE, scale. = FALSE)

pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

x_label <- paste0("PC1 (", pct_var[1], "%)")

y_label <- paste0("PC2 (", pct_var[2], "%)")



# Merge PCA coordinates with metadata

PCA_points <- data.frame(
  
  PC1 = pca_res$x[, 1],
  
  PC2 = pca_res$x[, 2],
  
  treatment = sample_info$treatment,
  
  day = sample_info$day,
  
  formal_cluster = sample_info$formal_cluster
  
)



# ── 5. Statistical Annotation Label ──────────────────────────────────────────



stats_label <- paste0(
  
  "PERMANOVA: p = ", global_perm$`Pr(>F)`[1], "\n",
  
  "PERMDISP:  p = ", disp_test$tab$`Pr(>F)`[1]
  
)



# ── 6. Final Plotting (Modified) ──────────────────────────────────────────────



cols_df <- c("C" = "#FC9B2D", "U" = "#7ACBC3") 

pch_df  <- c('0'=19, '7'=7, '14'=17, '21'=15, '35'=3)



# Reduced PDF size to 8x8 as a 10x10 square plot is quite large

pdf(file = paste0(figDir, "02-B_PCA_Cluster_CI_Treatment_Colors.pdf"), width = 8, height = 8)



ggplot(PCA_points, aes(x = PC1, y = PC2)) +
  
  # 95% CI Ellipses
  
  stat_ellipse(data = subset(PCA_points, !is.na(formal_cluster)),
               
               aes(group = formal_cluster, linetype = formal_cluster),
               
               color = "black", 
               
               level = 0.95, 
               
               linewidth = 0.8) +
  
  
  
  # Points
  
  geom_point(aes(color = treatment, shape = factor(day)), 
             
             size = 10, stroke = 1.5, alpha = 0.8) +
  
  
  
  # Styling & Legend Removal (guide = "none")
  
  scale_color_manual(values = cols_df, guide = "none") +
  
  scale_shape_manual(values = pch_df, guide = "none") +
  
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = "none") +
  
  
  
  # Force Square Aspect Ratio
  
  coord_fixed(ratio = 1) + 
  
  
  
  theme_minimal() +
  
  theme(
    
    aspect.ratio = 1, # Double reinforcement for square plot area
    
    axis.text    = element_text(size = 20),
    
    axis.title   = element_text(size = 22),
    
    legend.position = "none", # Removes all legend boxes
    
    panel.grid   = element_blank(),
    
    axis.line    = element_line(color = "black")
    
  )



dev.off()

# Base R approach (includes row names if your sample IDs are row names)
write.csv(PCA_points, file = "/users/PAS1573/riddell26/data/2B_PCA_points.csv", row.names = TRUE)

####################################
#    Indicator Species Analysis    #
####################################

library(indicspecies)

set.seed(42)

# 1. Define the 5 specific groups
isa_metadata <- sample_info %>%
  mutate(isa_group = case_when(
    treatment == "U" & day == 0  ~ "Group1_Unam_D0",
    day == 7                     ~ "Group2_D7_All", 
    treatment == "U" & day > 7   ~ "Group3_Unam_Late",
    treatment == "C" & day == 21 ~ "Group4_Cat_D21",       # New separate cluster
    treatment == "C" & day %in% c(14, 35) ~ "Group5_Cat_D14_D35", # Combined remaining Cat late
    TRUE                         ~ NA_character_
  )) %>%
  filter(!is.na(isa_group)) 

# 2. Align the abundance data (t_df) with our new groups
df_isa_input <- t_df[isa_metadata$row_names, ]

# --- ISA with duleg = TRUE (Single group indicators) ---

isa_duleg <- multipatt(
  x = df_isa_input,
  cluster = isa_metadata$isa_group,
  duleg = TRUE
)

summary(isa_duleg, alpha = 0.01)

assoc_matrix <- isa_duleg$sign
assoc_matrix <- assoc_matrix[!is.na(assoc_matrix$p.value),]
alpha <- 0.01

all_clusters_list <- list()

for (col_name in colnames(assoc_matrix)) {
  if (grepl("^s\\.", col_name)) {  
    group_name <- gsub("^s\\.", "", col_name)
    
    is_sig <- assoc_matrix[, col_name] == 1 & assoc_matrix$p.value <= alpha
    sig_vOTUs <- rownames(assoc_matrix)[is_sig]
    sig_pvals <- assoc_matrix$p.value[is_sig]
    
    if(length(sig_vOTUs) > 0) {
      all_clusters_list[[group_name]] <- data.frame(
        vOTU = sig_vOTUs,
        P_value = sig_pvals,
        cluster = group_name
      )
    }
  }
}

final_duleg_table <- do.call(rbind, all_clusters_list)
write.csv(final_duleg_table, '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_duleg.csv', row.names = FALSE)

# --- ISA with duleg = FALSE (Multiple group indicators) ---

isa_multi <- multipatt(
  x = df_isa_input,
  cluster = isa_metadata$isa_group,
  duleg = FALSE
)

assoc_multi <- isa_multi$sign
assoc_multi <- assoc_multi[!is.na(assoc_multi$p.value) & assoc_multi$p.value <= 0.01, ]

# Map indices to the 5 group names
group_levels <- levels(as.factor(isa_metadata$isa_group))
cluster_cols <- assoc_multi[, 1:length(group_levels)]

multi_results <- apply(cluster_cols, 1, function(x) {
  matches <- which(x == 1)
  paste(group_levels[matches], collapse = "+")
})

final_multi_table <- data.frame(
  vOTU = rownames(assoc_multi),
  significant_groups = unname(multi_results),
  p_value = assoc_multi$p.value,
  stringsAsFactors = FALSE
)

write.csv(final_multi_table, '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv', row.names = FALSE)
write.csv(final_multi_table, '/users/PAS1573/riddell26/data/2D_indicator_vOTUs_per_cluster_multi.csv', row.names = FALSE)

print("ISA Analysis Complete with 5 clusters. Files saved to results directory.")


# ── UpSet Plot Visualization for ISA Multi-grouping ──────────────────────────

# Load necessary libraries
library(ComplexUpset)
library(ggplot2)

# 1. Prepare data for UpSet plot
# We use the 'final_multi_table' generated from the ISA multi analysis
# It contains vOTU and significant_groups (e.g., "Group1+Group2")

# Define the five primary groups
groups <- c("Group1_Unam_D0", "Group2_D7_All", "Group3_Unam_Late", 
            "Group4_Cat_D21", "Group5_Cat_D14_D35")

# Create a binary matrix (0/1) for intersections
upset_df <- final_multi_table
for(g in groups) {
  upset_df[[g]] <- as.integer(grepl(g, upset_df$significant_groups))
}

# 2. Define custom colors based on the mapping provided
# These match the hex codes from your reference image
group_colors <- c(
  'Group2_D7_All'                                    = '#C994C7',
  'Group2_D7_All+Group3_Unam_Late'                   = '#6A51A3',
  'Group2_D7_All+Group3_Unam_Late+Group4_Cat_D21'    = '#737373', # Grey from image
  'Group2_D7_All+Group4_Cat_D21+Group5_Cat_D14_D35'  = '#FDBB84',
  'Group2_D7_All+Group5_Cat_D14_D35'                 = '#7BCCC4',
  'Group3_Unam_Late'                                 = '#08589E',
  'Group4_Cat_D21'                                   = '#DD0597',
  'Group4_Cat_D21+Group5_Cat_D14_D35'                = '#238B45'
)

# 3. Generate the UpSet Plot
pdf(file = paste0(figDir, "02-D_ISA_UpSet_Plot.pdf"), width = 10, height = 7)

upset(
  upset_df,
  groups,
  width_ratio = 0.1,
  min_size = 1,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = significant_groups)
    ) + 
      scale_fill_manual(values = group_colors, guide = "none") +
      ylab('vOTUs') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ),
  matrix = intersection_matrix(
    geom = geom_point(size = 3)
  ),
  set_sizes = FALSE, # Focus on intersection sizes
  sort_intersections = 'descending'
) +
  ggtitle("Significant Indicator vOTUs per Group Combination (p <= 0.01)")

dev.off()

print("UpSet Plot saved to figures directory.")


###########################################
## Bray curtis distances between samples ##
###########################################
# 1. Convert distance object to a matrix
bray_mat <- as.matrix(bray_dist)

# --- Comparison 1: Day 0 vs Day 7 (Irrespective of Treatment) ---

# Get indices based on the 'day' column in sample_info
idx_d0 <- which(sample_info$day == 0)
idx_d7 <- which(sample_info$day == 7)

# Extract distances between D0 and D7
# Note: If Day 0 was filtered out earlier, this vector will be empty
dist_0_vs_7 <- bray_mat[idx_d0, idx_d7]

avg_0_7 <- mean(dist_0_vs_7)
sd_0_7  <- sd(dist_0_vs_7)

# --- Comparison 2: Day 7 Catechin (C) vs Day 7 Unamended (U) ---

# Get indices for specific treatments at Day 7
idx_d7_c <- which(sample_info$day == 7 & sample_info$treatment == "C")
idx_d7_u <- which(sample_info$day == 7 & sample_info$treatment == "U")

# Extract distances between C_7 and U_7
dist_7_c_vs_u <- bray_mat[idx_d7_c, idx_d7_u]

avg_7_treats <- mean(dist_7_c_vs_u)
sd_7_treats  <- sd(dist_7_c_vs_u)

# --- Print Formatted Results for your Paper ---

cat("--- Results for Manuscript Section --- \n\n")

cat(paste0("Initial Shift (Day 0 to 7): \n", 
           "Avg Bray-Curtis = ", round(avg_0_7, 3), 
           " (SD = ", round(sd_0_7, 3), ")\n\n"))

cat(paste0("Treatment Difference at Day 7 (C vs U): \n", 
           "Avg Bray-Curtis = ", round(avg_7_treats, 3), 
           " (SD = ", round(sd_7_treats, 3), ")\n"))

###########################################
# Check Day 7 unamended vs Day 7 catechin #
###########################################

# 1. Subset the data for Day 7 only (Unamended vs Catechin)
# ---------------------------------------------------------
day7_idx <- sample_info$day == 7 & sample_info$treatment %in% c("U", "C")

# Ensure we have data for both groups
day7_metadata <- sample_info[day7_idx, ]
day7_dist     <- as.dist(as.matrix(bray_dist)[day7_idx, day7_idx])

# 2. PERMANOVA (Testing Centroid Position/Composition)
# ---------------------------------------------------------
day7_permanova <- adonis2(day7_dist ~ treatment, 
                          data = day7_metadata, 
                          permutations = 999)

print("=== PERMANOVA: Day 7 Unamended vs Catechin ===")
print(day7_permanova)

# 3. Beta-Dispersion with Wilcoxon (Testing Variance/Spread)
# ---------------------------------------------------------
# Step A: Calculate distances to centroids
disp_mod_day7 <- betadisper(day7_dist, day7_metadata$treatment)

# Step B: Extract distances and run Wilcoxon Rank-Sum test
dist_values <- disp_mod_day7$distances
group_labels <- day7_metadata$treatment

wilcox_disp <- wilcox.test(dist_values ~ group_labels)

print("=== Wilcoxon Test for Beta-Dispersion (Day 7) ===")
print(wilcox_disp)

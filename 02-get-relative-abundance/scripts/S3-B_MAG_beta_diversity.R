library(vegan)
library(apcluster)
library(dplyr)
library(ggplot2)
library(factoextra)
library(cluster)
library(readxl)



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

#################################################
# Test PERMANOVA catechin vs unamended day 7-35 #
#################################################

# ── 1. Prepare Metadata for Statistics ───────────────────────────────────────
# We need a metadata dataframe that matches the rows of 'log_df'
sample_info <- data.frame(row_names = rownames(log_df)) %>%
  mutate(
    treatment = sub("_.*", "", row_names),
    # Extract numeric day and handle suffixes like .1, .2
    day_num = as.numeric(regmatches(row_names, regexpr("\\d+", row_names)))
  )

# ── 1. Define and Filter Clusters ───────────────────────────────────────────

# Ensure formal_cluster is correctly mapped in the sample_info
cluster_data <- sample_info %>%
  mutate(
    # Create the three specific groups you requested
    test_group = case_when(
      day_num == 7                                ~ "Group1_Day7_Baseline",
      treatment == "U" & day_num %in% c(14,21,35) ~ "Group2_Unamended_Response",
      treatment == "C" & day_num %in% c(14,21,35) ~ "Group3_Catechin_Response",
      TRUE                                        ~ NA_character_
    )
  ) %>%
  # Remove Day 0 or any samples not in our 3 groups
  filter(!is.na(test_group))

# ── 2. Sync Distance Matrix with Filtered Metadata ──────────────────────────

# Subset the transformed data to match our filtered rows
log_df_filtered <- log_df[cluster_data$row_names, ]

# Calculate Bray-Curtis distance on the filtered subset
dist_filtered <- vegdist(log_df_filtered, method = "bray")

# ── 3. Run PERMANOVA on the Three Clusters ──────────────────────────────────

# We test 'test_group' which contains our 3 defined states
permanova_clusters <- adonis2(dist_filtered ~ test_group, 
                              data = cluster_data, 
                              permutations = 999)

print("=== PERMANOVA: Testing the 3 Biological Clusters ===")
print(permanova_clusters)

# ── 4. Pairwise Comparison (Post-hoc) ───────────────────────────────────────
# Since PERMANOVA only tells you if *any* group is different, 
# you might want to know if Group 2 is specifically different from Group 3.
# If you have the 'pairwise.adonis2' or 'RVAideMemoire' library:
# ── Pairwise PERMANOVA using adonis2 ──────────────────────────────────────────

# Function to run pairwise adonis2
run_pairwise_adonis <- function(dist_obj, groups) {
  groups <- as.factor(groups)
  levels <- levels(groups)
  pairs <- combn(levels, 2)
  results <- data.frame()
  
  for(i in 1:ncol(pairs)) {
    # Get labels for the pair
    g1 <- pairs[1, i]
    g2 <- pairs[2, i]
    
    # Identify samples belonging to this pair
    idx <- which(groups %in% c(g1, g2))
    sub_dist <- as.matrix(dist_obj)[idx, idx]
    sub_meta <- groups[idx]
    
    # Run adonis2
    ad <- adonis2(as.dist(sub_dist) ~ sub_meta, permutations = 999)
    
    # Extract R2 and p-value
    res <- data.frame(
      comparison = paste(g1, "vs", g2),
      R2 = ad$R2[1],
      p_value = ad$`Pr(>F)`[1]
    )
    results <- rbind(results, res)
  }
  
  # Adjust p-values for multiple comparisons (FDR)
  results$p_adj <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

# Run the test
cat("\n=== Pairwise adonis2 Results (3 Clusters) ===\n")
pairwise_results <- run_pairwise_adonis(dist_filtered, cluster_data$test_group)

print(pairwise_results)

# ── Interpretation ────────────────────────────────────────────────────────────
# If p_adj < 0.05, the two clusters are significantly distinct.
# R2 tells you the effect size (e.g., 0.40 = 40% of variation is due to group).

# ── 5. Check Dispersion (Spread) ────────────────────────────────────────────
# If dispersion is significant, differences might be due to variance, not just mean.
disp_clusters <- betadisper(dist_filtered, cluster_data$test_group)
print("=== Beta-Dispersion Test for 3 Clusters ===")
print(anova(disp_clusters))

# ──────────────────────────────────────────────────────────────────────────────
# TEST 2: Pre-Response Check (Day 7 ONLY)
# ──────────────────────────────────────────────────────────────────────────────

# 1. Prepare Subset
target_7_only <- which(sample_info$day_num == 7)
dist_7_only   <- vegdist(log_df[target_7_only, ], method = "bray")
meta_7_only   <- sample_info[target_7_only, ]

# 2. Run PERMANOVA 
# FIX: Remove 'day_num' from formula since it is now a constant (all samples = 7)
perm_7_only <- adonis2(dist_7_only ~ treatment, 
                       data = meta_7_only, 
                       permutations = 999)

print("=== PERMANOVA Results (Day 7 Only) ===")
print(perm_7_only)

# 3. Beta-Dispersion 
disp_7_only <- betadisper(dist_7_only, meta_7_only$treatment)
print("=== Beta-Dispersion (Day 7 Only) ===")
print(anova(disp_7_only))

# cols_df and pch_df from your original code
cols_df <- c("#FC9B2D", "#7ACBC3")
names(cols_df) <- c("C", "U")
pch_df <- c(4, 15, 17, 18, 19)
names(pch_df) <- c("0", "7", "14", "21", "35")

# ── 1. Calculate PCA ─────────────────────────────────────────────────────────

# Perform PCA on the Hellinger-transformed data
# we use scale = FALSE because Hellinger transformation already handles scaling
pca_res <- prcomp(log_df, center = FALSE, scale. = FALSE)

# Extract coordinates for plotting
PCA_points <- as.data.frame(pca_res$x)
PCA_points$row_names <- rownames(PCA_points)

# Calculate Percentage of Variance for Axis Labels
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# ── 2. Merge Metadata for Plotting ──────────────────────────────────────────

# Join with the sample_info created earlier to get treatment and day
PCA_points <- PCA_points %>%
  mutate(
    treatment = sub("_.*", "", row_names),
    day = as.character(regmatches(row_names, regexpr("\\d+", row_names))),
    # Define the same formal clusters used in your previous logic for ellipses
    formal_cluster = case_when(
      day == "7" ~ "Cluster 1 (Day 7)",
      treatment == "C" & day %in% c("14", "21", "35") ~ "Cluster 2 (C response)",
      treatment == "U" & day %in% c("14", "21", "35") ~ "Cluster 3 (U response)",
      TRUE ~ NA_character_
    )
  )

# ── 3. Final Plotting (Aesthetic Match) ──────────────────────────────────────

cols_df <- c("C" = "#FC9B2D", "U" = "#7ACBC3") 
pch_df  <- c('0'=19, '7'=7, '14'=17, '21'=15, '35'=3)

pdf(file = paste0(figDir, "S4-B_MAG_beta_diversity.pdf"), width = 8, height = 8)

ggplot(PCA_points, aes(x = PC1, y = PC2)) +
  # 95% CI Ellipses
  stat_ellipse(data = subset(PCA_points, !is.na(formal_cluster)),
               aes(group = formal_cluster, linetype = formal_cluster),
               color = "black", 
               level = 0.95, 
               linewidth = 0.8) +
  
  # Points: Size 10 and Stroke 1.5 match your requested style
  geom_point(aes(color = treatment, shape = factor(day)), 
             size = 10, stroke = 1.5, alpha = 0.8) +
  
  # Scaling
  scale_color_manual(values = cols_df, guide = "none") +
  scale_shape_manual(values = pch_df, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = "none") +
  
  # Labels and Force Square Aspect Ratio
  labs(x = pc1_label, y = pc2_label) +
  coord_fixed(ratio = 1) + 
  
  theme_minimal() +
  theme(
    aspect.ratio = 1, 
    axis.text    = element_text(size = 20),
    axis.title   = element_text(size = 22),
    legend.position = "none", 
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black")
  )

dev.off()
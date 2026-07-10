library(vegan)
library(apcluster)
library(dplyr)
library(ggplot2)
library(factoextra)
library(cluster)
# library(fpc)       # for Calinski-Harabasz

# ── 0. Load data ──────────────────────────────────────────────────────────────

#set location of your dataset and resulting data - change this for your own location
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

###################
###   Default   ###
###################

figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"

#read in and format abundance data
# >5 reads map to lytic gene from 001-identify-active-virus-structural-genes 
df <- read.csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')

# >1 gene / 10kb with a read mapped
# df <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/active_vOTUs_1_gene_per_10kb_relative_abundance.tsv", sep='\t')

# >1 read mapped
# df <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/active_vOTUs_1_read_mapped_relative_abundance.tsv", sep='\t')


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
t_df <- t(df) %>% as_tibble(rownames = NA)
# ── 1. Compare transformations for normality ──────────────────────────────────
trans <- list(
  raw       = t_df,
  hellinger = decostand(t_df, method = "hellinger"),
  log1p     = log1p(t_df)
)

# Shapiro-Wilk on column-wise means (one value per sample) as a proxy
norm_summary <- sapply(trans, function(x) {
  col_means <- rowMeans(x)
  sw <- shapiro.test(col_means)
  c(W = round(sw$statistic, 4), p = round(sw$p.value, 4))
})
print("Normality of row means per transformation:")
print(norm_summary)

# Visual: density plots of all values
trans_long <- bind_rows(lapply(names(trans), function(nm) {
  data.frame(value = as.vector(as.matrix(trans[[nm]])), 
             transformation = nm)
}))

ggplot(trans_long %>% filter(value > 0), aes(x = value, color = transformation)) +
  geom_density() +
  facet_wrap(~transformation, scales = "free") +
  labs(title = "Distribution comparison: raw vs hellinger vs log(n+1)") +
  theme_minimal()

# ── 2. Bray-Curtis distance (Hellinger-transformed) ───────────────────────────
log_df   <- trans$hellinger          # swap to log1p if normality test favours it
bray_dist <- vegdist(log_df, method = "bray")
bray_mat  <- as.matrix(bray_dist)

# Convert to data frame for easier manipulation
bray_df <- as.data.frame(as.table(bray_mat))
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

################################################################
# Plot PCA with PERMANOVA comparing cluster catchin day 14-35 #
# against cluster catechin day 7 + unamended day 7-35         #
################################################################

# ── 1. Create the Comparison Grouping (Base R) ───────────────────────────────

# Extract Treatment (everything before the first underscore)
treatments <- sub("_.*", "", rownames(log_df))

# Extract Day (grabs the first number found after the underscore)
# This handles "C_0.1", "C_35", and "C_35.1" correctly
day_matches <- regmatches(rownames(log_df), regexpr("\\d+\\.?\\d*", rownames(log_df)))
days <- as.numeric(day_matches)

sample_info <- data.frame(
  row_names = rownames(log_df),
  treatment = treatments,
  day = days
) %>%
  mutate(
    # Define the target group
    comparison_group = ifelse(treatment == "C" & day >= 14, "Catechin_Late", "Other")
  )

# ──PERMANOVA & Beta Dispersion ────────────────────────────────────────────

# PERMANOVA: Tests if the 'average' community differs
perm_target <- adonis2(bray_dist ~ comparison_group, 
                       data = sample_info, 
                       permutations = 999)

# Beta Dispersion: Tests if the 'variance/spread' differs
# A non-significant p-value here is actually GOOD—it means your 
# PERMANOVA result isn't just an artifact of unequal variance.
disp_target <- betadisper(bray_dist, sample_info$comparison_group)
disp_anova  <- anova(disp_target)

cat("\n── PERMANOVA RESULTS ──\n")
print(perm_target)

cat("\n── BETA DISPERSION RESULTS ──\n")
print(disp_anova)

# ANOSIM
# Group 1: Catechin Day 14-35 | Group 2: Everything Else
sample_info$is_target <- ifelse(sample_info$treatment == "C" & sample_info$day >= 14, 
                                "Catechin_Late", "Others")

# ── 2. ANOSIM (The "Tightness" Test) ──────────────────────────────────────
# R-value near 1 suggests the groups are completely distinct.
# R-value near 0 suggests no difference in grouping.
target_anosim <- anosim(bray_dist, sample_info$is_target, permutations = 999)

cat("\n── ANOSIM Results ──\n")
print(target_anosim)


# Compare for Unamended day 14-35 vs Catechin day 14-35
# Define 4 groups: Unamended Early, Unamended Late, Catechin Early, Catechin Late
sample_info <- sample_info %>%
  mutate(
    time_period = ifelse(day >= 14, "Late", "Early"),
    specific_group = paste(treatment, time_period, sep = "_")
  )

# Verify the new groups
print(table(sample_info$specific_group))

# Filter distance matrix to ONLY include "Late" samples for a direct comparison
late_indices <- which(sample_info$time_period == "Late")
late_dist    <- as.dist(as.matrix(bray_dist)[late_indices, late_indices])
late_meta    <- sample_info[late_indices, ]

# 2. PERMANOVA: Catechin Late vs Unamended Late
perm_late <- adonis2(late_dist ~ treatment, 
                     data = late_meta, 
                     permutations = 999)

# 3. Beta Dispersion: Catechin Late vs Unamended Late
disp_late <- betadisper(late_dist, late_meta$treatment)
disp_late_pval <- permutest(disp_late, permutations = 999)$tab$`Pr(>F)`[1]

cat("\n── Comparison: Catechin Late vs. Unamended Late ──\n")
print(perm_late)
cat(sprintf("Beta Dispersion p-value: %.4f\n", disp_late_pval))

# ── 4. Visualization ──────────────────────────────────────────────────────────

# PCoA Plot specifically highlighting this comparison
ord_res <- run_ordination(log_df, bray_dist, method = "PCoA")
plot_df <- cbind(ord_res$scores, sample_info)

ggplot(plot_df, aes(x = Axis1, y = Axis2, color = comparison_group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Catechin_Late" = "#D55E00", "Other" = "#56B4E9")) +
  labs(title = "PCoA: Catechin (Days 14-35) vs. Rest of Samples",
       subtitle = paste("PERMANOVA p =", perm_target$`Pr(>F)`[1]),
       x = ord_res$axis_labels[1], 
       y = ord_res$axis_labels[2]) +
  theme_minimal()








# ── 3. Optimal cluster number ─────────────────────────────────────────────────

# (a) Silhouette
sil_plot <- fviz_nbclust(bray_mat, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
print(sil_plot)

# (b) Elbow / WSS
wss_plot <- fviz_nbclust(bray_mat, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
print(wss_plot)

# (c) Gap statistic
set.seed(42)
gap_stat <- clusGap(bray_mat, FUN = kmeans, nstart = 25, K.max = 10, B = 500)
print(fviz_gap_stat(gap_stat))

# (d) Affinity propagation (biologically motivated baseline = 5 clusters)
ap_clust <- apcluster(negDistMat(r = 2), bray_mat, details = TRUE, seed = 42)
cat("Affinity propagation found", length(ap_clust@clusters), "clusters\n")

# Helper: extract AP labels at a specified K
ap_labels <- function(k) {
  ap_k <- apclusterK(negDistMat(r = 2), bray_mat, K = k, details = TRUE, seed = 42)
  lab  <- rep(NA_character_, nrow(bray_mat))
  for (i in seq_along(ap_k@clusters))
    lab[ap_k@clusters[[i]]] <- paste0("Cluster_", i)
  setNames(lab, rownames(bray_mat))
}

# ── 4 & 5. Evaluate K = 4, 5, 6 ──────────────────────────────────────────────
eval_k <- function(k) {
  labels <- ap_labels(k)
  
  # PERMANOVA
  perm   <- adonis2(bray_dist ~ labels, permutations = 999)
  R2     <- perm$R2[1]
  pval   <- perm$`Pr(>F)`[1]
  
  # Beta dispersion
  disp   <- betadisper(bray_dist, labels)
  disp_p <- anova(disp)$`Pr(>F)`[1]
  
  # Calinski-Harabasz
  km     <- kmeans(bray_mat, centers = k, nstart = 25)
  ch     <- calinhara(bray_mat, km$cluster)
  
  # Pairwise PERMANOVA
  grps   <- unique(labels)
  pairs  <- combn(grps, 2, simplify = FALSE)
  pw     <- lapply(pairs, function(g) {
    idx  <- labels %in% g
    sub_dist <- as.dist(bray_mat[idx, idx])
    fit  <- adonis2(sub_dist ~ labels[idx], permutations = 999)
    data.frame(g1 = g[1], g2 = g[2], R2 = fit$R2[1],
               F = fit[["F"]][1], p = fit$`Pr(>F)`[1])
  })
  pw_df  <- do.call(rbind, pw)
  pw_df$p_adj <- p.adjust(pw_df$p, method = "fdr")
  
  list(k = k, R2 = R2, perm_p = pval,
       disp_p = disp_p, CH = ch,
       pairwise = pw_df, labels = labels)
}

results <- lapply(3:9, eval_k)

# Summary table
summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(K = r$k, PERMANOVA_R2 = round(r$R2, 3),
             PERMANOVA_p = r$perm_p,
             BetaDisp_p  = round(r$disp_p, 3),
             Calinski_Harabasz = round(r$CH, 1))
}))

cat("\n── Cluster evaluation summary ──────────────────────────────\n")
print(summary_df)

cat("\n── Pairwise PERMANOVA per K ────────────────────────────────\n")
for (r in results) {
  cat(sprintf("\nK = %d\n", r$k))
  print(r$pairwise)
}
# cols_df and pch_df from your original code
cols_df <- c("#FC9B2D", "#7ACBC3")
names(cols_df) <- c("C", "U")
pch_df <- c(4, 15, 17, 18, 19)
names(pch_df) <- c("0", "7", "14", "21", "35")

# ── Ordination wrapper ─────────────────────────────────────────────────────────
# PCA operates on the transformed abundance matrix directly.
# PCoA (aka MDS) decomposes the Bray-Curtis distance matrix — more appropriate
#   for beta diversity since it preserves the ecological distance structure.
# NMDS is non-metric and best for stress-based visualization of Bray-Curtis.
# All three are worth comparing; PCoA is the most theoretically consistent with
#   the Bray-Curtis clustering you've already done.

run_ordination <- function(log_df, bray_dist, method = c("PCA", "PCoA", "NMDS")) {
  method <- match.arg(method)
  if (method == "PCA") {
    ord    <- prcomp(log_df, center = FALSE, scale. = FALSE)
    scores <- as.data.frame(ord$x[, 1:2])
    ev     <- summary(ord)$importance[2, ] * 100
    colnames(scores) <- c("Axis1", "Axis2")
    labs   <- c(sprintf("PC1 (%.1f%%)", ev[1]), sprintf("PC2 (%.1f%%)", ev[2]))
  } else if (method == "PCoA") {
    ord    <- wcmdscale(bray_dist, k = 2, eig = TRUE)
    scores <- as.data.frame(ord$points)
    ev     <- ord$eig / sum(ord$eig[ord$eig > 0]) * 100
    colnames(scores) <- c("Axis1", "Axis2")
    labs   <- c(sprintf("PCoA1 (%.1f%%)", ev[1]), sprintf("PCoA2 (%.1f%%)", ev[2]))
  } else {
    set.seed(42)
    ord    <- metaMDS(log_df, distance = "bray", k = 2, trymax = 200,
                      autotransform = TRUE, trace = 0)
    scores <- as.data.frame(scores(ord, display = "sites"))
    colnames(scores) <- c("Axis1", "Axis2")
    labs   <- c(sprintf("NMDS1 (stress = %.3f)", ord$stress), "NMDS2")
  }
  list(scores = scores, axis_labels = labs)
}

# Perform a PERMANOVA to test if catechin days 14, 21, and 35 as a sample cluster
# are significantly different from all other samples. Also compute beta dispersion.

# ── 1. Create the Comparison Grouping ──────────────────────────────────────────

# We need to extract the Treatment and Day from the row names of your matrix/distance object
# Your rownames look like "C_14", "U_0.1", etc.
sample_info <- data.frame(row_names = rownames(log_df)) %>%
  mutate(
    treatment = sub("_.*", "", row_names),
    day = as.numeric(sub(".*_", "", row_names)),
    # Define the target group: Catechin (C) AND Day 14, 21, or 35
    comparison_group = ifelse(treatment == "C" & day >= 14, "Catechin_Late", "Other")
  )

# Verify the grouping counts
print(table(sample_info$comparison_group))

# ── 2. PERMANOVA (adonis2) ────────────────────────────────────────────────────
# Testing: Is the 'centroid' of Catechin Day 14-35 different from others?

perm_target <- adonis2(bray_dist ~ comparison_group, 
                       data = sample_info, 
                       permutations = 999)

cat("\n── PERMANOVA: Catechin (Day 14-35) vs All Others ──\n")
print(perm_target)

# ── 3. Beta Dispersion (betadisper) ──────────────────────────────────────────
# Testing: Is the 'spread' (variance) within the groups different?
# High significance here means your PERMANOVA p-value might be driven by 
# differences in variance rather than just location.

disp_target <- betadisper(bray_dist, sample_info$comparison_group)
disp_anova  <- anova(disp_target)

cat("\n── Beta Dispersion Test (Permutation P-value) ──\n")
print(disp_anova)

# ── 4. Visualization ──────────────────────────────────────────────────────────

# PCoA Plot specifically highlighting this comparison
ord_res <- run_ordination(log_df, bray_dist, method = "PCoA")
plot_df <- cbind(ord_res$scores, sample_info)

ggplot(plot_df, aes(x = Axis1, y = Axis2, color = comparison_group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Catechin_Late" = "#D55E00", "Other" = "#56B4E9")) +
  labs(title = "PCoA: Catechin (Days 14-35) vs. Rest of Samples",
       subtitle = paste("PERMANOVA p =", perm_target$`Pr(>F)`[1]),
       x = ord_res$axis_labels[1], 
       y = ord_res$axis_labels[2]) +
  theme_minimal()



##################
## TEST catechin day 14-35 against all other samples
##################
# ── 1. Create the Comparison Grouping (Base R) ───────────────────────────────

# Extract Treatment (everything before the first underscore)
treatments <- sub("_.*", "", rownames(log_df))

# Extract Day (grabs the first number found after the underscore)
# This handles "C_0.1", "C_35", and "C_35.1" correctly
day_matches <- regmatches(rownames(log_df), regexpr("\\d+\\.?\\d*", rownames(log_df)))
days <- as.numeric(day_matches)

sample_info <- data.frame(
  row_names = rownames(log_df),
  treatment = treatments,
  day = days
) %>%
  mutate(
    # Define the target group
    comparison_group = ifelse(treatment == "C" & day >= 14, "Catechin_Late", "Other")
  )

# ──PERMANOVA & Beta Dispersion ────────────────────────────────────────────

# PERMANOVA: Tests if the 'average' community differs
perm_target <- adonis2(bray_dist ~ comparison_group, 
                       data = sample_info, 
                       permutations = 999)

# Beta Dispersion: Tests if the 'variance/spread' differs
# A non-significant p-value here is actually GOOD—it means your 
# PERMANOVA result isn't just an artifact of unequal variance.
disp_target <- betadisper(bray_dist, sample_info$comparison_group)
disp_anova  <- anova(disp_target)

cat("\n── PERMANOVA RESULTS ──\n")
print(perm_target)

cat("\n── BETA DISPERSION RESULTS ──\n")
print(disp_anova)

# ANOSIM
# Group 1: Catechin Day 14-35 | Group 2: Everything Else
sample_info$is_target <- ifelse(sample_info$treatment == "C" & sample_info$day >= 14, 
                                "Catechin_Late", "Others")

# ── 2. ANOSIM (The "Tightness" Test) ──────────────────────────────────────
# R-value near 1 suggests the groups are completely distinct.
# R-value near 0 suggests no difference in grouping.
target_anosim <- anosim(bray_dist, sample_info$is_target, permutations = 999)

cat("\n── ANOSIM Results ──\n")
print(target_anosim)


# Compare for Unamended day 14-35 vs Catechin day 14-35
# Define 4 groups: Unamended Early, Unamended Late, Catechin Early, Catechin Late
sample_info <- sample_info %>%
  mutate(
    time_period = ifelse(day >= 14, "Late", "Early"),
    specific_group = paste(treatment, time_period, sep = "_")
  )

# Verify the new groups
print(table(sample_info$specific_group))

# Filter distance matrix to ONLY include "Late" samples for a direct comparison
late_indices <- which(sample_info$time_period == "Late")
late_dist    <- as.dist(as.matrix(bray_dist)[late_indices, late_indices])
late_meta    <- sample_info[late_indices, ]

# 2. PERMANOVA: Catechin Late vs Unamended Late
perm_late <- adonis2(late_dist ~ treatment, 
                     data = late_meta, 
                     permutations = 999)

# 3. Beta Dispersion: Catechin Late vs Unamended Late
disp_late <- betadisper(late_dist, late_meta$treatment)
disp_late_pval <- permutest(disp_late, permutations = 999)$tab$`Pr(>F)`[1]

cat("\n── Comparison: Catechin Late vs. Unamended Late ──\n")
print(perm_late)
cat(sprintf("Beta Dispersion p-value: %.4f\n", disp_late_pval))

# ── 4. Visualization ──────────────────────────────────────────────────────────

# PCoA Plot specifically highlighting this comparison
ord_res <- run_ordination(log_df, bray_dist, method = "PCoA")
plot_df <- cbind(ord_res$scores, sample_info)

ggplot(plot_df, aes(x = Axis1, y = Axis2, color = comparison_group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Catechin_Late" = "#D55E00", "Other" = "#56B4E9")) +
  labs(title = "PCoA: Catechin (Days 14-35) vs. Rest of Samples",
       subtitle = paste("PERMANOVA p =", perm_target$`Pr(>F)`[1]),
       x = ord_res$axis_labels[1], 
       y = ord_res$axis_labels[2]) +
  theme_minimal()




##################
## EXTRA CODE   ##
##################

days <- sapply(strsplit(rownames(log_df), "day"), function(x) x[2])
days <- sapply(strsplit(days, '\\.'), function(x) x[1])
treatments <- sapply(strsplit(rownames(log_df), "_"), function(x) x[1])

df_types2 <- as.data.frame(cbind(days, treatments))
colnames(df_types2) <- c('day', 'treatment')
rownames(df_types2) <- rownames(log_df)

# ── Main plotting function ─────────────────────────────────────────────────────
plot_clusters <- function(log_df, bray_dist, df_types2, k_values = c(3, 4, 5, 6),
                          methods = c("PCA", "PCoA", "NMDS"),
                          figDir = "./") {
  
  for (method in methods) {
    ord  <- run_ordination(log_df, bray_dist, method)
    pts  <- ord$scores
    labs <- ord$axis_labels
    pts$treatment <- df_types2$treatment
    pts$day       <- as.character(df_types2$day)
    rownames(pts) <- rownames(df_types2)
    
    plot_list <- lapply(k_values, function(k) {
      # Get AP cluster labels for this K
      ap_k <- apclusterK(negDistMat(r = 2), as.matrix(bray_dist),
                         K = k, details = TRUE, seed = 42)
      cluster_labels <- rep(NA_character_, nrow(pts))
      for (i in seq_along(ap_k@clusters))
        cluster_labels[ap_k@clusters[[i]]] <- paste0("Cluster ", i)
      pts$cluster <- cluster_labels
      
      # Pull stats from summary_df for this K
      stats <- summary_df[summary_df$K == k, ]
      stats_label <- sprintf(
        "PERMANOVA R² = %.3f, p = %.3f\nBeta dispersion p = %.3f\nCalinski-Harabasz = %.1f",
        stats$PERMANOVA_R2, stats$PERMANOVA_p,
        stats$BetaDisp_p,   stats$Calinski_Harabasz
      )
      
      # Convex hulls per cluster
      hulls <- pts %>%
        group_by(cluster) %>%
        slice(chull(Axis1, Axis2))
      
      ggplot(pts, aes(x = Axis1, y = Axis2)) +
        geom_polygon(data = hulls,
                     aes(group = cluster, fill = cluster),
                     alpha = 0.12, color = "grey40",
                     linetype = "dashed", linewidth = 0.6) +
        geom_point(aes(color = treatment, shape = day), size = 5, stroke = 1.2) +
        geom_text(aes(label = day, color = treatment),
                  size = 3.5, vjust = -1, show.legend = FALSE) +
        scale_color_manual(values = cols_df) +
        scale_shape_manual(values = pch_df) +
        scale_fill_brewer(palette = "Set2") +
        annotate("text", x = Inf, y = Inf, label = stats_label,
                 hjust = 1.05, vjust = 1.3, size = 3.5,
                 color = "grey30", fontface = "italic") +
        labs(title   = sprintf("%s - K = %d clusters", method, k),
             x       = labs[1], y = labs[2],
             color   = "Treatment", shape = "Day", fill = "Cluster") +
        theme_minimal(base_size = 13) +
        theme(panel.grid      = element_blank(),
              axis.line       = element_line(color = "black", linewidth = 0.4),
              plot.title      = element_text(hjust = 0.5),
              legend.position = "right")
    })
    
    # One PDF per ordination method, one page per K
    pdf(file     = paste0(figDir, "cluster_comparison_", method, ".pdf"),
        width    = 10, height = 9)
    lapply(plot_list, print)
    dev.off()
    message("Saved: ", figDir, "cluster_comparison_", method, ".pdf")
  }
}

# ── Run ────────────────────────────────────────────────────────────────────────
plot_clusters(log_df, bray_dist, df_types2,
              k_values = c(3, 4, 5, 6),
              methods  = c("PCA", "PCoA", "NMDS"),
              figDir   = figDir)

# Check Beta dispersion without day 0 cluster:
# For K=4, get cluster labels
ap_4 <- apclusterK(negDistMat(r = 2), as.matrix(bray_dist),
                   K = 5, details = TRUE, seed = 42)
cluster_labels_4 <- rep(NA_character_, nrow(as.matrix(bray_dist)))
for (i in seq_along(ap_4@clusters))
  cluster_labels_4[ap_4@clusters[[i]]] <- paste0("Cluster_", i)
names(cluster_labels_4) <- rownames(as.matrix(bray_dist))

# Inspect which cluster is day0
table(cluster_labels_4, df_types2$day)

# Identify and remove the day0 cluster
day0_cluster <- names(which(table(cluster_labels_4, df_types2$day)[, "0"] > 0 &
                              rowSums(table(cluster_labels_4, df_types2$day)[, colnames(table(cluster_labels_4, df_types2$day)) != "0"]) == 0))

keep <- cluster_labels_4 != day0_cluster
sub_labels <- cluster_labels_4[keep]
sub_dist   <- as.dist(as.matrix(bray_dist)[keep, keep])

# Retest beta dispersion without day0 cluster
disp_sub  <- betadisper(sub_dist, sub_labels)
anova(disp_sub)

# Retest PERMANOVA without day0 cluster
adonis2(sub_dist ~ sub_labels, permutations = 999)

# ── PCA Plot matched to 12x12 PDF Scaling ─────────────────────────────────────
plot_pca_k5_final <- function(log_df, bray_dist, df_types2, k_values = 5, 
                              figDir = "./") {
  
  method <- "PCA"
  ord  <- run_ordination(log_df, bray_dist, method)
  pts  <- ord$scores
  labs <- ord$axis_labels
  
  pts$treatment <- df_types2$treatment
  pts$day       <- as.character(df_types2$day)
  rownames(pts) <- rownames(df_types2)
  
  # Preserving your requested PCH mapping
  pch_mapping <- c("0" = 19, "7" = 7, "14" = 17, "21" = 15, "35" = 3)
  cols_df <- c("C" = "#FC9B2D", "U" = "#7ACBC3")
  
  plot_list <- lapply(k_values, function(k) {
    # Cluster labels
    ap_k <- apclusterK(negDistMat(r = 2), as.matrix(bray_dist),
                       K = k, details = TRUE, seed = 42)
    cluster_labels <- rep(NA_character_, nrow(pts))
    for (i in seq_along(ap_k@clusters))
      cluster_labels[ap_k@clusters[[i]]] <- paste0("Cluster ", i)
    pts$cluster <- cluster_labels
    
    # PERMANOVA stats
    stats <- summary_df[summary_df$K == k, ]
    stats_label <- sprintf("PERMANOVA R² = %.3f\np = %.3f",
                           stats$PERMANOVA_R2, stats$PERMANOVA_p)
    
    # Hulls
    hulls <- pts %>% group_by(cluster) %>% slice(chull(Axis1, Axis2))
    
    ggplot(pts, aes(x = Axis1, y = Axis2)) +
      # Hulls matched to your polygon code: black solid line, width 1
      geom_polygon(data = hulls, aes(group = cluster),
                   fill = NA, color = "black", linetype = "solid", linewidth = 1) +
      
      # Points scaled to size 15 to match your previous text labels
      geom_point(aes(color = treatment, shape = day), size = 15, stroke = 2) +
      
      scale_color_manual(values = cols_df) +
      scale_shape_manual(values = pch_mapping) +
      
      # Corner Annotation scaled for 12x12
      annotate("text", x = Inf, y = Inf, label = stats_label,
               hjust = 1.1, vjust = 1.5, size = 10, fontface = "bold") +
      
      labs(x = labs[1], y = labs[2], title = "") +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        legend.position = "none",
        panel.grid = element_blank(), # Clean background like your example
        axis.line = element_line(color = "black", linewidth = 0.5)
      )
  })
  
  # Exporting as 12x12 PDF
  output_path <- paste0(figDir, "02-B_PCA_treatment_day_polygon_shapes.pdf")
  pdf(file = output_path, width = 12, height = 12)
  lapply(plot_list, print)
  dev.off()
  
  message("Saved 12x12 PCA to: ", output_path)
}

# ── Run ──────────────────────────────────────────────────────────────────────
plot_pca_k5_final(log_df, bray_dist, df_types2, k_values = 5, figDir = figDir)

# Define values for K to analyze
k_to_check <- c(4, 5)

# List to store plots for the final output
disp_plots <- list()

for (k in k_to_check) {
  
  # 1. Generate Cluster Labels for K
  ap_res <- apclusterK(negDistMat(r = 2), as.matrix(bray_dist), 
                       K = k, details = TRUE, seed = 42)
  
  cluster_labels <- rep(NA_character_, nrow(as.matrix(bray_dist)))
  for (i in seq_along(ap_res@clusters)) {
    cluster_labels[ap_res@clusters[[i]]] <- paste0("Cluster_", i)
  }
  names(cluster_labels) <- rownames(as.matrix(bray_dist))
  
  # 2. Identify the Day 0 Cluster
  # Find which cluster index has the most '0' days
  day_tab <- table(cluster_labels, df_types2$day)
  day0_cluster <- rownames(day_tab)[which.max(day_tab[, "0"])]
  
  # --- STATS: WITH DAY 0 ---
  disp_full <- betadisper(bray_dist, cluster_labels)
  a_full <- anova(disp_full)
  p_full <- a_full$`Pr(>F)`[1]
  
  # --- STATS: WITHOUT DAY 0 ---
  keep <- cluster_labels != day0_cluster
  sub_dist <- as.dist(as.matrix(bray_dist)[keep, keep])
  sub_labs <- cluster_labels[keep]
  
  disp_sub <- betadisper(sub_dist, sub_labs)
  a_sub <- anova(disp_sub)
  p_sub <- a_sub$`Pr(>F)`[1]
  
  # 3. Print Results to Console
  cat("\n--- Results for K =", k, " ---\n")
  cat("Day 0 Cluster Identified as:", day0_cluster, "\n")
  cat("ANOVA (Full) p-value:   ", round(p_full, 5), "\n")
  cat("ANOVA (No Day 0) p-value:", round(p_sub, 5), "\n")
  
  # 4. Create Comparison Plot
  df_full <- data.frame(Dist = disp_full$distances, Group = disp_full$group, Type = "With Day 0")
  df_sub  <- data.frame(Dist = disp_sub$distances, Group = disp_sub$group, Type = "Without Day 0")
  plot_df <- rbind(df_full, df_sub)
  
  disp_plots[[as.character(k)]] <- ggplot(plot_df, aes(x = Group, y = Dist, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
    facet_wrap(~Type, scales = "free_x") +
    labs(title = paste("Beta Dispersion: K =", k),
         subtitle = sprintf("Full p = %.4f | Sub p = %.4f", p_full, p_sub),
         y = "Distance to Centroid", x = "Cluster") +
    theme_minimal() +
    theme(legend.position = "none")
}


# Print each plot individually
print(disp_plots[["4"]])
print(disp_plots[["5"]])

# ── Targeted Pairwise Comparisons ─────────────────────────────────────────────
# Comparison 1: Unamended Day 7  vs  Catechin Day 14
# Comparison 2: Unamended Day 14 vs  Catechin Day 14

run_targeted_comparison <- function(bray_dist, label1, label2, 
                                    group1_name, group2_name,
                                    permutations = 999) {
  
  bray_mat <- as.matrix(bray_dist)
  all_names <- rownames(bray_mat)
  
  # Grab indices for each group (handles replicates via grepl)
  idx1 <- grep(paste0("^", label1), all_names)
  idx2 <- grep(paste0("^", label2), all_names)
  
  if (length(idx1) == 0) stop(paste("No samples found matching:", label1))
  if (length(idx2) == 0) stop(paste("No samples found matching:", label2))
  
  cat(sprintf("\n══════════════════════════════════════════════════════\n"))
  cat(sprintf("Comparison: %s  vs  %s\n", group1_name, group2_name))
  cat(sprintf("  Samples in group 1 (%s): %s\n", group1_name,
              paste(all_names[idx1], collapse = ", ")))
  cat(sprintf("  Samples in group 2 (%s): %s\n", group2_name,
              paste(all_names[idx2], collapse = ", ")))
  cat(sprintf("══════════════════════════════════════════════════════\n"))
  
  # Subset distance matrix
  idx_all  <- c(idx1, idx2)
  sub_dist <- as.dist(bray_mat[idx_all, idx_all])
  sub_labs <- c(rep(group1_name, length(idx1)),
                rep(group2_name, length(idx2)))
  names(sub_labs) <- all_names[idx_all]
  
  # ── PERMANOVA ──────────────────────────────────────────────────────────────
  set.seed(42)
  perm_res <- adonis2(sub_dist ~ sub_labs, permutations = permutations)
  
  cat("\n── PERMANOVA ─────────────────────────────────────────────\n")
  print(perm_res)
  cat(sprintf("  R²  = %.4f\n", perm_res$R2[1]))
  cat(sprintf("  F   = %.4f\n", perm_res[["F"]][1]))
  cat(sprintf("  p   = %.4f\n", perm_res$`Pr(>F)`[1]))
  
  # ── Beta Dispersion ────────────────────────────────────────────────────────
  disp_res  <- betadisper(sub_dist, sub_labs)
  disp_anov <- anova(disp_res)
  
  cat("\n── Beta Dispersion (ANOVA) ───────────────────────────────\n")
  print(disp_anov)
  cat(sprintf("  Dispersion p = %.4f\n", disp_anov$`Pr(>F)`[1]))
  
  # Per-group distances to centroid
  cat("\n── Distance to Centroid (mean ± SD) ──────────────────────\n")
  dist_df <- data.frame(
    sample = names(disp_res$distances),
    dist   = disp_res$distances,
    group  = disp_res$group
  )
  cent_summary <- dist_df %>%
    group_by(group) %>%
    summarise(
      n    = n(),
      mean = round(mean(dist), 4),
      sd   = round(sd(dist),   4),
      .groups = "drop"
    )
  print(cent_summary)
  
  # ── Permutation test for dispersion ───────────────────────────────────────
  set.seed(42)
  disp_perm <- permutest(disp_res, permutations = permutations, pairwise = TRUE)
  cat("\n── Beta Dispersion (permutation test) ────────────────────\n")
  print(disp_perm)
  
  # Return everything for downstream use
  invisible(list(
    group1       = group1_name,
    group2       = group2_name,
    samples_g1   = all_names[idx1],
    samples_g2   = all_names[idx2],
    sub_dist     = sub_dist,
    sub_labs     = sub_labs,
    permanova    = perm_res,
    betadisper   = disp_res,
    disp_anova   = disp_anov,
    disp_perm    = disp_perm,
    cent_summary = cent_summary
  ))
}

# ── Run the two comparisons ────────────────────────────────────────────────────

# Comparison 1: Unamended Day 7 vs Catechin Day 14
comp1 <- run_targeted_comparison(
  bray_dist   = bray_dist,
  label1      = "U_day7",
  label2      = "C_day7",
  group1_name = "U_day7",
  group2_name = "C_day7"
)

# Comparison 2: Unamended Day 14 vs Catechin Day 14
comp2 <- run_targeted_comparison(
  bray_dist   = bray_dist,
  label1      = "U_day14",
  label2      = "C_day14",
  group1_name = "U_day14",
  group2_name = "C_day14"
)

# ── Summary table across both comparisons ─────────────────────────────────────
comparison_summary <- data.frame(
  Comparison       = c("U_day7 vs C_day7", "U_day14 vs C_day14"),
  PERMANOVA_R2     = c(comp1$permanova$R2[1],          comp2$permanova$R2[1]),
  PERMANOVA_F      = c(comp1$permanova[["F"]][1],      comp2$permanova[["F"]][1]),
  PERMANOVA_p      = c(comp1$permanova$`Pr(>F)`[1],    comp2$permanova$`Pr(>F)`[1]),
  BetaDisp_p_anova = c(comp1$disp_anova$`Pr(>F)`[1],  comp2$disp_anova$`Pr(>F)`[1])
) %>%
  mutate(
    PERMANOVA_R2     = round(PERMANOVA_R2,     4),
    PERMANOVA_F      = round(PERMANOVA_F,      4),
    PERMANOVA_p      = round(PERMANOVA_p,      4),
    BetaDisp_p_anova = round(BetaDisp_p_anova, 4),
    # FDR correction across the two PERMANOVA tests
    PERMANOVA_p_fdr  = round(p.adjust(PERMANOVA_p, method = "fdr"), 4),
    # Interpretation flag
    Sig_PERMANOVA    = ifelse(PERMANOVA_p      < 0.05, "YES", "no"),
    Homogeneous_Disp = ifelse(BetaDisp_p_anova > 0.05, "YES (ok)", "NO (unequal)")
  )

cat("\n══════════════════════════════════════════════════════════════\n")
cat("SUMMARY: Targeted Pairwise Comparisons\n")
cat("══════════════════════════════════════════════════════════════\n")
print(comparison_summary)


library(vegan)
library(dplyr)

run_anosim_comparison <- function(bray_dist, label1, label2, group1_name, group2_name) {
  # 1. Subset the distance matrix
  bray_mat  <- as.matrix(bray_dist)
  idx1      <- grep(paste0("^", label1), rownames(bray_mat))
  idx2      <- grep(paste0("^", label2), rownames(bray_mat))
  
  idx_all   <- c(idx1, idx2)
  sub_dist  <- as.dist(bray_mat[idx_all, idx_all])
  sub_labs  <- c(rep(group1_name, length(idx1)), rep(group2_name, length(idx2)))
  
  # 2. Run ANOSIM
  set.seed(42)
  res <- anosim(sub_dist, sub_labs, permutations = 999)
  
  # 3. Output summary
  cat(sprintf("\n--- %s vs %s ---\n", group1_name, group2_name))
  cat(sprintf("R statistic: %.4f | p-value: %.4f\n", res$statistic, res$signif))
  
  return(res)
}

# --- Execute Comparisons ---
comp1 <- run_anosim_comparison(bray_dist, "U_day7",  "C_day7",  "U_day7",  "C_day7")
comp2 <- run_anosim_comparison(bray_dist, "U_day14", "C_day14", "U_day14", "C_day14")

# --- Summary Table ---
anosim_summary <- data.frame(
  Comparison = c("U_day7 vs C_day7", "U_day14 vs C_day14"),
  R_stat     = c(comp1$statistic, comp2$statistic),
  p_value    = c(comp1$signif,    comp2$signif)
) %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"),
         Significant = ifelse(p_fdr < 0.05, "Yes", "No"))

print(anosim_summary)
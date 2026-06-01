library(vegan)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)

# ── SETTINGS ──────────────────────────────────────────────────────────────────
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")
figDir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/"
metadata_path <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv"

# Helper function to extract and filter distances
get_uc_distances <- function(dist_obj, source_name) {
  mat <- as.matrix(dist_obj)
  df_dist <- as.data.frame(as.table(mat))
  colnames(df_dist) <- c("S1", "S2", "Distance")
  
  # Helper to parse names like "C_35.1"
  parse_name <- function(x) {
    parts <- strsplit(as.character(x), "_")[[1]]
    treat <- parts[1]
    day_raw <- parts[2]
    day_num <- gsub("\\..*", "", day_raw) # Remove make.unique suffixes
    return(data.frame(Treat = treat, Day = day_num))
  }
  
  m1 <- do.call(rbind, lapply(df_dist$S1, parse_name))
  m2 <- do.call(rbind, lapply(df_dist$S2, parse_name))
  
  df_dist <- cbind(df_dist, m1, m2)
  colnames(df_dist)[4:7] <- c("T1", "D1", "T2", "D2")
  
  # Filter: Same day, U vs C only, skip day0
  final <- df_dist %>%
    filter(D1 == D2, D1 != "day0", T1 == "U", T2 == "C") %>%
    mutate(Source = source_name, Day = factor(D1, levels = c("day7", "day14", "day21", "day35")))
  
  return(final)
}

# ── 1. PROCESS vOTU DATA ──────────────────────────────────────────────────────
votu_raw <- read.csv('../results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
rownames(votu_raw) <- votu_raw$vOTU
votu_raw$vOTU <- NULL

# Load Metadata to get treatment order
meta_tab <- read.csv(metadata_path) %>%
  mutate(treatment = case_when(treatment == "unamended" ~ "U", treatment == "catechin" ~ "C", TRUE ~ "T")) %>%
  filter(treatment != "T")

# Sync names (assuming column order in df matches metadata rows after filtering T)
colnames(votu_raw) <- make.unique(paste(meta_tab$treatment, meta_tab$timepoint, sep="_"))
votu_hel <- decostand(t(votu_raw), method = "hellinger")
votu_dist <- vegdist(votu_hel, method = "bray")
votu_final <- get_uc_distances(votu_dist, "vOTU")

# ── 2. PROCESS MAG DATA ───────────────────────────────────────────────────────
mag_raw <- read_excel('../results/MAGs/metaT/MAG_abundance_table.xlsx')

# Cleanup CT samples
ct_samples <- c('STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127', 'GTDB')

mag_clean <- mag_raw %>% select(-any_of(ct_samples))
mag_df <- as.data.frame(mag_clean)
rownames(mag_df) <- mag_df[[1]]
mag_df <- mag_df[,-1]
mag_df <- mag_df[rowSums(mag_df) > 0, ]

colnames(mag_df) <- make.unique(paste(meta_tab$treatment, meta_tab$timepoint, sep="_"))
mag_hel <- decostand(t(mag_df), method = "hellinger")
mag_dist <- vegdist(mag_hel, method = "bray")
mag_final <- get_uc_distances(mag_dist, "MAG")

# ── 3. COMBINE AND PLOT ───────────────────────────────────────────────────────
combined_df <- rbind(votu_final, mag_final)

# Side-by-Side Plot
ggplot(combined_df, aes(x = Day, y = Distance, fill = Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Source), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 1.5, alpha = 0.5) +
  facet_wrap(~Source) +
  scale_fill_manual(values = c("vOTU" = "#D55E00", "MAG" = "#0072B2")) +
  scale_color_manual(values = c("vOTU" = "#882255", "MAG" = "#117733")) +
  labs(title = "Community Divergence: Unamended vs. Catechin",
       subtitle = "Bray-Curtis Distances (Hellinger Transformed)",
       y = "Dissimilarity (Bray-Curtis)",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

library(ggpubr)

# 1. Define the comparisons you want to show with brackets
# Since you want to show Day 7 is different from the rest:
my_comparisons <- list( 
  c("day7", "day14"), 
  c("day7", "day21"), 
  c("day7", "day35") 
)

# 2. Define the significance symbols
symnum_args <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
  symbols = c("***", "**", "*", "ns")
)

# 3. Clean the dataframe names to ensure no duplicates exist
combined_df_clean <- combined_df[, !duplicated(colnames(combined_df))]

# 4. Generate the Plot
ggplot(combined_df_clean, aes(x = Day, y = Distance, fill = Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Source), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1.5, alpha = 0.5) +
  facet_wrap(~Source) +
  
  # Add brackets and p-values
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    symnum.args = symnum_args,
    step.increase = 0.1, # Shifts brackets upward so they don't overlap
    tip.length = 0.02
  ) +
  
  # Styling
  scale_fill_manual(values = c("vOTU" = "#D55E00", "MAG" = "#0072B2")) +
  scale_color_manual(values = c("vOTU" = "#882255", "MAG" = "#117733")) +
  labs(title = "Community Divergence: Unamended vs. Catechin",
       subtitle = "Pairwise T-tests against Day 7 baseline",
       y = "Dissimilarity (Bray-Curtis)",
       x = "Timepoint") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

library(multcompView)
library(dplyr)
library(ggplot2)

# ── 1. Enhanced CLD Function using Wilcoxon Rank Sum ────────────────────────
get_ordered_cld <- function(data, source_label) {
  df_sub <- data %>% filter(Source == source_label)
  
  # ── Order factor levels by mean Distance (low → high) ────────────────────
  # This ensures 'a' corresponds to the groups with the lowest dissimilarity
  day_order <- df_sub %>%
    group_by(Day) %>%
    summarise(mean_dist = mean(Distance), .groups = "drop") %>%
    arrange(mean_dist) %>%
    pull(Day)
  
  df_sub$Day <- factor(df_sub$Day, levels = day_order)
  
  # ── Pairwise Wilcoxon Test (Non-Parametric) ──────────────────────────────
  # Using FDR (Benjamini-Hochberg) to correct for multiple comparisons
  wilcox_res <- pairwise.wilcox.test(df_sub$Distance, 
                                     df_sub$Day, 
                                     p.adjust.method = "fdr")
  
  # ── Convert Pairwise Matrix to Vector for multcompLetters ────────────────
  # We extract the p-values and format them as a named vector (e.g., "Day0-Day7")
  p_matrix <- wilcox_res$p.value
  
  # Create a named vector from the matrix
  p_vec <- as.vector(p_matrix)
  names(p_vec) <- outer(rownames(p_matrix), colnames(p_matrix), paste, sep="-")
  p_vec <- p_vec[!is.na(p_vec)] # Remove empty comparisons from the triangle
  
  # Generate the Compact Letter Display
  cld_obj <- multcompLetters(p_vec)
  
  # Extract letters in the sorted factor order (low mean → high mean)
  cld_letters <- cld_obj$Letters[day_order]
  
  return(data.frame(
    Day     = names(cld_letters),
    Letters = as.character(cld_letters),
    Source  = source_label
  ))
}

# ── 2. Re-calculate CLDs (Non-Parametric) ───────────────────────────────────
all_letters_fixed <- bind_rows(
  get_ordered_cld(combined_df_clean, "MAG"),
  get_ordered_cld(combined_df_clean, "vOTU")
)

# ── 3. Label positions (top of whiskers) ─────────────────────────────────────
y_label_pos <- combined_df_clean %>%
  group_by(Source, Day) %>%
  summarise(y_max = max(Distance) + 0.05, .groups = "drop") %>%
  left_join(all_letters_fixed, by = c("Source", "Day"))

# ── 4. Final Publication-Ready Plot ──────────────────────────────────────────
# We set the family to "Helvetica" at the device level
pdf(file = paste(figDir, "SX_CLD_bray_distance_nonparametric.pdf", sep=''), 
    width = 10, height = 10, 
    family = "Helvetica")

ggplot(combined_df_clean, aes(x = Day, y = Distance, fill = Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Source),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 2.5, alpha = 0.4) + # Increased point size for visibility
  facet_wrap(~Source) +
  geom_text(data = y_label_pos, aes(x = Day, y = y_max, label = Letters),
            size = 8,             # Increased from 5 (Letters above boxes)
            fontface = "bold", 
            family = "Helvetica") +
  scale_fill_manual(values  = c("vOTU" = "#D55E00", "MAG" = "#0072B2")) +
  scale_color_manual(values = c("vOTU" = "#882255", "MAG" = "#117733")) +
  labs(
    title    = "Statistical Divergence of Viral and Microbial Communities",
    subtitle = "Shared letters indicate no significant difference (Pairwise Wilcoxon Rank Sum, FDR-adjusted p > 0.05)",
    y        = "Bray-Curtis Dissimilarity (U vs. C)",
    x        = "Day"
  ) +
  # base_size = 18 bumps up all axes and labels globally
  theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(
    strip.text       = element_text(face = "bold", size = 20), # Facet titles
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 22),
    plot.subtitle    = element_text(hjust = 0.5, size = 14),
    axis.title       = element_text(size = 20),                # X and Y axis titles
    axis.text        = element_text(size = 16, color = "black"), # Tick labels
    panel.spacing    = unit(2, "lines")                        # More space between facets
  )

dev.off()
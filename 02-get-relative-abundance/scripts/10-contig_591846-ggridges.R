library(ggplot2)
library(ggridges)
library(dplyr)
library(gggenomes)
library(cowplot)

# --- 1. Load and Process Data ---
# Ridgeline Data
pileup_data <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/pileup_data/combined_pileup.csv")
pileup_data$sample <- gsub("_90FILTERED", "", pileup_data$sample)

metadata <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv")
metadata$sample <- metadata$Sample

pileup_data <- pileup_data %>%
  left_join(metadata, by = "sample") %>%
  mutate(sample = treatment_day_replicate) %>%
  filter(grepl("^catechin", sample, ignore.case = TRUE))
pileup_data$day <- factor(pileup_data$day)

# BLAST Alignment Data
blast_file <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/results/contig_591846_STM_0716_E_M_E034_A_bin.10_blastn.txt"
# Skip comment lines and define headers based on the "Fields" line in your file
blast_data <- read.table(blast_file, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
colnames(blast_data) <- c("qacc", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Genome Data
wdir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/JAGFXR01_vOTUs"
seqs <- read_seqs(paste(wdir, "contig_591846.fna", sep='/'))
genes <- read_feats(paste(wdir, "contig_591846.gff", sep='/'))
genome_length <- max(seqs$length, na.rm = TRUE)

# --- 2. Create Plots ---

# A. Ridgeline Plot (Top)
# Calculate relative heights for gridlines based on ggridges internal scaling
max_cov <- max(pileup_data$coverage, na.rm = TRUE)
grid_values <- c(2000, 4000, 6000, 8000)

# Create a dataframe for the gridlines
gridlines <- expand.grid(
  day = unique(pileup_data$day),
  cov = grid_values
) %>%
  mutate(
    # Map the day factor to its integer position on the y-axis
    y_base = as.numeric(factor(day, levels = rev(levels(pileup_data$day)))), 
    # Calculate the relative height based on your scale = 0.9
    rel_height = (cov / max_cov) * 0.9,
    y_pos = y_base + rel_height
  )

# A. Ridgeline Plot (Top) - With Manual Gridlines
ridgeline_plot <- ggplot(pileup_data, aes(x = position, y = day, height = coverage, group = sample)) +
  # Add the custom gridlines behind the ridges
  geom_segment(data = gridlines, 
               aes(x = 0, xend = genome_length, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, color = "gray70", linetype = "dashed") +
  geom_density_ridges(
    stat = "identity", 
    scale = 0.9, 
    alpha = 0.2, 
    fill = "grey", 
    color = "black"
  ) +
  scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black"), 
    axis.ticks.y = element_line(color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 0, 5)
  ) +
  labs(y = "Time Point")

# B. BLAST Alignment Plot (Middle)
# This shows bars where the phage (query) matches the host (subject)
blast_plot <- ggplot(blast_data) +
  geom_segment(aes(x = qstart, xend = qend, y = 0, yend = 0), 
               color = "red", size = 3) +
  # Add black outlines/ends for better visibility if they overlap
  geom_segment(aes(x = qstart, xend = qstart + 10, y = 0, yend = 0), color = "black", size = 3) +
  scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 0, 5),
        axis.title.y = element_text(size = 10, angle = 90, vjust = 0.5)) +
  labs(y = "Host Align")

# C. Genome Plot (Bottom)
genome_plot <- gggenomes(genes = genes, seqs = seqs) +
  geom_seq() +
  geom_gene(aes(fill = strand)) +
  scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10),
        axis.line.x = element_line(), axis.ticks.x = element_line(),
        plot.margin = margin(0, 5, 5, 5)) +
  labs(x = "Genomic Position (bp)")

# --- 3. Align and Combine ---
combined_plot <- plot_grid(
  ridgeline_plot,
  blast_plot,
  genome_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(3, 0.4, 1) # Adjust height of blast bars here
)

# --- 4. Save ---
ggsave("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/03-C_contig_591846_ridgeline_blast_genome.png", 
       plot = combined_plot, width = 12, height = 7, dpi = 300)

# --- Compute Unique Coverage ---

unique_coverage <- blast_data %>%
  # 1. Ensure start is always less than end (BLAST can occasionally flip these)
  mutate(
    start = pmin(qstart, qend),
    end = pmax(qstart, qend)
  ) %>%
  # 2. Sort by start position
  arrange(start) %>%
  # 3. Merge overlapping intervals
  # We identify a new "island" whenever a start position is greater than 
  # the maximum end position seen so far.
  mutate(new_island = start > lag(cummax(end), default = first(end) - 1)) %>%
  mutate(island_id = cumsum(new_island)) %>%
  # 4. Collapse the islands to find the min start and max end of each
  group_by(island_id) %>%
  summarise(
    island_start = min(start),
    island_end = max(end),
    .groups = "drop"
  ) %>%
  # 5. Calculate width of each unique island and sum them up
  mutate(width = (island_end - island_start) + 1) %>%
  summarise(total_unique_bases = sum(width))

# Print results
total_bases <- unique_coverage$total_unique_bases
percent_coverage <- (total_bases / genome_length) * 100

cat("Total Unique Bases Covered:", total_bases, "bp\n")
cat("Percent of Genome Covered:", round(percent_coverage, 2), "%\n")

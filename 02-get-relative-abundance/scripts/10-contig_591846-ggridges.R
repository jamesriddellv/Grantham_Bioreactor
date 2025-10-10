library(ggplot2)
library(ggridges)
library(dplyr)

# Read the combined data
pileup_data <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/pileup_data/combined_pileup.csv")
pileup_data$sample <- gsub("_90FILTERED", "", pileup_data$sample)

metadata <- read.csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv")
metadata$sample <- metadata$Sample

pileup_data <- pileup_data %>%
  left_join(metadata, by = "sample") %>%
  mutate(sample = treatment_day_replicate) %>%
  filter(grepl("^catechin", sample, ignore.case = TRUE)) %>%
  select(-treatment_day_replicate)

pileup_data$day <- factor(pileup_data$day)

ridgeline_plot <- ggplot(pileup_data, aes(x = position, y = coverage, height = coverage, group = day)) +
  geom_density_ridges(
    stat = "identity",
    scale = 1,
    alpha = 0.4,
    fill = "lightgrey",   # lighter fill for shading
    color = "black",
    rel_min_height = 0    # ensures the fill reaches the baseline
  ) +
  facet_wrap(~ day, ncol = 1, scales = "fixed") +  # <-- fixed y-axis range across days
  labs(
    x = "Position",
    y = "Coverage"
  ) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # tighten vertical gaps
    strip.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

print(ridgeline_plot)

# Save plot
ggsave("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/03-C_contig_591846_ridgeline.pdf", plot = ridgeline_plot,width = 20, height = 10, dpi = 300)
ggsave("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/03-C_contig_591846_ridgeline.png", plot = ridgeline_plot,width = 20, height = 10, dpi = 300)

# Load required libraries
library(gggenomes)
library(dplyr)
library(ggplot2)
library(ggridges)
library(cowplot)

# Your existing genome plot code
wdir <- "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/JAGFXR01_vOTUs"

seqs <- read_seqs(c(
  paste(wdir, "contig_591846.fna", sep='/')
))

genes <- read_feats(c(
  paste(wdir, "contig_591846.gff", sep='/')
))

# Get genome length for consistent x-axis
genome_length <- max(seqs$length, na.rm = TRUE)

# Create ridgeline plot (on top)
ridgeline_plot <- pileup_data %>%
  ggplot(aes(x = position, y = day, height = coverage, group = sample)) +
  geom_density_ridges(
    stat = "identity",
    scale = 0.9,
    alpha = 0.2,
    fill = "grey",
    color = "black"
  ) +
  scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  theme_ridges() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 1, 5)  # Reduce bottom margin
  ) +
  labs(y = "Time Point")

# Create genome plot (on bottom)
genome_plot <- gggenomes(genes = genes, seqs = seqs) +
  geom_seq() +
  geom_gene(aes(fill = strand)) +
  scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    plot.margin = margin(1, 5, 5, 5)  # Reduce top margin
  ) +
  labs(x = "Genomic Position (bp)")

# Use cowplot to align the plots
combined_plot <- plot_grid(
  ridgeline_plot,
  genome_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(2, 1)
)

print(combined_plot)

ggsave("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/03-C_contig_591846_ridgeline_genome.pdf", plot = combined_plot,width = 10, height = 5, dpi = 150)
ggsave("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/03-C_contig_591846_ridgeline_genome.png", plot = combined_plot,width = 10, height = 5, dpi = 150)

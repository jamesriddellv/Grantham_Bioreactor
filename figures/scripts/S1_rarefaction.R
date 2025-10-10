library(vegan)
library(viridis)
library(dplyr)

setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures")

# read in manually curated vOTU coverm mean stats (>=90% ANI, >=75% hcov)
data <- read.csv('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/results/coverm_mean/combined_vOTU_mean_90_75.tsv', sep='\t')

# drop zero rows
data <- data[rowSums(data[-1]) > 0,]

sample_metadata <- read.csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaG_sample_metadata.csv')

# rename and sort columns

# clean up column names to match sample_metadata
for (col in 1:ncol(data)) {
  colnames(data)[col] <- sub("_90FILTERED_SORTED\\.Mean|_MG.*", "", colnames(data)[col])
}

# sort columns alphabetically to match sample_metadata
data <- data[, order(colnames(data))]

# rename with sample_metadata
sample_metadata$col.label <- paste0(sample_metadata$treatment, "_", gsub("day", "", sample_metadata$timepoint))
colnames(data)[1] <- "vOTU"
colnames(data)[-1] <- sample_metadata$col.label[match(colnames(data)[-1], sample_metadata$Sample)]

data.mtx <- t(as.matrix(data))

colnames(data.mtx) <- data.mtx[1, ]
data.mtx <- data.mtx[-1, ]

data.mtx <- apply(data.mtx, 2, as.numeric)



ceiling.mtx <- ceiling(data.mtx)

S <- specnumber(ceiling.mtx)

(raremax <- min(rowSums(ceiling.mtx)))

Srare <- rarefy(ceiling.mtx, raremax)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

pdf("S1_rarefaction_curve.pdf", width = 8, height = 6)

samples <- colnames(data)[-1]
colors <- rainbow(length(samples))  # Make sure colors match the number of samples

# Increase right margin
par(mar = c(5, 4, 4, 12))  

# Create the rarefaction curves
rarecurve(ceiling.mtx, step = 20, sample = raremax, col = colors, cex = 0.6)

# Add the legend inside the margin instead of fully outside
legend("topright", inset = c(-0.45, 0), legend = samples, col = colors, pch = 1, cex = 1, xpd=TRUE)

dev.off()

## MetaT ##

data <- read.csv('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/results/metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED_no0s.tsv', header=F, sep='\t')

colnames(data)=c(
  "gene",
  "STM_0716_E_M_E002",
  "STM_0716_E_M_E003",
  "STM_0716_E_M_E004",
  "STM_0716_E_M_E025",
  "STM_0716_E_M_E027",
  "STM_0716_E_M_E029",
  "STM_0716_E_M_E030",
  "STM_0716_E_M_E031",
  "STM_0716_E_M_E033",
  "STM_0716_E_M_E034",
  "STM_0716_E_M_E035",
  "STM_0716_E_M_E050",
  "STM_0716_E_M_E051",
  "STM_0716_E_M_E052",
  "STM_0716_E_M_E054",
  "STM_0716_E_M_E055",
  "STM_0716_E_M_E056",
  "STM_0716_E_M_E058",
  "STM_0716_E_M_E059",
  "STM_0716_E_M_E060",
  "STM_0716_E_M_E062",
  "STM_0716_E_M_E063",
  "STM_0716_E_M_E064",
  "STM_0716_E_M_E066",
  "STM_0716_E_M_E067",
  "STM_0716_E_M_E068",
  "STM_0716_E_M_E070",
  "STM_0716_E_M_E071",
  "STM_0716_E_M_E072",
  "STM_0716_E_M_E121",
  "STM_0716_E_M_E122",
  "STM_0716_E_M_E123",
  "STM_0716_E_M_E125",
  "STM_0716_E_M_E126",
  "STM_0716_E_M_E127",
  "STM_0716_E_M_E129",
  "STM_0716_E_M_E130",
  "STM_0716_E_M_E131"
)

# drop last two rows of metadata
data <- head(data, -2)

vOTU_gene_ids <- read.delim('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt')
vOTU_gene_ids <- vOTU_gene_ids %>% select(-length)

data <- merge(x=data, y=vOTU_gene_ids, by='gene', all.x=TRUE)


data <- data %>%
  group_by(vOTU) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

rownames(data) <- data$vOTU


# filter out vOTUs with <50% genes mapped in at least one sample
votus_to_keep <- read.delim("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/active_vOTUs_50_percent_genes_mapped.tsv", header = TRUE, sep = "\t")

# Filter your main dataframe based on those vOTUs
data <- data[data$vOTU %in% votus_to_keep$vOTU, ]

data$vOTU <- NULL

metadata <- read.csv("/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")
metadata <- metadata %>%
  mutate(treatment = case_when(
    treatment == "unamended" ~ "U",
    treatment == "CT" ~ "T",
    treatment == "catechin" ~ "C",
    TRUE ~ treatment
  ))

treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')

colnames(data) <- treatment_timepoint

# make unique column names for df to investigate separately
colnames(data) <- make.unique(names(data))

data <- data %>%
  select(-starts_with("T"))

data.mtx <- t(as.matrix(data))


data.mtx <- apply(data.mtx, 2, as.numeric)

S <- specnumber(data.mtx)

(raremax <- min(rowSums(data.mtx)))

Srare <- rarefy(data.mtx, raremax)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

pdf("SX_rarefaction_curve_metaT.pdf", width = 8, height = 6)

samples <- colnames(data)[-1]
colors <- rainbow(length(samples))  # Make sure colors match the number of samples

# Increase right margin
par(mar = c(5, 4, 4, 12))  

# Create the rarefaction curves
rarecurve(data.mtx, step = 20, sample = raremax, col = colors, cex = 0.6)

# Add the legend inside the margin instead of fully outside
legend("topright", inset = c(-0.45, 0), legend = samples, col = colors, pch = 1, cex = 1, xpd=TRUE)

dev.off()

# Plot rarefied species richness over time and treatment
rarefied_richness <- metadata
rarefied_richness$rarefied_richness <- Srare

# Load required package
library(ggplot2)

# Drop CT treatment
rarefied_richness <- subset(rarefied_richness, treatment != "T")

cols_df <- c("C" = "#FC9B2D", "U" = "#7ACBC3")

rarefied_richness <- rarefied_richness %>% mutate(timepoint = factor(timepoint, levels = c("day0", "day7", "day14", "day21", "day35")))

ggplot(rarefied_richness, aes(x = timepoint, y = rarefied_richness, group = treatment, color = treatment)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +  # Line for mean richness
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # Error bars (optional)
  stat_summary(fun = mean, geom = "point", size = 3) +   # Points for mean
  labs(
    title = "Rarefied Richness over Time by Treatment",
    x = "Timepoint",
    y = "Rarefied Richness",
    color = "Treatment"
  ) +
  theme_minimal()


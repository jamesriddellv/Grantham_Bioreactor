# excerpted from https://www.reneshbedre.com/blog/expression_units.html

# Load required libraries
library(dplyr)
library(EdgeR)

# set working
setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts")

# Read data into data frames
readcounts_df <- read.csv('../results/ALL-SAMPLES-coverm_readcounts.txt', sep=',')
names(readcounts_df) <- gsub("\\.Mean", "", colnames(readcounts_df))
readcounts_df <- readcounts_df[, order(names(readcounts_df))]

total_reads <- read.csv('../results/coverm_stats.txt', sep=',')[c('Sample', 'Total.Mapped')]
total_reads$Sample <- gsub("'", "", total_reads$Sample)

length <- read.table('../results/vOTU_length.txt', sep='\t', header=TRUE)
names(length) <- c('Contig', 'Length')

# Sort data frames
readcounts_df <- readcounts_df %>%
  mutate(Contig = factor(Contig, levels = unique(Contig))) %>%
  arrange(Contig)

length <- length %>%
  mutate(Contig = factor(Contig, levels = unique(Contig))) %>%
  arrange(Contig)

# Merge data frames
readcounts_length <- merge(readcounts_df, length, by='Contig', all.x=TRUE)

# Assuming you want to set the 'Contig' column as row names
row.names(readcounts_length) <- readcounts_length$Contig

# If 'Contig' is not needed as a separate column in the data frame
readcounts_length <- readcounts_length[, -which(names(readcounts_length) == 'Contig')]

# gene length is in bp in expression dataset and converted to Kbp
rpk <- readcounts_length[,1:38] / (readcounts_length$Length / 1000)

# load metadata
metadata <- read.csv('../data/metadata.txt', sep='\t')

# comparing groups
group <- factor(metadata$treatment)
y <- DGEList(counts=rpk, group=group)

# normalize for library size by cacluating scaling factor using TMM (default method)
y <- calcNormFactors(y)
y$samples


# count per million read (normalized count)
norm_counts <- cpm(y)
head(norm_counts)

# write to csv

write.csv(norm_counts, file = "../results/GeTMM_normalized_counts_readcounts.csv", row.names = TRUE)
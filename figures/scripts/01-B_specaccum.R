### Species Accumulation Curve ###

# Load necessary libraries
library(vegan)  # For ecological data analysis
library(ggplot2)  # For plotting
library(dplyr)

setwd("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures")

# read in manually curated vOTU metaG coverm mean stats (>=90% ANI, >=75% hcov)
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

rownames(data) <- data$vOTU
data$vOTU <- NULL

data.mtx <- t(as.matrix(data))

data.mtx <- apply(data.mtx, 2, as.numeric)

# Create a species accumulation curve
spec_accum <- specaccum(data.mtx > 0, method = "random")

# Plot the species accumulation curve
ggplot() +
  geom_ribbon(aes(x = seq_along(spec_accum$richness), ymin = spec_accum$richness - spec_accum$sd, ymax = spec_accum$richness + spec_accum$sd), fill = "lightblue") +
  geom_line(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "red", linewidth = 1) +
  geom_point(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "black", size = 2) +
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +  # Adjust x-axis to show whole numbers from 1 to 10
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),   # Adjust axis title size
    axis.text = element_text(size = 14),    # Adjust axis text size
  )

# Extract the richness values and sample indices
richness <- spec_accum$richness
n_samples <- length(richness)

# Define the last 10 samples
last_4_indices <- (n_samples - 3):n_samples
last_4_richness <- richness[last_4_indices]

# Perform linear regression to get the slope
slope_model <- lm(last_4_richness ~ last_4_indices)

# Extract the slope (coefficient) and intercept
slope <- coef(slope_model)[2]
intercept <- coef(slope_model)[1]

# Compute the fitted line across the full range of samples
full_line <- data.frame(
  x = 1:n_samples,
  y = intercept + slope * (1:n_samples)  # y = mx + b
)

# Plot the species accumulation curve with the extended slope
ggplot() +
  # Ribbon for standard deviation
  geom_ribbon(aes(x = seq_along(richness), 
                  ymin = spec_accum$richness - spec_accum$sd, 
                  ymax = spec_accum$richness + spec_accum$sd), 
              fill = "lightblue") +
  # Main species accumulation curve
  geom_line(aes(x = seq_along(richness), y = richness), 
            color = "red", size = 1) +
  geom_point(aes(x = seq_along(richness), y = richness), 
             color = "black", size = 2) +
  # Extended slope line across the full x-axis
  geom_line(data = full_line, aes(x = x, y = y), 
            color = "grey", linetype = "dashed", size = 0.5) +
  # Slope label in the upper left
  annotate("text", 
           x = min(full_line$x) + 1,  # Slightly offset from the left
           y = max(full_line$y) + 2,  # Slightly offset from the top
           label = paste0("Slope = ", round(slope, 2)), 
           color = "black", 
           size = 3, 
           hjust = 0) +  # Left-align text
  # Labels and adjustments
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  # X-axis labeled 1 through 10
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    axis.title = element_text(size = 10),   
    axis.text = element_text(size = 8)    
  )

ggsave("S2_vOTU_specaccum_metaG.pdf", width = 8, height = 6, dpi = 300)

### Species accumulation metaG no presence cutoff ### 

# read in manually curated vOTU metaG coverm mean stats (>=90% ANI, >=75% hcov)
data <- read.csv('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/results/coverm_mean/no_minfrac/metaG_vOTU_mean_90_75_0.tsv', sep='\t')

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

rownames(data) <- data$vOTU
data$vOTU <- NULL

data.mtx <- t(as.matrix(data))

data.mtx <- apply(data.mtx, 2, as.numeric)

# Create a species accumulation curve
spec_accum <- specaccum(data.mtx > 0, method = "random")

# Plot the species accumulation curve
ggplot() +
  geom_ribbon(aes(x = seq_along(spec_accum$richness), ymin = spec_accum$richness - spec_accum$sd, ymax = spec_accum$richness + spec_accum$sd), fill = "lightblue") +
  geom_line(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "red", linewidth = 1) +
  geom_point(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "black", size = 2) +
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +  # Adjust x-axis to show whole numbers from 1 to 10
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),   # Adjust axis title size
    axis.text = element_text(size = 14),    # Adjust axis text size
  )

# Extract the richness values and sample indices
richness <- spec_accum$richness
n_samples <- length(richness)

# Define the last 10 samples
last_4_indices <- (n_samples - 3):n_samples
last_4_richness <- richness[last_4_indices]

# Perform linear regression to get the slope
slope_model <- lm(last_4_richness ~ last_4_indices)

# Extract the slope (coefficient) and intercept
slope <- coef(slope_model)[2]
intercept <- coef(slope_model)[1]

# Compute the fitted line across the full range of samples
full_line <- data.frame(
  x = 1:n_samples,
  y = intercept + slope * (1:n_samples)  # y = mx + b
)

ggplot() +
  # Ribbon for standard deviation
  geom_ribbon(aes(x = seq_along(richness), 
                  ymin = spec_accum$richness - spec_accum$sd, 
                  ymax = spec_accum$richness + spec_accum$sd), 
              fill = "lightblue") +
  # Main species accumulation curve
  geom_line(aes(x = seq_along(richness), y = richness), 
            color = "red", size = 1) +
  geom_point(aes(x = seq_along(richness), y = richness), 
             color = "black", size = 2) +
  # Extended slope line across the full x-axis
  geom_line(data = full_line, aes(x = x, y = y), 
            color = "grey", linetype = "dashed", size = 0.5) +
  # Slope label in the upper left
  annotate("text", 
           x = min(full_line$x) + 1,  # Slightly offset from the left
           y = max(full_line$y) + 2,  # Slightly offset from the top
           label = paste0("Slope = ", round(slope, 2)), 
           color = "black", 
           size = 3, 
           hjust = 0) +  # Left-align text
  # Labels and adjustments
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  # X-axis labeled 1 through 10
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    axis.title = element_text(size = 10),   
    axis.text = element_text(size = 8)    
  )

ggsave("S2_vOTU_specaccum_metaG_no_presence_cutoff.pdf", width = 8, height = 6, dpi = 300)

### Species accumulation metaT ###

library(dplyr)
library(tidyr)

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

data <- data %>%
  mutate(vOTU = sub("_.*", "", gene))

data$gene <- NULL

df <- data %>%
  group_by(vOTU) %>%
  summarise(across(starts_with("STM"), sum)) %>%
  data.frame()



rownames(df) <- df$vOTU
df$vOTU <- NULL

data.mtx <- t(as.matrix(df))

data.mtx <- apply(data.mtx, 2, as.numeric)

# Create a species accumulation curve
spec_accum <- specaccum(data.mtx > 0, method = "random")

# Plot the species accumulation curve
ggplot() +
  geom_ribbon(aes(x = seq_along(spec_accum$richness), ymin = spec_accum$richness - spec_accum$sd, ymax = spec_accum$richness + spec_accum$sd), fill = "lightblue") +
  geom_line(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "red", size = 1) +
  geom_point(aes(x = seq_along(spec_accum$richness), y = spec_accum$richness), color = "black", size = 2) +
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  scale_x_continuous(breaks = seq(1, 38, by = 5), limits = c(1, 38)) +  # Adjust x-axis to show whole numbers from 1 to 38
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),   # Adjust axis title size
    axis.text = element_text(size = 14),    # Adjust axis text size
  )

diff(spec_accum$richness)

# Load required libraries
library(ggplot2)
library(vegan)  # For specaccum

# Generate the species accumulation curve
spec_accum <- specaccum(data.mtx > 0, method = "random")

# Extract the richness values and sample indices
richness <- spec_accum$richness
n_samples <- length(richness)

# Define the last 10 samples
last_10_indices <- (n_samples - 9):n_samples
last_10_richness <- richness[last_10_indices]

# Perform linear regression to get the slope
slope_model <- lm(last_10_richness ~ last_10_indices)

# Extract the slope (coefficient) and intercept
slope <- coef(slope_model)[2]
intercept <- coef(slope_model)[1]

# Compute the fitted line across the full range of samples
full_line <- data.frame(
  x = 1:n_samples,
  y = intercept + slope * (1:n_samples)  # y = mx + b
)

# Plot the species accumulation curve with the extended slope
ggplot() +
  # Ribbon for standard deviation
  geom_ribbon(aes(x = seq_along(richness), 
                  ymin = spec_accum$richness - spec_accum$sd, 
                  ymax = spec_accum$richness + spec_accum$sd), 
              fill = "lightblue") +
  # Main species accumulation curve
  geom_line(aes(x = seq_along(richness), y = richness), 
            color = "red", size = 1) +
  geom_point(aes(x = seq_along(richness), y = richness), 
             color = "black", size = 2) +
  # Extended slope line across the full x-axis
  geom_line(data = full_line, aes(x = x, y = y), 
            color = "grey", linetype = "dashed", size = 1) +
  # Slope label in the upper left
  annotate("text", 
           x = min(full_line$x) + 1,  # Slightly offset from the left
           y = max(full_line$y) + 2,  # Slightly offset from the top
           label = paste0("Slope = ", round(slope, 2)), 
           color = "black", 
           size = 3, 
           hjust = 0) +  # Left-align text
  # Labels and adjustments
  labs(x = "Samples", y = "vOTUs (cumulative)") +
  # X-axis labeled 1 through 10
  scale_x_continuous(breaks = seq(1, 38, by = 5), limits = c(1, 38)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White background
    axis.title = element_text(size = 10),   
    axis.text = element_text(size = 8)    
  )


ggsave("01-B_vOTU_specaccum_metaT.pdf", width = 8, height = 6, dpi = 300)


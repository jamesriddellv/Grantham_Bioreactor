library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# load relative abundance table
df <- read.csv('data/GeTMM_normalized_counts_readcounts.csv')
names(df)[names(df) == 'X'] <- 'Contig'

# load length table
length_df <- read.csv('data/vOTU_length.txt', sep='\t')

# filter to 10kb only
length_df_10kb <- length_df %>% filter(STM_0716_E_M_E002.Length > 10000)
contigs_10kb <- length_df_10kb$Contig
df_10kb <- df[df$Contig %in% contigs_10kb,]
abundance_table <- df_10kb[,-1]
non_zero_rows <- rowSums(abundance_table != 0) > 0
filtered_df <- df_10kb[non_zero_rows,]

# compute total abundance
filtered_df$total_abundance <- rowSums(abundance_table[non_zero_rows,])

total_df <- filtered_df[,c('Contig', 'total_abundance')]

# compute rank
total_df$rank <- rank(-total_df$total_abundance, ties.method='random')

top_100 <- total_df[total_df$rank >= 1 & total_df$rank <= 20,]

# plot
rankplot <- ggplot(top_100, aes(x=rank, y=total_abundance)) +  # set x and y axis
  geom_point(shape=21, size=2.25, aes(color=Contig, fill=Contig)) + # plot data points and colour them according to species.
  geom_line(colour = "grey") +  # plot line to connect data points
  theme_bw() +  # set theme for plot
  labs(x="Species Rank",y="Total relative abundance (%)", # label the axis
       title = "Top 100 rank abundance viruses") +  # give your plot a title
  theme(panel.grid.major = element_blank(),  # get rid of grid lines
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),  # set axis font size
        legend.position = "top",  # position the legend at the top
        legend.background = element_rect(colour="black"))  # outline the legend
ggsave("total_rank_plot.png", plot=rankplot, width=8, height=6, units="in", dpi=300)



## MEAN ##

filtered_df$mean_abundance <- rowMeans(abundance_table[non_zero_rows,])

mean_df <- filtered_df[,c('Contig', 'mean_abundance')]

# compute rank
mean_df$rank <- rank(-mean_df$mean_abundance, ties.method='random')
top_20 <- mean_df[mean_df$rank >= 1 & mean_df$rank <= 20,]




## Histogram ##

histogram <- ggplot() +
  geom_histogram(aes(x = mean_df$mean_abundance), fill = "skyblue", color = "black", bins = 100) +
  labs(x = "Mean rel. abundance", y = "Frequency", title = "Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.00     0.13     2.21    61.45    16.44 39488.84")

ggsave("figures/mean_abundance_histogram.png", plot=histogram, width=8, height=6, units="in", dpi=300)
# Print the histogram
print(histogram)




# plot
mean_rankplot <- ggplot(top_20, aes(x=rank, y=mean_abundance)) +  # set x and y axis
  geom_point(shape=21, size=2.25, aes(color=Contig, fill=Contig)) + # plot data points and colour them according to species.
  geom_line(colour = "grey") +  # plot line to connect data points
  theme_bw() +  # set theme for plot
  labs(x="Species Rank",y="Mean relative abundance", # label the axis
       title = "Top 20 rank mean rel. abundance across samples viruses") +  # give your plot a title
  theme(panel.grid.major = element_blank(),  # get rid of grid lines
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),  # set axis font size
        legend.position = "top",  # position the legend at the top
        legend.background = element_rect(colour="black"))  # outline the legend
ggsave("figures/top_20_mean_rank_plot.png", plot=mean_rankplot, width=8, height=6, units="in", dpi=300)

# get 99th% 
percentile_99 <- quantile(mean_df$mean_abundance, probs = 0.99)
percentile_99

common_taxa <- mean_df[mean_df$rank <= 25,]$Contig
write.csv(common_taxa, "data/common_vOTUs_via_rank_abundance.csv", row.names=FALSE)

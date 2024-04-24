################ df Dean Way ##################

#install and load all of our packages
#install.packages("tidyverse", dependencies = TRUE )
#install.packages('vegan')
#install.packages('pheatmap')
#install.packages('apcluster')
#install.packages('ggplot2')
#install.packages('corrplot')
library(BiocManager)
#BiocManager::install("microbiome")
library(vegan)
library(tidyverse)
library(apcluster)
library(corrplot)

#set location of your dataset and resulting data - change this for your own location
setwd("/Users/riddellj/Documents/Graduate_school/Sullivan_Lab/PhageTherapyOfEcosystems/Grantham_Bioreactor")

### CONTIG LEVEL ANALYSIS ###

# we have summed the abundance of all genes to produce a contig-level analysis

#read in and format abundance data
df <- read.csv('data/5kb_vOTU_geTMM_relative_abundance.tsv', sep='\t')
rownames(df) <- df$Contig
df <- df[!(names(df) == "Contig") & !(names(df) == "X")]

# remove contigs with all zero abundance values
df <- subset(df, !apply(df, 1, function(row) all(row == 0)))

metadata <- read.csv('data/sample_metadata.csv')
metadata <- metadata %>%
  mutate(treatment = case_when(
    treatment == "unamended" ~ "U",
    treatment == "CT" ~ "T",
    treatment == "catechin" ~ "C",
    TRUE ~ treatment
  ))

treatment_timepoint <- paste(metadata$treatment, metadata$timepoint, sep='_')

colnames(df) <- treatment_timepoint

# group by treatment and timepoint
df_avg <- as.data.frame(t(rowsum(t(df), names(df)) / c(table(names(df)))))

# make unique column names for df to investigate separately
colnames(df) <- make.unique(names(df))

#log transformation
#note this is ln
t_df <- t(df) %>% as_tibble(rownames = NA)
log_df <- log(t_df +1)


#bray curtis
bray_df <- vegdist(log_df, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)


#some data will need to be a matrix or dataframe later
bray_matrix_df <- as.matrix(bray_df)
write.table(bray_matrix_df, file = "data/df_bray_curtis.csv")

#my prefered way of generating figures - faster, works in base R
pdf(file = "figures/eco_stats/df_bray.pdf", width = 10, height = 10);
heatmap(bray_matrix_df, main = "Bray Curtis")
dev.off()


#affinity propagation
df_ap_clust <- apcluster(negDistMat(r=3), bray_matrix_df, details = TRUE)

#r here is a power used to scale the data and allow for more easily recognizable separations
pdf(file = "figures/eco_stats/df_ap_cluster.pdf", width = 10, height = 10);
heatmap(df_ap_clust)
dev.off()
df_ap_clust

#writing to a tibble a list of samples and their respective cluster
#set the number of clusters generated in apcluster() here
df_clusters <- 4

#create blank tibble
df_types <- tibble()

types <- tibble()
#add each cluster to tibble
for (i in 1:df_clusters) {
  df_temp_types <- data.frame("site" = names(df_ap_clust[[i]]),
                               "type" = paste0("Cluster ", i))
  
  df_types <- bind_rows(df_types, df_temp_types)
}



#load the type data
df_type_no_change <- df_types
sites_df <- df_types$site
df_types <- as.data.frame(df_types)
df_types <- as.data.frame(df_types[, -1])
rownames(df_types) <- sites_df
names(df_types)[1] <- "type"
fac_df <-factor(df_types$type)
cols_df <- c("red", "blue", "orange", "green", "violet", "gray", "cyan")

#ordination
NMDS_df <- metaMDS(log_df, distance = "bray", k=3, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_df

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "figures/eco_stats/NMDS_stress_df.pdf", width = 12, height = 9);
stress_df <- stressplot(NMDS_df)
plot(stress_df)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "figures/eco_stats/NMDS_df.pdf", width = 10, height = 10);
plot(NMDS_df, type="p", display="sites")
points(NMDS_df, display = "sites", col = cols_df[fac_df], pch = 19)
legend("topright", legend=levels(fac_df), col=cols_df, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.1, y = 0.06, "Stress = 0.1129278")
text(x = -0.1, y = 0.05, "Linear Fit = 0.941")
dev.off()


### NMDS by treatment ###

#TODO: WORK HERE MORE figure this out to plot NMDS by treatment
# Separate and remove "day" from the "day" column and cast it to an integer

treatments = sapply(strsplit(rownames(log_df), "_"), function(x) x[1])
df_types = as.data.frame(treatments)
colnames(df_types) <- 'type'
rownames(df_types) <- rownames(log_df)

#load the type data
fac_df <-factor(df_types$type)
cols_df <- c("red", "blue", "orange", "green", "violet", "gray", "cyan")
#ordination
NMDS_df <- metaMDS(log_df, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_df

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "figures/eco_stats/NMDS_treatment_stress_df.pdf", width = 12, height = 9);
stress_df <- stressplot(NMDS_df)
plot(stress_df)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "figures/eco_stats/NMDS_treatment_df.pdf", width = 10, height = 10);
plot(NMDS_df, type="p", display="sites")
points(NMDS_df, display = "sites", col = cols_df[fac_df], pch = 19)
legend("topright", legend=levels(fac_df), col=cols_df, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.5, y = 0.50, "Stress = 0.1129278")
text(x = -0.5, y = 0.45, "Linear Fit = 0.941")
dev.off()


### by day ###

days <- sapply(strsplit(rownames(log_df), "day"), function(x) x[2])
treatments = sapply(strsplit(days, '\\.'), function(x) x[1])

df_types = as.data.frame(treatments)
colnames(df_types) <- 'type'
rownames(df_types) <- rownames(log_df)

#load the type data
fac_df <-factor(df_types$type)
cols_df <- c("red", "blue", "orange", "green", "violet", "gray", "cyan")
#ordination
NMDS_df <- metaMDS(log_df, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_df

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "figures/eco_stats/NMDS_day_stress_df.pdf", width = 12, height = 9);
stress_df <- stressplot(NMDS_df)
plot(stress_df)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "figures/eco_stats/NMDS_day_df.pdf", width = 10, height = 10);
plot(NMDS_df, type="p", display="sites")
points(NMDS_df, display = "sites", col = cols_df[fac_df], pch = 19)
legend("topright", legend=levels(fac_df), col=cols_df, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.5, y = 0.50, "Stress = 0.1129278")
text(x = -0.5, y = 0.45, "Linear Fit = 0.941")
dev.off()


#envfit, vector fitting for relating environmental data to ordination
# dfface_env_table <-read.table("df_data_time_scale.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
# dfface_env_table
# attach(dfface_env_table)
# NMDS_meta_df = envfit(NMDS_df~Time+Salinity+Temp, permutations=9999)
# NMDS_meta_df
# pdf(file = "figures/eco_stats/envfit_df.pdf", width = 10, height = 10);
# plot(NMDS_df, type="p", display="sites")
# points(NMDS_df, display = "sites", col = cols_df[fac_df], pch = 19)
# legend("topright", legend=levels(fac_df), col=cols_df, pch = 19)
# text(x = -0.1, y = 0.06, "Stress = 0.07051535")
# text(x = -0.1, y = 0.05, "Linear Fit = 0.988")
# # please fill in the correct values above
# plot(NMDS_meta_df)
# dev.off()



#arrows represent strength of r2 values

#for the ambitious - give beta-dispersion a try here https://rdrr.io/cran/vegan/man/betadisper.html

# inverse simpson, simple alpha diversity
inv_simp_df <- diversity(log_df, index = "invsimpson", groups = df_types$type)
inv_simp_whole_df <- diversity(log_df, index = "invsimpson")

inv_simp_df
inv_simp_whole_df


#tibble now with site, diversity and cluster
inv_simp_whole_tib_df <- tibble("site" = names(inv_simp_whole_df),
                                 "inv_simp" = inv_simp_whole_df,
                                 "cluster" = df_types)
inv_simp_whole_tib_df <- data.frame(site = names(inv_simp_whole_df),
                                     inv_simp = inv_simp_whole_df,
                                     cluster = df_types); inv_simp_whole_tib_df
colnames(inv_simp_whole_tib_df)[colnames(inv_simp_whole_tib_df) == 'type'] <- 'cluster'

write.table(inv_simp_whole_tib_df, file = "inv_simpson_t_df.csv")



alphaD_tib_df <- inner_join(x = df_type_no_change, y = inv_simp_whole_tib_df, by = "site")
pdf(file = "figures/eco_stats/gginv_simp_by_cluster_df.pdf", width = 10, height = 10);
alphaD_plot_df <- ggplot(alphaD_tib_df, aes(x=cluster, y=inv_simp)) +
  geom_boxplot()
alphaD_plot_df
dev.off()

# get alpha diversity by treatment over number of days

# Use separate() to split the "site" column
alphaD_tib_df <- separate(alphaD_tib_df, site, into = c("treatment", "day"), sep = "_|\\.")

# Remove "day" from the "day" column and cast it to an integer
alphaD_tib_df <- alphaD_tib_df %>%
  mutate(day = as.integer(sub("day", "", day)))

# Calculate the average inverse Simpson index for each treatment group and day
average_data <- aggregate(inv_simp ~ treatment + day, data = alphaD_tib_df, FUN = mean)

# Plot the line plot with average line
pdf(file = "figures/eco_stats/gginv_simp_by_treatment_df.pdf", width = 10, height = 10)
alphaD_plot_df <- ggplot(alphaD_tib_df, aes(x = day, y = inv_simp, group = treatment, color = treatment)) +
  geom_point() +
  stat_summary(geom = "line", fun.y = mean)
  labs(x = "Day", y = "Average Inverse Simpson Index", title = "Average Inverse Simpson Index over Time by Treatment") +
  theme_minimal()
print(alphaD_plot_df)
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
kw_test_df <- kruskal.test(inv_simp_whole_df ~ cluster, data = inv_simp_whole_tib_df);kw_test_df
kw_test_df


#mantel test, are two matrices of different types related?
dfface_env_table
df_env_data_mantel <- dfface_env_table
rownames(df_env_data_mantel) <- df_env_data_mantel$Sample_id
df_env_data_mantel <- df_env_data_mantel[, -c(1,2,3,4,5)]

env_man_df <- vegdist(df_env_data_mantel, method = "manhattan", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE)
env_man_df
bray_df
mantel_df <- mantel(bray_df, env_man_df, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_df

#spearman correlations
env_man_m_df<-as.matrix(env_man_df)
env_data_m_df<-as.matrix(df_env_data_mantel)
cor_bray_env_man_df <- cor(bray_matrix_df, env_man_m_df, method = "spearman", use = "pairwise.complete.obs")
cor_bray_env_man_df
cor_bray_env_df <- cor(bray_matrix_df, env_data_m_df, method = "spearman", use = "pairwise.complete.obs")
cor_bray_env_df
cor_abund_env_df<-cor(log_df, env_data_m_df, method = "spearman", use = "pairwise.complete.obs")
cor_abund_env_df

#plotting the correlations
pdf(file = "figures/eco_stats/cor_df.pdf", width = 10, height = 10);
cor_plt<-corrplot(cor_bray_env_df, cl.ratio = 0.5)
dev.off()




########################### Shannon's H and making boxplots with kruskal wallis


# shannon's h, simple alpha diversity
shannons_df <- diversity(log_df, index = "shannon", groups = df_types$type)
shannons_whole_df <- diversity(log_df, index = "shannon")


shannons_df
shannons_whole_df


#tibble now with site, diversity and cluster
shannons_whole_tib_df <- tibble("site" = names(shannons_whole_df),
                                 "shannon_h" = shannons_whole_df,
                                 "cluster" = df_types)
shannons_whole_tib_df <- data.frame(site = names(shannons_whole_df),
                                     shannon_h = shannons_whole_df,
                                     cluster = df_types); shannons_whole_tib_df
colnames(shannons_whole_tib_df)[colnames(shannons_whole_tib_df) == 'type'] <- 'cluster'

write.table(shannons_whole_tib_df, file = "shannons_t_df.csv")



shannon_alphaD_tib_df <- inner_join(x = df_type_no_change, y = shannons_whole_tib_df, by = "site")
pdf(file = "figures/eco_stats/ggshannons_df.pdf", width = 10, height = 10);
shannon_alphaD_plot_df <- ggplot(shannon_alphaD_tib_df, aes(x=cluster, y=shannon_h)) +
  geom_boxplot()
shannon_alphaD_plot_df
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
library("ggpubr")
library(ggplot2)
df_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"),
                        c("Cluster 1", "Cluster 6"), c("Cluster 1", "Cluster 7"), c("Cluster 2", "Cluster 3"), 
                        c("Cluster 2", "Cluster 4"), c("Cluster 2", "Cluster 5"), c("Cluster 2", "Cluster 6"), c("Cluster 2", "Cluster 7"), 
                        c("Cluster 3", "Cluster 4"), c("Cluster 3", "Cluster 5"), c("Cluster 3", "Cluster 6"), 
                        c("Cluster 3", "Cluster 7"), c("Cluster 4", "Cluster 5"), c("Cluster 4", "Cluster 6"), 
                        c("Cluster 4", "Cluster 7"), c("Cluster 5", "Cluster 6"), c("Cluster 5", "Cluster 7"),
                        c("Cluster 6", "Cluster 7"))

kw_shannons_df <- kruskal.test(shannons_whole_df ~ cluster, data = shannons_whole_tib_df)
p_value_shannons_df <- kw_shannons_df$p.value

Shannon_H_df_Boxplot <- ggplot(shannons_whole_tib_df, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = df_comparisons) +
  stat_compare_means(label.y = 10.6, label.x = 2)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H dfface Data")
ggsave("figures/eco_stats/shannons_H_df_comparisons.png", plot = Shannon_H_df_Boxplot, width = 6, height = 4)



df_box_shannon <- ggplot(shannons_whole_tib_df, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  geom_text(aes(label = sprintf("Kruskal-Wallis, p = %.3f", p_value_shannons_df)), x = 2.5, y = 10.12, hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H df"); df_box_shannon

ggsave("figures/eco_stats/shannon_h_df.png", plot = df_box_shannon, width = 6, height = 4)




########################## Inverse Simpson


Inv_Simp_df_Boxplot <- ggplot(inv_simp_whole_tib_df, aes(x = cluster, y = inv_simp, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Inv_simp") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = df_comparisons) +
  stat_compare_means(label.y = 36000, label.x = 1.25)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Inverse Simpson's df Data")
ggsave("figures/eco_stats/inv_simp_df_comparisons.png", plot = Inv_Simp_df_Boxplot, width = 6, height = 4)



kw_df <- kruskal.test(inv_simp_whole_df ~ cluster, data = inv_simp_whole_tib_df)
p_value_df <- kw_df$p.value

df_box_inv_simp <- ggplot(inv_simp_whole_tib_df, aes(x = cluster, y = inv_simp, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Inv_simp") +
  xlab("Cluster")+
  geom_text(aes(label = sprintf("Kruskal-Wallis, p = %.3f", p_value_df)), x = 4.5, y = max(inv_simp_whole_tib_df$inv_simp), hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Inverse Simpson's df"); df_box_inv_simp

ggsave("figures/eco_stats/inv_simp_df.png", plot = df_box_inv_simp, width = 6, height = 4)

##############################################  MRPP test ############################################################
grouping <- df_types$type


mrpp_result <- mrpp(bray_df, grouping, permutations = 999, distance = "bray"); mrpp_result

cols_SUR_diel <- c("goldenrod3", "cyan")
SUR_types_diel <- SUR_types
SUR_types_diel$diel <- grouping_sur_day_night
fac_SUR_diel <- factor(SUR_types_diel$diel)

plot(NMDS_SUR, type="p", display="sites")
points(NMDS_df, display = "sites", col = cols_SUR_diel[fac_SUR_diel], pch = 19)
legend("topright", legend=levels(fac_SUR_diel), col=cols_SUR_diel, pch = 19)
text(x = -0.1, y = 0.08, "Stress = 0.07051535")
text(x = -0.1, y = 0.07, "Linear Fit = 0.988")
text(x = -0.1, y = 0.06, "MRPP = 0.314")

ordihull(NMDS_SUR, fac_SUR_diel, display = "sites", draw = c("polygon"), col = cols_SUR_diel, alpha = 100, label = FALSE)

Surface_env_table_timescale <-read.table("SUR_data_time_scale.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
Surface_env_table_timescale <- Surface_env_table_timescale[, -c(1,2,4,5)]
NMDS_meta_SUR_diel = envfit(NMDS_SUR~Time+Salinity+Temp, permutations=9999)
NMDS_meta_SUR_diel
#pdf(file = "envfit_SUR.pdf", width = 10, height = 10);
plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR_diel[fac_SUR_diel], pch = 19)
legend("topright", legend=levels(fac_SUR_diel), col=cols_SUR_diel, pch = 19)
text(x = -0.1, y = 0.08, "Stress = 0.07051535")
text(x = -0.1, y = 0.07, "Linear Fit = 0.988")
text(x = -0.1, y = 0.06, "MRPP = 0.314")

ordihull(NMDS_SUR, fac_SUR_diel, display = "sites", draw = c("polygon"), col = cols_SUR_diel, alpha = 100, label = FALSE)
# please fill in the correct values above
plot(NMDS_meta_SUR_diel)
dev.off()







plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR[fac_SUR], pch = 19)
legend("topright", legend=levels(fac_SUR), col=cols_SUR, pch = 19)
text(x = -0.1, y = 0.08, "Stress = 0.07051535")
text(x = -0.1, y = 0.07, "Linear Fit = 0.988")
text(x = -0.1, y = 0.06, "MRPP = 0.314")

ordihull(NMDS_SUR, fac_SUR, display = "sites", draw = c("polygon"), col = cols_SUR, alpha = 100, label = FALSE)




########## Plotting the PCoA ##################


# Dean method
#read in and format abundance data
df <- t(vOTU_df)
df <- as_tibble(df, rownames = NA)
#log transformation
#note this is ln
log_df <- log(df +1)
#bray curtis
bray_df <- vegdist(log_df, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
#ordination
NMDS_df <- metaMDS(log_df, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_df



##### testing ####
NMDS_df <- metaMDS(log_df, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
pcoa_df <- metaMDS(Bray.chron.df, distance = "bray", k=2, trymax = 000, autotransform = TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
pcoa_df_2 <- metaMDS(test_veg, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
pcoa_df_3 <- metaMDS(test_dean_veg, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)

fac_df <-factor(df_types$type)
fac_pcoa_df <- factor(pcoa.chron.bray.plotting$group)
col_pcoa_df <- c("darkgoldenrod2", "cyan3")

plot(pcoa_df, type="p", display="sites")
points(pcoa_df, display = "sites", col = cols_df_diel[fac_df_diel], pch = 19)
legend("topright", legend=levels(fac_df_diel), col=cols_df_diel, pch = 19)
ordihull(pcoa_df, fac_pcoa_df, display = "sites", draw = c("polygon"), col = col_pcoa_df, alpha = 100, label = FALSE)


plot(NMDS_df, type = "p", display = "sites")
plot(pcoa_df_2, type = "p", display = "sites")
plot(pcoa_df_3, type = "p", display = "sites")

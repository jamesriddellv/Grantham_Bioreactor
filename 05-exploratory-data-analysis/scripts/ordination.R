library(phyloseq)
library(ggplot2)
library(plyr)

GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)
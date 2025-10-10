library(edgeR)
rpk = read.delim('/users/PAS1573/riddell26/data_5x_rpk.csv', sep=',')
rownames(rpk) <- rpk$gene
rpk$gene <- NULL
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df)

write.table(getmms_df, '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs.tsv', sep="\t")

library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)

setwd("../03-mapped-reads-per-gene")

##################
## starting with reverse stranded
##################
# read in htseq table
# this table has rows with all zeros removed
counts = read.delim("htseq_vOTUs_100M_95FILTERED_REVSTRANDED_no0s.tsv",sep="\t",header = FALSE)
colnames(counts)=c("gene","STM_0716_E_M_E002","STM_0716_E_M_E003","STM_0716_E_M_E004","STM_0716_E_M_E025","STM_0716_E_M_E027","STM_0716_E_M_E029","STM_0716_E_M_E030","STM_0716_E_M_E031","STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035","STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052","STM_0716_E_M_E054","STM_0716_E_M_E055","STM_0716_E_M_E056","STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060","STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064","STM_0716_E_M_E066","STM_0716_E_M_E067","STM_0716_E_M_E068","STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072","STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123","STM_0716_E_M_E125","STM_0716_E_M_E126","STM_0716_E_M_E127","STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131")
# get lib size
lib=as.data.frame(colSums(counts[,2:39]))
lib$no_feature=t(counts[2022482,2:39])
lib$ambig=t(counts[2022483,2:39])
colnames(lib)=c("sum","no_feature","ambig")

#plot mapping stats
lib$mapped=(lib$sum-lib$no_feature-lib$ambig)
lib$sample=row.names(lib)
lib%>%filter(sample!="sum")%>%gather(-sample,-sum,key="status",value="counts")%>%
  ggplot()+
  geom_bar(aes(x=sample,y=counts,fill=status),stat="identity",position="stack")+
  theme_classic()+
  geom_vline(xintercept=50000000)+
  coord_flip()

#remove bottom two stats rows
counts=counts[1:2022481,]
row.names(counts)=counts[,1]
counts=counts[,-1]

# read in gene lengths
len=read.delim("gene_lengths.txt",sep="\t",header = FALSE)
colnames(len)=c("gene","length")

# moving on to filtering
counts$sum=rowSums(counts)
range(counts$sum) # range 1 to 975112

filtered_counts=counts %>%
  filter(sum>0)
range(filtered_counts$sum) # range 1 to 975112, because i already removed rows of all zero

#  making read counts less than 5 = 0
counts5 = filtered_counts
counts5[counts5==1|counts5==2|counts5==3|counts5==4] <-0

pa_counts=ifelse(counts5>0,1,0)
pa_counts=as.data.frame(pa_counts)
pa_counts=pa_counts[,-39] # get rid of old sum column
pa_counts$n=rowSums(pa_counts)

# just has to be in 1 sample
n2 = pa_counts %>%
  filter(n>=1)

x = row.names(n2)
x=as.data.frame(x)
x$n=n2$n
filtered_counts$name=row.names(filtered_counts)

filtered_1x = inner_join(filtered_counts,x,by=c("name"="x"))
colnames(filtered_1x)
filtered_1x = filtered_1x[,c(40,1:38)]
colnames(filtered_1x)

#write_xlsx(filtered_1x,"filtered_5counts_1sample_Rev.xlsx") # too big to write
write.csv(filtered_1x,"filtered_5counts_1sample_Rev.txt",sep = "\t")

# convert to geTMM
count=filtered_1x
count=left_join(count,len,by=c("name"="gene"))
count5=count%>%distinct()
rpk = ((count5[,2:39]*10^3)/count5$length)
row.names(rpk)=count5$name
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
lib
lib_size=c(100000000,98819300,95400980,78691532,100000000,88540314,98819300,100000000,83168470,100000000,99235256,91608098,93570944,80222604,81306340,100000000,100000000,100000000,84557472,100000000,88283948,100000000,100000000,100000000,100000000,96316002,63972406,61512094,100000000,100000000,100000000,100000000,100000000,100000000,100000000,78973250,100000000,100000000)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:39]))
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples

rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df) #0.000 8082.522

getmms_df$gene=row.names(getmms)
colnames(getmms_df)
getmms_df = getmms_df[,c(39,1:38)]
colnames(getmms_df)
#write_xlsx(getmms_df, 'getmms_ge5_1X.xlsx') # too big
write.table(getmms_df, 'getmms_ge5_1X_REV.txt', sep="\t")


##################
## now with unstranded
##################
# read in htseq table
# this table has rows with all zeros removed
counts = read.delim("htseq_2304MAGs_100M_97_UNSTRANDED_no0s.txt",sep="\t",header = FALSE)
colnames(counts)=c("gene","STM_0716_E_M_E002","STM_0716_E_M_E003","STM_0716_E_M_E004","STM_0716_E_M_E025","STM_0716_E_M_E027","STM_0716_E_M_E029","STM_0716_E_M_E030","STM_0716_E_M_E031","STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035","STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052","STM_0716_E_M_E054","STM_0716_E_M_E055","STM_0716_E_M_E056","STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060","STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064","STM_0716_E_M_E066","STM_0716_E_M_E067","STM_0716_E_M_E068","STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072","STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123","STM_0716_E_M_E125","STM_0716_E_M_E126","STM_0716_E_M_E127","STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131")
# get lib size
lib=as.data.frame(colSums(counts[,2:39]))
lib$no_feature=t(counts[2224883,2:39])
lib$ambig=t(counts[2224884,2:39])
colnames(lib)=c("sum","no_feature","ambig")

#plot mapping stats
lib$mapped=(lib$sum-lib$no_feature-lib$ambig)
lib$sample=row.names(lib)
lib%>%filter(sample!="sum")%>%gather(-sample,-sum,key="status",value="counts")%>%
  ggplot()+
  geom_bar(aes(x=sample,y=counts,fill=status),stat="identity",position="stack")+
  theme_classic()+
  geom_vline(xintercept=50000000)+
  coord_flip()
  


#remove bottom two stats rows
counts=counts[1:2224882,]
row.names(counts)=counts[,1]
counts=counts[,-1]

# read in gene lengths
len=read.delim("gene_lengths.txt",sep="\t",header = FALSE)
colnames(len)=c("gene","length")

# moving on to filtering
counts$sum=rowSums(counts)
range(counts$sum) # range 1 to 975839

filtered_counts=counts %>%
  filter(sum>0)
range(filtered_counts$sum) # range 1 to 975839, because i already removed rows of all zero

#  making read counts less than 5 = 0
counts5 = filtered_counts
counts5[counts5==1|counts5==2|counts5==3|counts5==4] <-0

pa_counts=ifelse(counts5>0,1,0)
pa_counts=as.data.frame(pa_counts)
pa_counts=pa_counts[,-39] # get rid of old sum column
pa_counts$n=rowSums(pa_counts)

# just has to be in 1 sample
n2 = pa_counts %>%
  filter(n>=1)

x = row.names(n2)
x=as.data.frame(x)
x$n=n2$n
filtered_counts$name=row.names(filtered_counts)

filtered_1x = inner_join(filtered_counts,x,by=c("name"="x"))
colnames(filtered_1x)
filtered_1x = filtered_1x[,c(40,1:38)]
colnames(filtered_1x)

#write_xlsx(filtered_1x,"filtered_5counts_1sample_Rev.xlsx") # too big to write
write.table(filtered_1x,"filtered_5counts_1sample_UNSTRANDED.txt",sep = "\t")

# convert to geTMM
count=filtered_1x
count=left_join(count,len,by=c("name"="gene"))
count5=count%>%distinct()
rpk = ((count5[,2:39]*10^3)/count5$length)
row.names(rpk)=count5$name
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
lib
lib_size=c(100000000,98819300,95400980,78691532,100000000,88540314,98819300,100000000,83168470,100000000,99235256,91608098,93570944,80222604,81306340,100000000,100000000,100000000,84557472,100000000,88283948,100000000,100000000,100000000,100000000,96316002,63972406,61512094,100000000,100000000,100000000,100000000,100000000,100000000,100000000,78973250,100000000,100000000)
lib_size=as.data.frame(lib_size)
lib_size$sample=c(colnames(filtered_1x[,2:39]))
row.names(lib_size)=lib_size$sample
lib_size=select(lib_size,-sample)

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]
rpk.norm$samples$lib.size
rpk.norm$samples

rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df) #0.000 9036.287

getmms_df$gene=row.names(getmms)
colnames(getmms_df)
getmms_df = getmms_df[,c(39,1:38)]
colnames(getmms_df)
#write_xlsx(getmms_df, 'getmms_ge5_1X.xlsx') # too big
write.table(getmms_df, 'getmms_ge5_1X_UNSTRANDED.txt', sep="\t")
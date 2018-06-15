#!/usr/bin/env Rscript
# Plot cumulative histogram
# Kenneth Condon
# Jan 2018

# clear environment
rm(list=ls())

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(plyr)
  library(reshape2)
})

# biomart
biomart <- read.delim(file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/ensembl91.txt", header = TRUE, dec = ".", fill = FALSE)
biomart$Gene.stable.ID <- as.character(biomart$Gene.stable.ID)

#zt3
zt3 <- read.delim(file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/p1_zt3_all.anno", header = FALSE, dec = ".", fill = TRUE)
colnames(zt3)[12] <- colnames(biomart)[1]
zt3$V20 <- NULL
zt3$time <- "zt3"
zt3$motif <- substring(zt3$V10, 15,16)
zt3<-subset(zt3, motif=="M1" | motif=="M2" | motif=="M3")
zt3.counts <- as.data.frame(count(zt3$Gene.stable.ID)) # requires plyr package
colnames(zt3.counts)[1] <- colnames(biomart)[1]
colnames(zt3.counts)[2] <- "zt3"
zt3.counts$Gene.stable.ID <- as.character(zt3.counts$Gene.stable.ID)

# zt15
zt15<- read.delim(file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/p1_zt15_all.anno", header = FALSE, dec = ".", fill = TRUE)
colnames(zt15)[12] <- colnames(biomart)[1]
zt15$V20 <- NULL
zt15$time <- "zt15"
zt15$motif <- substring(zt15$V10, 16,17)
zt15<-subset(zt15, motif=="M1" | motif=="M2" | motif=="M3")
zt15.counts <- as.data.frame(count(zt15$Gene.stable.ID)) # requires plyr package
colnames(zt15.counts)[1] <- colnames(biomart)[1]
colnames(zt15.counts)[2] <- "zt15"
zt15.counts$Gene.stable.ID <- as.character(zt15.counts$Gene.stable.ID)

# get counts and sort by difference between timepoints
all <- rbind(zt3, zt15)
temp <- full_join(zt3.counts,zt15.counts)
temp2 <- left_join(temp, biomart)
all.counts <- temp2[,c(1,4,2,3)]
suppressWarnings(all.counts[is.na(all.counts)] <- 0)
all.counts$diff<- abs((all.counts$zt3)-(all.counts$zt15))
all.counts.sorted <- all.counts[order(-all.counts$diff),]
write.table(all.counts.sorted, file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/all.counts.txt", sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
all.counts.sorted<-head(all.counts.sorted, n=60)

# using reshape2 to covert to long format
data.m <- melt(all.counts.sorted, 
                     variable.name = "variable",
                     value.names = "value",
                     id.vars = c("Gene.stable.ID", "Gene.name", "diff"))

# Plot histogram
target.hist<-ggplot(data.m, aes(x = reorder(Gene.name,-diff),y = value,fill=as.factor(variable)))+
  geom_bar(position="stack",stat="identity",colour = "black")+
  #geom_text(aes(label= value), size = 4, vjust = -0.2, hjust = 0.5)+
  labs(title= "Zfhx3 Target genes") +
  scale_y_continuous(limits = c(0,10), breaks = round(seq(min(0), max(10), by = 1),1), name = "Diff in number of binding sites")+
  theme_bw()+
  theme(
    legend.title=element_blank(),
    plot.title = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=12, angle=35, hjust = 1),#, face="bold"),   
    axis.text.y = element_text(size=12),
    panel.border = element_rect(colour="black",size=0.5),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect()
  )
target.hist

png("~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/targetHist.png", width = 1800, height = 800)
target.hist
dev.off()

# get all info about a named gene
mygeneName <- "Prcc"
mygene <- subset(all.counts.sorted, all.counts.sorted$Gene.name == mygeneName)
mygene <- subset(all, all$Gene.stable.ID == mygene$Gene.stable.ID)
mygene$Gene.name<-mygeneName
mygene[order(mygene$motif),]



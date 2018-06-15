#!/usr/bin/env Rscript
# Scatter plot of genomic locations
# Kenneth Condon
# Jan 2018

# clear environment
rm(list=ls())

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# store color blind friendly palette for plotting
cbPalette <- c("#D55E00", "#56B4E9", "#009E73")

# create scatter plot
filename = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/p1zt15.anno"
anno <- read.delim(file = filename, header = TRUE, dec = ".", fill = TRUE)
anno.filtered<-subset(anno, qval<0.01)
anno.loc<-subset(anno.filtered, select=c("seqnames", "start", "end", "width", "motif"))
anno.loc$mid<-abs((anno.loc$start)+((anno.loc$width)/2))
anno.min.mid<-signif(min(anno.loc$mid), 1)
anno.max.mid<-signif(max(anno.loc$mid), 1)
anno.loc$seqnames <- factor(anno.loc$seqnames, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM"))
plotScatter<-ggplot(anno.loc, aes(seqnames, mid, color=motif))+
  geom_point(size = 4, alpha = 0.4)+
  #labs(title= "anno Motif Distribution") +
  scale_y_continuous(limits = c(anno.min.mid,anno.max.mid), breaks = round(seq(min(anno.min.mid), max(anno.max.mid), by = 2e+07),1), name = "Motif Location")+
  #scale_y_continuous(limits = c(1.23e+08,1.43e+08), breaks = round(seq(min(1.23e+08), max(1.43e+08), by = 1e+06),1), name = "motif position")+
  scale_colour_manual(values=c(cbPalette[1],cbPalette[2],cbPalette[3]), name = "")+
  theme_bw()+
  theme(
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    #plot.title = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=30, face="bold"),
    axis.text.x = element_text(size=18, angle=35, hjust = 1, face="bold"),   
    axis.text.y = element_text(size=12),
    panel.border = element_rect(colour="black",size=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect()
  )
plotScatter

png(paste(filename,".scatter.png", sep=""), width = 1800, height = 800)
plotScatter
dev.off()


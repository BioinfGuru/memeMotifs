#!/usr/bin/env Rscript
# Plot venn diagrams
# Kenneth Condon
# Jan 2018

# clear environment
rm(list=ls())

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(VennDiagram)
  library(gplots)
})

# biomart
biomart <- read.delim(file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/ensembl91.txt", header = TRUE, dec = ".", fill = FALSE)
biomart$Gene.stable.ID <- as.character(biomart$Gene.stable.ID)

###################################################################################################
################                                ZT3                                ################
###################################################################################################
# prepare dataframe
zt3 <- read.table(text = system("cut -f 1-19 ~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/p1_zt3_all.anno | intersectBed -wa -wb -a - -b ~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/mESC_tads.mm10.txt", intern=TRUE))
zt3$motif <- substring(zt3$V10, 15,16)
colnames(zt3)[12] <- colnames(biomart)[1]
colnames(zt3)[23]<-"tad"
zt3$Gene.stable.ID <- as.character(zt3$Gene.stable.ID)
zt3<-subset(zt3, motif=="M1" | motif=="M2" | motif=="M3")
zt3 <- left_join(zt3, biomart)
zt3<-na.omit(zt3)

# plot venn of gene targets
zt3.targets<-list(unique(zt3[zt3$motif=="M1", 'Gene.name']),
                  unique(zt3[zt3$motif=="M2", 'Gene.name']),
                  unique(zt3[zt3$motif=="M3", 'Gene.name']))

class(zt3.targets)
length(zt3.targets)
class(zt3.targets[[1]])

zt3.targets.venn <- venn.diagram(zt3.targets,
                         NULL,
                         fill=c("darkmagenta", "darkblue", "darkgreen"),
                         alpha=c(0.5,0.5,0.5),
                         cex = 1.5,
                         lab.cex=5,
                         cat.fontface=4,
                         cat.pos=c(0,0,0),
                         category.names=c("ZT3.M1.targets", "ZT3.M2.targets", "ZT3.M3.targets"))
grid.newpage() # to reset image
grid.draw(zt3.targets.venn)

# Write targets to files
write.table(unique(zt3[zt3$motif=="M1", 'Gene.name']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M1.geneTargets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt3[zt3$motif=="M2", 'Gene.name']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M2.geneTargets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt3[zt3$motif=="M3", 'Gene.name']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M3.geneTargets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

# plot venn of TADs
zt3.tads<-list(unique(zt3[zt3$motif=="M1", 'tad']),
               unique(zt3[zt3$motif=="M2", 'tad']),
               unique(zt3[zt3$motif=="M3", 'tad']))
zt3.tads.venn <- venn.diagram(zt3.tads,
                                 NULL,
                                 fill=c("darkmagenta", "darkblue", "darkgreen"),
                                 alpha=c(0.5,0.5,0.5),
                                 cex = 1.5,
                                 lab.cex=5,
                                 cat.fontface=4,
                                 category.names=c("ZT3.M1.tads", "ZT3.M2.tads", "ZT3.M3.tads")
                                 #main="Gene Targets"
)
grid.newpage() # to reset image
grid.draw(zt3.tads.venn)

# Write tads to files
write.table(unique(zt3[zt3$motif=="M1", 'tad']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M1.tads.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt3[zt3$motif=="M2", 'tad']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M2.tads.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt3[zt3$motif=="M3", 'tad']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt3.M3.tads.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)


####################################################################################################
################                                ZT15                                ################
####################################################################################################
# prepare dataframe
zt15 <- read.table(text = system("cut -f 1-19 ~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/p1_zt15_all.anno | intersectBed -wa -wb -a - -b ~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/mESC_tads.mm10.txt", intern=TRUE))
zt15$motif <- substring(zt15$V10, 16,17)
colnames(zt15)[12] <- colnames(biomart)[1]
colnames(zt15)[23]<-"tad"
zt15$Gene.stable.ID <- as.character(zt15$Gene.stable.ID)
zt15<-subset(zt15, motif=="M1" | motif=="M2" | motif=="M3")
zt15 <- left_join(zt15, biomart)
zt15<-na.omit(zt15)

# plot venn of gene targets
zt15.targets<-list(unique(zt15[zt15$motif=="M1", 'Gene.name']),
                   unique(zt15[zt15$motif=="M2", 'Gene.name']),
                   unique(zt15[zt15$motif=="M3", 'Gene.name']))
zt15.venn <- venn.diagram(zt15.targets,
                         NULL,
                         fill=c("darkmagenta", "darkblue", "darkgreen"),
                         alpha=c(0.5,0.5,0.5),
                         cex = 1.5,
                         lab.cex=5,
                         cat.fontface=4,
                         category.names=c("ZT15.M1.targets", "ZT15.M2.targets", "ZT15.M3.targets"))
grid.newpage() # to reset image
grid.draw(zt15.venn)

# Write targets to files
write.table(unique(zt15[zt15$motif=="M1", 'Gene.name']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt15.M1.geneTargets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt15[zt15$motif=="M2", 'Gene.name']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt15.M2.geneTargets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

# plot venn of TADs
zt15.tads<-list(unique(zt15[zt15$motif=="M1", 'tad']),
                unique(zt15[zt15$motif=="M2", 'tad']),
                unique(zt15[zt15$motif=="M3", 'tad']))
zt15.tads.venn <- venn.diagram(zt15.tads,
                               NULL,
                               fill=c("darkmagenta", "darkblue", "darkgreen"),
                               alpha=c(0.5,0.5,0.5),
                               cex = 1.5,
                               lab.cex=5,
                               cat.fontface=4,
                               category.names=c("ZT15.M1.tads", "ZT15.M2.tads", "ZT15.M3.tads")
                               #main="Gene Targets"
)
grid.newpage() # to reset image
grid.draw(zt15.tads.venn)

# Write tads to files
write.table(unique(zt15[zt15$motif=="M1", 'tad']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt15.M1.tads.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(unique(zt15[zt15$motif=="M2", 'tad']), file = "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/annotations/zt15.M2.tads.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)


####################################################################################################
















#!/usr/bin/env Rscript
# Plot boxplot
# Kenneth Condon
# Jan 2018

# clear environment
rm(list=ls())

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# set working directory
wd <- "~/NGS/working_projects/Zfhx3_ChipSeq_v2/FASTQ_2_filtered/Macs_results/motifs/results/"

# read files
zt3 <- read.delim(file = paste(wd,"p1zt3.unfiltered.anno", sep=""), header = TRUE, dec = ".", fill = TRUE)
zt15 <- read.delim(file = paste(wd,"p1zt15.unfiltered.anno", sep=""), header = TRUE, dec = ".", fill = TRUE)

# choose motif
mymotif <- "M3"

# subset by motif
zt3<-subset(zt3, motif == mymotif)
zt15<-subset(zt15, motif == mymotif)

############
# STATISTICS
############

# create contingency table
zt3$hit <- 0
zt3$hit[zt3$qval <0.01] <- 1
zt15$hit <- 0
zt15$hit[zt15$qval<0.01] <- 1
zt3.t<-table(zt3$hit)
zt15.t<-table(zt15$hit)
all.t<-rbind(zt3.t,zt15.t)

# fisher test
result<-fisher.test(all.t)
pvalue<-signif(as.numeric(result$p.value),4)
pvalue.str<-paste("p = ",pvalue, sep="")
pvalue.str

#########
# BOXPLOT
#########

# transform qvalues
zt3$log <- -log(zt3$qval)
zt15$log <- -log(zt15$qval)

# create boxplot input
zt3$time <-"zt3"
zt15$time <- "zt15"
zt3<-subset(zt3, select=c("log", "time"))
zt15<-subset(zt15, select=c("log", "time"))
boxplot.in<-rbind(zt3,zt15)

# Create the boxplot
myBoxPlot<- ggplot(boxplot.in, aes(x= time, y = log, fill = time))+
  geom_hline(aes(yintercept=-log(0.01)), linetype="dashed", colour = "red")+
  stat_boxplot(geom ='errorbar')+
  geom_boxplot(outlier.size = 0)+
  guides(fill=FALSE)+
  scale_y_continuous(limits = c(0,18), breaks = round(seq(min(0), max(18), by = 2),1), name = "-log(qvalue)")+
  annotate("text", x=1.5, y=0, label = pvalue.str, size = 7)+
  theme_bw() +
  theme(
    plot.title = element_text(size=25,face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=18, angle=35, hjust = 1, face="bold"),  
    axis.text.y = element_text(size=18),
    panel.border = element_rect(colour="BLACK",size=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect()
  )

myBoxPlot

png(paste(wd,mymotif,"diff.box.png", sep=""), width = 600, height = 800)
myBoxPlot
dev.off()

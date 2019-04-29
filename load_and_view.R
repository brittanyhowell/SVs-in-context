## Usage: Plot mutliple annotation tracks given a set of tables
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 24th January 2019

#### Inputs:
# Coordinates to plot on chromosome
# Table of SVs
# Table of "annotations"
# Table of SNPs (Ideally this would have correlations)

#### Outputs: A beautiful plot with coordinates along the bottom, and separate tracks displaying geom_rects of the annotations


# setwd("~/Documents/scripts/brie_scripts/sv-view-multi-track/reference/ucsc/")
# ucsc.file <- read.table(file = "refSeq_funcElems.bed4", sep = "\t", stringsAsFactors = F)
# colnames(ucsc.file) <- c("chr", "start", "end", "description")
# 
# 
# segdup.file <- read.table(file = "seg_dups", sep = "\t", stringsAsFactors = F)
# segdup <- segdup.file[,1:4]
# colnames(segdup) <- c("chr", "start", "end", "description") 
# 
# sv.file <- read.table(file = "SVs.tmp", numerals="allow.loss", sep = "\t", stringsAsFactors = F, header = T)
# sv <- sv.file[,c(2,3,11,9)]
# colnames(sv) <- c("chr", "start", "end", "description") 
# 
# 
# sr.file <- read.table(file = "simple_repeats", sep = "\t", stringsAsFactors = F, header = T)
# sv <- sv.file[,c(2,3,11,9)]
# colnames(sv) <- c("chr", "start", "end", "description") 
# 
# 
# chroms <- paste("chr", c(1:22,"X","Y","M"), sep="")
# 
# fdat <- NULL
# for (i in 1:length(chroms))  {
#   print(c("filtering for chromosome", i))
#   indchr <- ucsc.file[,1]==chroms[i] 
#   fdat <- rbind(fdat, ucsc.file[indchr,])
# }
# ucsc <- fdat
# 
# fdat <- NULL
# for (i in 1:length(chroms))  {
#   print(c("filtering for chromosome", i))
#   indchr <- segdup[,1]==chroms[i] 
#   fdat <- rbind(fdat, segdup[indchr,])
# }
# segdup <- fdat
# 
# fdat <- NULL
# for (i in 1:length(chroms))  {
#   print(c("filtering for chromosome", i))
#   indchr <- sv[,1]==chroms[i] 
#   fdat <- rbind(fdat, sv[indchr,])
# }
# sv <- fdat
# 
# 
# 
# 
# ucsc.chr2 <- ucsc[which(ucsc$chr=="chr2" & ucsc$start>start.coord & ucsc$end<end.coord),]
# segdup.chr2 <- segdup[which(segdup$chr=="chr2" & segdup$start>start.coord & segdup$end<end.coord),]
# sv.chr2 <- sv[which(sv$chr=="chr2" & sv$start > start.coord & sv$end < end.coord),]
# 
# ucsc <- ucsc.chr2
# segdup <- segdup.chr2
# sv <- sv.chr2
# 
# start.coord <- 77850000
# end.coord <- 78000000


###### THIS IS THE STAGE WHERE THE ACTUAL PROGRAM WILL START

args = commandArgs(trailingOnly = TRUE)
library("ggplot2")#, lib.loc="/nfs/team151/software/Rlibs/")

sv <- read.table(file = args[1],  sep = "\t", stringsAsFactors = F, header = F)
ucsc <- read.table(file =args[2], sep = "\t", stringsAsFactors = F, header = F)
segdup <- read.table(file = args[3], sep = "\t", stringsAsFactors = F, header = F)
sr <- read.table(file = args[4], sep = "\t", stringsAsFactors = F, header = F)
dgv <- read.table(file = args[8], sep = "\t", stringsAsFactors = F, header = F)
clinv <- read.table(file = args[9], sep = "\t", stringsAsFactors = F, header = F)

# which( chrB

sv <- sv[which(sv[,1]!="chrB"),]
ucsc <- ucsc[which(ucsc[,1]!="chrB"),]
segdup <- segdup[which(segdup[,1]!="chrB"),]
sr <- sr[which(sr[,1]!="chrB"),]
dgv <- dgv[which(dgv[,1]!="chrB"),]
clinv <- clinv[which(clinv[,1]!="chrB"),]

start.coord <- as.integer(args[5])
end.coord <- as.integer(args[6])



## Read in other tracks
## Add levels to plots
## Place in groups and draw lines in between them
## 

# size = 8, aes(x= ALT, y= num, label = num),

png(filename = paste(args[7],paste("sv",start.coord,end.coord,"png",sep="."),sep = "/"), width = 700, height = 500)


p <- ggplot() +
 geom_rect(aes(xmin=sv[,2], xmax=sv[,3]+1000, ymin=0,ymax=1, fill=sv[,4] )) +
  
  scale_x_continuous(limits= c(start.coord-50000,end.coord+10000)) +
  # facet_grid(segdup$chr ~.)  +
  theme_bw() +
  theme(strip.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12)
        ) + 
  xlab(paste("\nchromosome coordinates -",sv[,1],sep=" "))  + 
  guides(fill=guide_legend(title="annotations:")) +
  geom_text(size=6,aes(x=start.coord-50000,y=10.5,label="clinvar"), hjust = 0)  +
  geom_text(size=6,aes(x=start.coord-50000,y=8.5,label="DGV"), hjust = 0)  +
  geom_text(size=6,aes(x=start.coord-50000,y=6.5,label="simple_repeats"), hjust = 0)  +
  geom_text(size=6,aes(x=start.coord-50000,y=4.5,label="func_elem"), hjust = 0)  +
  geom_text(size=6,aes(x=start.coord-50000,y=2.5,label="seg_dup"), hjust = 0)  +
  geom_text(size=6,aes(x=start.coord-50000,y=0.5,label="SV"), hjust = 0)  +
  scale_fill_brewer(palette="Dark2")

  
if (nrow(ucsc) > 0 ) {
  p <- p +  geom_rect(aes(xmin=ucsc[,2], xmax=ucsc[,3], ymin=4,ymax=5, fill=ucsc[,4] ))
  }
if (nrow(segdup) > 0 ) {
  p <- p + geom_rect(alpha=0.8, aes(xmin=segdup[,2], xmax=segdup[,3], ymin=2,ymax=3), fill="black")
  }
if (nrow(sr) > 0 ) {
  p <- p + geom_rect(aes(xmin=sr[,2], xmax=sr[,3], ymin=6,ymax=7 ), fill="grey")
  }
if (nrow(dgv) > 0 ) {
  p <- p + geom_rect(aes(xmin=dgv[,2], xmax=dgv[,3], ymin=8,ymax=9, fill=dgv[,4] )) 
  }
 
if (nrow(clinv) > 0 ) {
  p <- p + geom_rect(aes(xmin=clinv[,2]-100, xmax=clinv[,3]+100, ymin=10,ymax=11)) 
  p <- p + geom_text(position=position_jitter(width=50000, height=1),size=4,aes(x=clinv[,2],y=10,label=clinv[,4]), hjust = 0)  
  }
print(p)
  # 
  #  +
  # +
  
  # geom_hline(aes(yintercept=1.5)) +
  

dev.off()
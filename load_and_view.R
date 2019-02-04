## Usage: Plot mutliple annotation tracks given a set of tables
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 24th January 2018

#### Inputs:
# Coordinates to plot on chromosome
# Table of SVs
# Table of "annotations"
# Table of SNPs (Ideally this would have correlations)

#### Outputs: A beautiful plot with coordinates along the bottom, and separate tracks displaying geom_rects of the annotations

setwd("~/Documents/scripts/brie_scripts/sv-view-multi-track/reference/ucsc/")
ucsc.file <- read.table(file = "refSeq_funcElems.bed4", sep = "\t", stringsAsFactors = F)
colnames(ucsc.file) <- c("chr", "start", "end", "description")

# ucsc$chr <- as.character(ucsc$chr)
# ucsc <- head(ucsc)

segdup.file <- read.table(file = "seg_dups", sep = "\t", stringsAsFactors = F)
segdup <- segdup.file[,1:4]
colnames(segdup) <- c("chr", "start", "end", "description") 

sv.file <- read.table(file = "../../../../../data/sv-analyse/QC/sorted_all.txt", sep = "\t", stringsAsFactors = F, header = T)
sv <- sv.file[,c(2,3,11,9)]
colnames(sv) <- c("chr", "start", "end", "description") 


chroms <- paste("chr", c(1:22,"X","Y","M"), sep="")

fdat <- NULL
for (i in 1:length(chroms))  {
  print(c("filtering for chromosome", i))
  indchr <- ucsc.file[,1]==chroms[i] 
  fdat <- rbind(fdat, ucsc.file[indchr,])
}
ucsc <- fdat

fdat <- NULL
for (i in 1:length(chroms))  {
  print(c("filtering for chromosome", i))
  indchr <- segdup[,1]==chroms[i] 
  fdat <- rbind(fdat, segdup[indchr,])
}
segdup <- fdat

fdat <- NULL
for (i in 1:length(chroms))  {
  print(c("filtering for chromosome", i))
  indchr <- sv[,1]==chroms[i] 
  fdat <- rbind(fdat, sv[indchr,])
}
sv <- fdat




ucsc.chr2 <- ucsc[which(ucsc$chr=="chr2"),]
segdup.chr2 <- segdup[which(segdup$chr=="chr2"),]
sv.chr2 <- sv[which(sv$chr=="chr2"),]

start.coord <- 77850000
end.coord <- 78000000


###### THIS IS THE STAGE WHERE THE ACTUAL PROGRAM WILL START

ucsc <- ucsc.chr2[which(ucsc.chr2$start>start.coord & ucsc.chr2$start<end.coord ),]
segdup <- segdup.chr2[which(segdup.chr2$start>start.coord & segdup.chr2$start<end.coord),]
sv <- sv.chr2[which(sv.chr2$start>start.coord & sv.chr2$start<end.coord),]


## Read in other tracks
## Add levels to plots
## Place in groups and draw lines in between them
## 

# size = 8, aes(x= ALT, y= num, label = num),

ggplot() +
  geom_rect(aes(xmin=ucsc$start, xmax=ucsc$end, ymin=4,ymax=5, fill=ucsc$description )) +
  # geom_point(aes(x=segdup$start,y=1.5))  +
  geom_hline(aes(yintercept=1.5)) +
  geom_rect(alpha=0.8, aes(xmin=segdup$start, xmax=segdup$end, ymin=2,ymax=3), fill="black") +
  geom_rect(aes(xmin=sv$start, xmax=sv$end, ymin=0,ymax=1, fill=sv$description)) +
  scale_x_continuous(limits= c(start.coord,end.coord)) +
  # facet_grid(segdup$chr ~.)  +
  theme_bw() +
  theme(strip.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12)
        ) + 
  xlab(paste("\nchromosome coordinates -",sv$chr,sep=" "))  + 
  guides(fill=guide_legend(title="refSeq functional element:")) +
  geom_text(size=4,aes(x=start.coord,y=4.5,label="UCSC"))  +
  geom_text(size=4,aes(x=start.coord,y=2.5,label="SegDup"))  +
  # geom_text(size=4,aes(x=start.coord,y=4.5,label="UCSC"))  +
  geom_text(size=4,aes(x=start.coord,y=0.5,label="SV"))  +
  scale_fill_brewer(palette="Dark2")


############################################################################### 
#
# COMETgazer for methylation analysis
# 
# By: Emanuele Libertini & Rifat Hamoudi
#
# COMETvintage : template R script for DMC analysis using edgeR
#
# Input : count distributions over genomic windows of 100,000 bp in size (one for each sample)
#
# Output : text file with the coordinates of DMCs
#          chr | start | end 
#
################################################################################


# Count distributions (output of OORTcloud) for test samples need to be read in
# in this example these are represented by the R objects: sample1_r1, sample1_r2, sample1_r3, sample1_r4, sample2_r1, sample2_r2


library(Repitools)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
library(rtracklayer)
library(edgeR)

genome.windows.100K <- genomeBlocks(Hsapiens, chrs=paste("chr", c(1:22), sep=""), width=100000)

group <- factor(c(1,1,1,1,2,2))

# Coefficient of variation : this is a user-defined variable, it can be estimated by edgeR
# hCOMETs have typically higher dispersion 
bcv = 0.4 	 

counts.low <- data.frame(sample1_r1, sample1_r2, sample1_r3, sample1_r4, sample2_r1, sample2_r2)
counts.low <- subset(counts.low, rowSums(counts.low) >  4) # this is a also user-defined parameter; lCOMETs have typically low counts

y.low.100K <- DGEList(counts=counts.low, group=group)
et.low <- exactTest(y.low.100K, dispersion=bcv^2)
summary(de.low <- decideTestsDGE(et.low, p=0.05, adjust="BH"))

detags.low <- rownames(y.low.100K)[as.logical(de.low)]
DMC.low <- as.data.frame(genome.windows.100K)[detags.low,]

png("plotSmear.png")
plotSmear(et.low, de.tags=detags.low)
abline(h = c(-2, 2), col = "blue")  
dev.off()

write.table(summary(de.low <- decideTestsDGE(et.low, p=0.05, adjust="BH")), "summary.low.methylation.DMB.100K.txt", sep="\t", quote=F)                                                                         
write.table(length(detags.low), "detags.low.methylation.DMB.100K.txt", sep="\t", quote=F)

DMC.low <- as.data.frame(genome.windows.100K)[detags.low,]

write.table(DMC.low, "DMC.in.lCOMETs.txt", sep="\t", quote=F, row.names=F, col.names=F)



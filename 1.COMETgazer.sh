#!/bin/bash

############################################################################ 
#
# COMETgazer suite for methylation analysis
# 
# By: Emanuele Libertini & Rifat Hamoudi
#
# COMETgazer : methylome segmentation into COMETs based on OM values
#
# Input : whole genome bisulfite sequencing methylation data that have been smoothed 
#         The input file should have the following format:
#
#         Chr No | Begin | End | Smoothed methylation value
#
# Output : block | chr | start | stop | meth | oscillator | rounded OM value | size | max | min | average | range
#
############################################################################


# OM estimation and assignment of CpGs to COMETs

for chrom in `seq 1 22` 
do
	R --no-save --verbose <<EOF

		library(quantmod)
		library(TTR) 	# This calls the Technical Trading Rules library in R which enhances quantmod library

		chr$chrom.meth.df <- read.table("chr$chrom.txt", header=F)

		del1 <- Delt(chr$chrom.meth.df[,4])
		chr$chrom.df <- data.frame(chr$chrom.meth.df, del1)
		chr$chrom.df[1,5] <- 0
		g <- factor(round(chr$chrom.df[,5], 1))
		chr$chrom.df[,6] <- g
		names(chr$chrom.df)[6] <- "g"
		write.table(chr$chrom.df, "chr$chrom.df.txt", sep="\t", row.names=F, col.names=F, quote=F)
EOF

	sed 's/chr//g'  chr$chrom.df.txt > chr.df.txt
	touch output.chr.df.txt

	./blocks	# Calls a program written in C++ for block segmentation/assignment

	mv output.chr.df.txt chr$chrom.blocks.df.txt  # File containing the oscillator of methylation (OM) distribution for every CpG 
done



# Reducing and summarising the data to COMET coordinates

for chrom in `seq 1 22` 
do
	awk '($6 != 0) { print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' chr$chrom.blocks.df.txt > chr$chrom.run.10.percent.blocks.txt
	awk 'NR>1' chr$chrom.run.10.percent.blocks.txt > chr$chrom.blocks.10.txt
	awk 'NR==2' chr$chrom.blocks.df.txt > chr$chrom.blocks.head
	cat chr$chrom.blocks.head chr$chrom.blocks.10.txt > chr$chrom.blocks.10.final.txt

	awk -F"\t" '
	{
    		sum[$7]+=$4;
    		if(!($7 in min) || (min[$7]>$4))
        		min[$7]=$4;
    		if(!($7 in max) || (max[$7]<$4))
        		max[$7]=$4;
    		count[$7]++
	}
	END 
	{
		for(element in sum)
      			printf("%s, %01f, %01f, %01f\n", element, max[element], min[element], sum[element]/count[element])
	}' chr$chrom.blocks.df.txt > chr$chrom.block.sizes.txt
	

	sort -gk 1 chr$chrom.block.sizes.txt > chr$chrom.block.sizes.sorted.txt

	R --no-save --verbose <<EOF
		library(GenomicFeatures)
		chr$chrom.block.sizes <- read.csv("chr$chrom.block.sizes.sorted.txt", header=F, skip=1)
		chr$chrom.block.sizes[,5] <- (chr$chrom.block.sizes[,2] - chr$chrom.block.sizes[,3])/2*100
		names(chr$chrom.block.sizes) <- c("block", "max", "min", "average", "range")
		write.table(as.matrix(summary(chr$chrom.block.sizes[,5])), "chr$chrom.summary.block.sizes.range.txt", sep="\t")

		chr$chrom.blocks <- read.table("chr$chrom.blocks.10.final.txt", header=F)
		names(chr$chrom.blocks) <- c("chr", "start", "stop", "meth", "oscillator", "g", "block")
		chr$chrom.tail <- read.table(text=(system("tail -1 chr$chrom.blocks.df.txt", intern = TRUE)), sep = "\t")

		chr$chrom.blocks[1:dim(chr$chrom.blocks)[1]-1,3] <- chr$chrom.blocks[2:dim(chr$chrom.blocks)[1],2]
		chr$chrom.blocks[dim(chr$chrom.blocks)[1],3] <- chr$chrom.tail[2]
		chr$chrom.blocks[,8] <- chr$chrom.blocks[,3] - chr$chrom.blocks[,2]
		names(chr$chrom.blocks)[8] <- "size"

		chr$chrom.blocks.total <- merge(chr$chrom.blocks, chr$chrom.block.sizes, by ="block")

		chr$chrom.blocks.head <- subset(chr$chrom.blocks.total, chr$chrom.blocks[,8] > 300 )
		chr$chrom.blocks.trimmed <- subset(chr$chrom.blocks.head, chr$chrom.blocks.head[,8] < quantile(chr$chrom.blocks.head[,8], 0.9999))
		write.table(chr$chrom.blocks.trimmed, "chr$chrom.blocks.verified.txt", sep="\t", row.names=F)
		write.table(as.matrix(summary(chr$chrom.blocks.trimmed[,8])), "chr$chrom.block.sizes.summary.txt", sep="\t")
EOF
done


# I/O cleanup

rm *.head
rm *run*
rm *blocks.10*
rm *sizes*
mkdir COMETs
mv *verified* COMETs
mkdir data.frame
mkdir summaries
mv *summary* summaries/
mv chr*txt data.frame/



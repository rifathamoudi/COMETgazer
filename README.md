
COMETgazer suite
-----------------


COMETgazer mehylation analysis software suite consists of 3 software:


1) COMETgazer : bash script for methylome segmentation into COMETs based on OM values

Input			: whole genome bisulfite sequencing methylation data that have been smoothed

Part of the bash script is blocks.cpp which is a C++ program assigning CpGs to COMETs


Example Input file	: chr22.txt as an example of Chromosome 22 only. The software loops through all 22 chromosomes
Example Output file	: chr22_blocks_verified.txt example of COMET segmentation of Chromosome 22 only



2) OORTcloud : bash script for counting COMET distributions according to methylation level

Input			: this will be the output from COMETgazer e.g chr22_blocks_verified.txt. The software loops through all 22 chromosomes

Example Input file	: chr22_blocks_verified.txt example of COMET segmentation of Chromosome 22 only. The software loops through all 22 chromosomes
Example Output file	: low.txt low methylation level COMET distribution genome wide
			  medium.txt medium methylation level COMET distribution genome wide
			  high.txt high methylation level COMET distribution genome wide



3) COMETvintage : template R script for DMC analysis using edgeR

Example Input file	: low.txt low methylation level COMET distribution genome wide for one methylome
			  medium.txt medium methylation level COMET distribution genome wide for one methylome
			  high.txt high methylation level COMET distribution genome wide for one methylome
Example Output file	: text file with the coordinates of DMCs


#Create tdf files from a list of bams (no header):

library(methods)
source("/xchip/cga_home/mara/scripts/functions_mara.R")

args <- commandArgs(trailingOnly=TRUE)

bam_list.file <- args[1]

bam_list = read.table(bam_list.file, sep="\t", header=F, quote='', stringsAsFactors=F)

for (i in seq(1, dim(bam_list)[1])) {
	cmd = paste("use igvtools_2.3; igvtools count -w 5000000 ", bam_list$V1[i]," ",bam_list$V1[i], ".tdf hg19", sep="")
	bcmd = bsub_cmd(cmd, jname=paste(bam_list$V1[i], ".tdf", sep=""), queue="hour")
	system(bcmd)
}

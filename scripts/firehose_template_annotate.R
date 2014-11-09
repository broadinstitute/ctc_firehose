#################################
#CTC: QC Analysis
#Match walk up seq bam files to the upload table provided for firehose upload
#8 Oct 2014
#################################

#input: 
#[1] a list of fill path bam files
#[2] the template for filling in with full bam file paths. 
#[3] sample set up for upload

#source library
library(Rsamtools)

#load in files
args <- commandArgs(trailingOnly=TRUE)
bams.file = args[1]
template.file = args[2]
sset_id = args[3]

print("Loading template and bam paths...")
template = read.table(template.file, header=T, quote='', stringsAsFactors=F, sep="\t")
bams = read.table(bams.file, header=F, quote='', stringsAsFactors=F, sep="\t")

print("Adding bam files to template...")
for (i in seq(1,dim(bams)[1])) {
	bam.path = bams[i,"V1"]
	a = scanBamHeader(bam.path)
	sm = a[[1]][[2]]$`@RG`[grep("SM:",a[[1]][[2]]$`@RG`)]
	samp = substr(sm,start=4, stop=nchar(sm))
	if (samp %in% template$sample_id) {
		print(paste("Adding this sample:",samp))
		template[template$sample_id %in% samp,"bam_file_lowpass"] <- bam.path
	} else {
		print(paste("This sample does not exist in the master template:", samp))
	}
}

print("Writting out files...")
#write out full template file even if missing things:
write.table(template,paste(template.file, ".bam_lowpass.annotated", sep=""), sep="\t", quote=F, row.names=F)
#write out file of missing samples still:
temp_miss = subset(template, is.na(bam_file_lowpass))
write.table(temp_miss,paste(template.file, ".bam_lowpass.annotated", sep=""), sep="\t", quote=F, row.names=F)
#write out current firehose upload file:
temp_complete = subset(template, !(is.na(bam_file_lowpass)))
write.table(temp_complete,"upload_lowpass_bams.txt", sep="\t", quote=F, row.names=F)

sample_set_id = rep(sset_id, 1, dim(temp_complete)[1])
temp_sset = cbind(temp_complete$sample_id, sample_set_id)
colnames(temp_sset) <- c("sample_id", "sample_set_id")
write.table(temp_sset, "upload_sset.txt", sep="\t", quote=F, row.names=F)

print("Done!")
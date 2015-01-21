##########
#CTC: Firehose combine metrics across samples
#MR
#Nov 17,2014
##############

args <- commandArgs(trailingOnly=TRUE)
hybrid_files.file <- args[1]
corr_coeff.file <- args[2]
sset <- args[3]
#wig_xls.file <- args[4]	

#to add: correlation length, % contamination from bacteria? - chandra. 

#libraries:

#testing:
#hybrid_files.file = "/xchip/cga_home/mara/projects/ctc_sigma/bam_file_lowpass_metrics.txt"
#corr_coeff.file = "/xchip/cga_home/mara/projects/ctc_sigma/corrCoeff.txt"
#wgs_xls.file = "/xchip/cga_home/mara/projects/ctc_sigma/lowpass_21Oct2014.xls"
#sset = "lowpass_21Oct2014"

print("reading in firehose hybrid metrics annotation file..")
hybrid_files = read.table(hybrid_files.file, sep="\t", header=F, quote='', stringsAsFactors=F)
corr_coeff = read.table(corr_coeff.file, sep="\t", header=F, quote='', stringsAsFactors=F)
#wgs_xls_workbook = loadWorkbook(wgs_xls.file)
#wgs_xls_wgs = readWorksheet(wgs_xls_workbook, sheet="Coverage", header=T)

print("looping through metrics file...")
final = data.frame()
all_names = NULL
count = 0
for (i in seq(1,dim(hybrid_files)[1])) {
	print(paste("reading:",hybrid_files[i,"V2"]))
	a = t(read.table(hybrid_files[i,"V2"], sep="\t", header=T, comment="#", quote='', stringsAsFactors=F))
	name_file = hybrid_files[i,"V1"]
	#colnames(a) = name_file
	CORR_COEFF = corr_coeff[corr_coeff[,1] %in% hybrid_files[i,"V1"],"V2"]
	a = rbind(CORR_COEFF,a)[,3]
	if (count %in% 0) {
		final = a
		all_names = c(all_names,name_file)
	} else {
		final = cbind(final, a)
		all_names = c(all_names,name_file)
	}
	count = count + 1
}
colnames(final) = all_names

print("writting out file...")
write.table(final, file=paste(sset, ".metrics.txt", sep="\t", quote=F, quoate=F, row.names=F)
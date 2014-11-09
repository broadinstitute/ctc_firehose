#################
#CTC Firehose
#Generate a master table (ignore power calculations - just raw maf input)
#
#MR Aug 11, 2014
##################

#Load Arguments:
args <- commandArgs(trailingOnly=TRUE)
libdir <- args[1]
maf_files.file <- args[2] #list of maf files with column A as pair_id
sample_type.file <- args[3] #list of the cell_type annotation (WBC, CTC, TP, TM, TR)
pair_set_id <- args[4] #for naming of output file
number_overlap <- args[5] #number of cells to overlap for calling of mutation

maf_files = read.table(maf_files.file, sep="\t", quote='', stringsAsFactors=F, header=F)
sample_type = read.table(sample_type.file, sep="\t", quote='', stringsAsFactors=F, header=F)
colnames(maf_files) <- c("pair_id", "maf_file")
colnames(sample_type) <- c("sample_id", "cell_type")
#colnames(map_sample) <- c("sample_id", "case_sample")

#Source necessary libraries:
source(paste(libdir, "heatmap_functions_v3.R", sep=""))

#Load Data:
data_file_maf = paste(pair_set_id, ".RData", sep="")
print("generating maf file combined file...")
all_maf = merge_mutation_file(maf_files, "maf")
save(file=data_file_maf, all_maf)

#Fancy mapping set up:
map = all_maf[,c("Tumor_Sample_Barcode", "pair_id")]
colnames(map) = c("sample_id","pair_id")

combo = merge(maf_files, map, by="pair_id")
combo = merge(combo,sample_type, by="sample_id")

#Rearrange Mutations so only include >=X CTCs and >=X WBCs
ctc_samples = subset(combo, cell_type %in% "CTC")
wbc_samples = subset(combo, cell_type %in% "WBC")
bulk_samples = subset(combo, !(cell_type %in% "CTC") & !(cell_type %in% "WBC"))

if (!(dim(ctc_samples)[1] %in% 0)) {
	ctc = 1
	print("calculating ctc maf")
	ctc_maf = subset(all_maf, Tumor_Sample_Barcode %in% ctc_samples$sample_id)
	print("extracting sites to consider...")
	ctc_stats = data.frame(table(ctc_maf$master_key))
	ctc_sites = as.character(subset(ctc_stats, Freq>=number_overlap)[,"Var1"])
} else {
	ctc = 0
	print("no ctcs")
	ctc_sites = NULL
}

if (!(dim(wbc_samples)[1] %in% 0)) {
	wbc = 1
	print("calculating wbc maf")
	wbc_maf = subset(all_maf, Tumor_Sample_Barcode %in% wbc_samples$sample_id)
	print("extracting sites to consider...")
	wbc_stats = data.frame(table(wbc_maf$master_key))
	wbc_sites = as.character(subset(wbc_stats, Freq>=number_overlap)[,"Var1"])
} else {
	wbc = 0
	print("no wbcs")
	wbc_sites = NULL
}

if (!(dim(bulk_samples)[1] %in% 0)) {
	print("calculating bulk maf")
	bulk_maf = subset(all_maf, Tumor_Sample_Barcode %in% bulk_samples$sample_id)
	print("extracting sites to consider...")
	bulk_sites = unique(bulk_maf$master_key)
} else {
	print("no bulk")
	bulk_sites = NULL
}

#Sites to Consider in Master Diagram:
all_sites = unique(c(ctc_sites, wbc_sites, bulk_sites))
if (length(all_sites) %in% 0) {
	print("ERROR: No mutation sites for output master table. Check input maf files.")
}

#Generate binary table of called or not:
print("generating binary table of called or not...")
all_samples = unique(all_maf$Tumor_Sample_Barcode)
master = matrix(0, length(all_sites), length(all_samples))
rownames(master) <- as.character(all_sites)
colnames(master) = as.character(all_samples)
for (i in seq(1,length(all_samples))) {
	print(paste("adding", all_samples[i], "to master table..."))
	single_maf = subset(all_maf, Tumor_Sample_Barcode %in% all_samples[i])
	master[intersect(rownames(master), single_maf$master_key),all_samples[i]]=1
}

master = data.frame(master)
master$key = rownames(master)

#Add to binary an all_ctc column
if (ctc) {
	print("adding combined ctcs to master binary table...")
	master$ctc_combined <- 0
	master$ctc_combined[master$key %in% ctc_sites] = 1
}

if (wbc) {
	print("adding combined wbcs to master binary table...")
	master$wbc_combined <-0
	master$wbc_combined[master$key %in% wbc_sites] = 1
}

#Adding useful mutation columns from key column
print("master made, turning key into columns...")
tmp = data.frame(do.call(rbind, strsplit(master$key, split=":")))
colnames(tmp) = c("Hugo_Symbol", "Chromosome", "Start_position", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
master = cbind(master, tmp)

#Write out table:
print("writing out master table...")
write.table(master, file=paste(pair_set_id, ".combined_mutation_table.txt", sep=""), quote=F, sep="\t", row.names=F)

print("writing out interval file...")
interval_file = paste(tmp$Chromosome, paste(tmp$Start_position, tmp$Start_position, sep="-"), sep=":")
write.table(interval_file, file=paste(pair_set_id, ".forcecall.intervals", sep=""), quote=F, sep="\t", row.names=F, col.names=F)

print("done!")

##############################
#TESTING:
#Rscript fh_combine_mutations.R /xchip/cga_home/mara/projects/ctc_sigma/scripts/ /xchip/cga_home/mara/projects/ctc_sigma/debug/maf_file.input.tsv /xchip/cga_home/mara/projects/ctc_sigma/debug/tumor_cell_type.input.tsv 10448EV-TEST 3

#testing from session:
#libdir = "/xchip/cga_home/mara/projects/ctc_sigma/scripts/"
#maf_files.file = "/xchip/cga_home/mara/projects/ctc_sigma/debug/maf_file.input.tsv"
#sample_type.file = "/xchip/cga_home/mara/projects/ctc_sigma/debug/tumor_cell_type.input.tsv"
#pair_set_id = "10448EV-TEST"
#number_overlap = 3
#############################
#################
#CTC Firehose
#Power Function
#Generate Table with Counts for Depth for all sites not called. (Also, generate one with arbitrary input cutoff.)
#
#MR Aug 13, 2014
##################

args <- commandArgs(trailingOnly=TRUE)
libdir <- args[1]
master.file <- args[2] #list of maf files with column A as pair_id
call_stats.file <- args[3] # list of force called files (to use for power)
#sample_type.file <- args[4] #list of the cell_type annotation (WBC, CTC, TP, TM, TR)
pair_set_id <- args[4] #for naming of output file
#depth_cutoff <- args[5] #minimum depth of coverage required - not using right now

###################
#Load in files:
print("loading in files...")
master = read.table(master.file, sep="\t", header=T, quote='', stringsAsFactors=F)
call_stats = read.table(call_stats.file, sep="\t", header=T, quote='', stringsAsFactors=F)
colnames(call_stats) <- c("pair_id", "call_stats_file")
master$key = paste(master$Chromosome, master$Start_position, sep=":")
#Source necessary libraries:
#source(paste(libdir, "heatmap_functions_v4.R", sep=""))

#Generate sample to pair map:
print("generate sample to pair map...")
pair_map = data.frame(stringsAsFactors=F)
for (i in seq(1,dim(call_stats)[1])) {
	cs.file = as.character(call_stats$call_stats_file[i])
	if (file.exists(cs.file)) {
		header = read.table(header=FALSE, nrow=1,file=cs.file, sep="\t")
		take <- c("tumor_name")
  		takecols <- ifelse(t(header) %in% take, NA, 'NULL')
		tname <- unique(read.table(header = TRUE, colClasses=takecols, file=cs.file, sep="\t", quote='', stringsAsFactors=FALSE, nrows=3))
		tmp = t(data.frame(c(call_stats$pair_id[i],as.character(tname))))
		pair_map = rbind(pair_map, tmp)
	}
}
colnames(pair_map) = c("pair_id","tumor_name")

#Add depth to master file
print("merging in tumor depth on each sample...")
final = master[,c("key","Hugo_Symbol","Chromosome","Start_position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")]
colna=c("key","Hugo_Symbol","Chromosome","Start_position","Variant_Classification","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")
for (i in seq(1,dim(pair_map)[1])) {
	if (pair_map[i,"tumor_name"] %in% colnames(master)) {
		print(as.character(pair_map[i,"tumor_name"]))
		colna = c(colna,as.character(pair_map[i,"tumor_name"]))
		cs.file = as.character(call_stats[pair_map[i,"pair_id"],"call_stats_file"])
		cs = read.table(cs.file, sep="\t", header=T, quote='', stringsAsFactors=F)
		cs$key = paste(cs$contig, cs$position, sep=":")
		cs$depth = paste(cs$t_alt_count + cs$t_ref_count)
		final = merge(final, cs[,c("key","depth")], all.x=T, by="key")
		colnames(final) = colna
	} else {
		print("this sample not in master table:", pair_map[i,"tumor_name"])
	}
}

#Add in called information:
print("add whether or not called to master table...")
for (i in seq(1,dim(pair_map)[1])) {
	if (pair_map[i,"tumor_name"] %in% colnames(master)) {
		pair = as.character(pair_map[i,"tumor_name"])
		print(pair)
		sites_called = master[master[,pair]%in% 1,"key"]
		final[final$key %in% sites_called,pair] = "called"
	} else {
		print("this sample not in master table:", pair_map[i,"tumor_name"])
	}
}

#Write out table:
print("write out final table...")
write.table(final, file=paste(pair_set_id, ".combined_mutation_table.depth.txt", sep=""), sep="\t", quote=F, row.names=F)
print("done!")

#########################
#TEST:
#Rscript fh_annotate_master_with_depth.R /xchip/cga_home/mara/projects/ctc_sigma/scripts/ /xchip/cga/gdac-prod/cga/jobResults/CTC_MasterMafIntervals/10448VPM-All/10379049/10448VPM-All.combined_mutation_table.txt /xchip/cga_home/mara/projects/ctc_sigma/debug/10448VPM.call_stats.txt 10448VPM





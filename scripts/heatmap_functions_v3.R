#Functions for Generating Heatmap:

merge_mutation_file <- function(maf_info, type) {
	#Merge Maf files based on the following four main columns
	combined_maf = NULL
	maf_column = maf_info[,"maf_file"]

	if (type=="maf") {
		for (i in seq(1,length(maf_column)[1])) {
			if (file.exists(maf_column[i])) {
			    header = read.table(header=FALSE, nrow=1,file=as.character(maf_column[i]), sep="\t", stringsAsFactors=FALSE)
			    take <- c("Hugo_Symbol", "Chromosome", "Start_position", "Tumor_Sample_Barcode", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
			    takecols <- ifelse(t(header) %in% take, NA, 'NULL')
			    maf <- read.table(header = TRUE, colClasses=takecols, file=as.character(maf_column[i]), sep="\t", quote='')
			   	maf$key = paste(maf$Chromosome, maf$Start_position, sep=":")
			   	maf$master_key = paste(maf$Hugo_Symbol, maf$Chromosome, maf$Start_position, maf$Variant_Classification, maf$Reference_Allele, maf$Tumor_Seq_Allele1, maf$Tumor_Seq_Allele2, sep=":")
		    	maf$pair_id = maf_info$pair_id[i]
		    }
		    if (dim(maf)[1]==0) {
		        print("maf file empty")
		    } else {
		    	combined_maf = rbind(combined_maf, maf)
		    }
		}
	}
	if (type=="absolute") {
		for (i in seq(1,length(maf_column)[1])) {
		    header = read.table(header=FALSE, nrow=1,file=as.character(maf_column[i]), sep="\t")
		    take <- c("Hugo_Symbol","Chromosome", "Start_position", "Tumor_Sample_UUID", "ref", "alt", "Variant_Classification", "i_judgement", "cancer_cell_frac", "Tumor_Sample_Barcode")
			takecols <- ifelse(t(header) %in% take, NA, 'NULL')
		    maf <- read.table(header = TRUE, colClasses=takecols, file=as.character(maf_column[i]), sep="\t", quote='', stringsAsFactors=FALSE)
		   	maf$key = paste(maf$Chromosome, maf$Start_position, sep=":")
		   	maf$master_key = paste(maf$Hugo_Symbol, maf$Chromosome, maf$Start_position, maf$Variant_Classification, sep=":")
		    if (dim(maf)[1]==0) {
		        print("maf file empty")
		    } else {
		    	combined_maf = rbind(combined_maf, maf)
		    }
		}
	}
	if (type=="callstats") {
		for (i in seq(1,length(maf_column)[1])) {
		    header = read.table(header=FALSE, nrow=1,file=as.character(maf_column[i]), sep="\t")
		    take <- c("contig","position","tumor_name", "judgement", "t_alt_count", "t_ref_count")
  			takecols <- ifelse(t(header) %in% take, NA, 'NULL')
		    maf <- read.table(header = TRUE, colClasses=takecols, file=as.character(maf_column[i]), sep="\t", quote='', stringsAsFactors=FALSE)
		   	maf$key = paste(maf$contig, maf$position, sep=":")
		   	print(unique(as.character(maf$tumor_name)))
		    if (dim(maf)[1]==0) {
		        print("maf file empty")
		    } else {
		    	combined_maf = rbind(combined_maf, maf)
		    }
		}
	}
	return(combined_maf)
}	

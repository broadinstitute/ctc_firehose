############# PUT INTO FIREHOSE (once add fraction dropped sites and run as pairset too)
# MR
# CTC Allelic Distortion (run on a pairset)
# INPUT: List of Callstats Files (het pull down of normal), and pair_type_ctc
##################

source("functions.R")
library(RColorBrewer)
library(TeachingDemos)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

callstats.file <- args[1]
pair_type.file <- args[2]

#callstats.file = "/xchip/cga_home/mara/projects/ctc_sigma/analysis_raw_bam_files/call_stats_capture_hetsites.wex_protocol_investigation.txt"
#callstats.file = "/xchip/cga_home/mara/projects/ctc_sigma/analysis_raw_bam_files/call_stats_capture_hetsites.wex_wbc.txt"
#pair_type.file = "/xchip/cga_home/mara/projects/ctc_sigma/analysis_raw_bam_files/pair_type_ctc.wex_protocol_investigation.txt"

#callstats = read.table(callstats.file, sep="\t", header=T, quote='', stringsAsFactors=F)
#shortletter = read.table(pair_type.file, sep="\t", header=T, quote='', stringsAsFactors=F)


#######################
#Load in call stats files:
#######################
callstats = read.table(callstats.file, sep="\t", col.names=c("pair_id", "call_stats_capture_hetsites"), quote='', stringsAsFactors=F)
shortletter = read.table(shortletter, sep="\t", col.names=c("sample_id", "pair_type_ctc"), quote='', stringsAsFactors=F)

cs = data.frame(stringsAsFactors=F)
for (i in seq(1, dim(callstats)[1])) {
		filecs = as.character(callstats[i,"call_stats_capture_hetsites"])
		#pairid = as.character(shortletter[which(callstats$pair_id[i] %in% shortletter[,"pair_id"]),"pair_type_ctc"]) #need to map sample to name - how does this firehose output look?
		if (length(filecs) > 0) {
		if (file.exists(filecs)) {
			print(i)
			header = read.table(header=FALSE, nrow=1,file=as.character(filecs), sep="\t")
			take <- c("contig", "position","tumor_name","t_alt_count", "t_ref_count","n_alt_count","n_ref_count","failure_reasons","judgement", "dbsnp_site")
			takecols <- ifelse(t(header) %in% take, NA, 'NULL')
			callstats_all <- read.table(header = TRUE, colClasses=takecols, file=filecs, sep="\t")	
			tumor_name = as.character(unique(callstats_all$tumor_name))		
			callstats_all$type <- shortletter[shortletter$sample_id %in% tumor_name,"pair_type_ctc"]
			cs <- rbind(cs, callstats_all)
		} else {
			print(paste("call stats het file does not exist for:", callstats$pair_id))
		}
	}
}
#create general statistics to use later. 
cs$site = sprintf("%s:%s", cs$contig, cs$position)
cs$total_cov <- cs$t_alt_count + cs$t_ref_count
cs$af = cs$t_alt_count / cs$total_cov
cs[is.na(cs$af),"af"] = 0
cs = subset(cs, n_alt_count > 1 & dbsnp_site %in% "DBSNP" & t_alt_count+t_ref_count >=40)

print("min_alt to call a site covered is 3")
min_alt = 3

cs$ado_type = "both alleles present"
cs$ado_type[cs$t_alt_count==0 & cs$t_ref_count>=min_alt] <- "lost alt allele"
cs$ado_type[cs$t_alt_count>=min_alt & cs$t_ref_count==0] <- "lost ref allele"
cs$ado_type[cs$total_cov < min_alt] <- "lost both alleles"
#cs[is.na(cs$ado_type), "ado_type"] = "both alleles present"

pdf("allelic_distortion_plots.pdf")

#plot of allelic distortion for all samples: (density and histogram):
#ggplot(cs, aes(af, fill=tumor_name)) + geom_histogram()
print(ggplot(cs, aes(af, group=tumor_name, colour=type)) + geom_line(stat="density"))

#plot of fraction of het sites covered by indpendent libraries:
cs$tumor_name = as.factor(cs$tumor_name)
df <- ddply(cs, .(tumor_name), summarize, fraction = length(ado_typo[ado_type == "both alleles present"]) / length(ado_type) )
x = factor(cs$tumor_name, levels(cs$tumor_name)[order(df$fraction)])
print(ggplot(cs, aes(x, fill=ado_type)) + geom_bar(position="fill") + scale_y_continuous(limits = c(0,1)) +  ylab("Fraction Het Sites") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))

dev.off()

#make lorentz curves into a separate module. 

#######################
#TESTING

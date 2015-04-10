########################
#CTC - Summarize the individuals in FH
#Mara Rosenberg
#April 10, 2015
############################

firehose_sample.file = "/xchip/cga_home/mara/projects/ctc_sigma/CTC_Samples_10April2015.txt"

firehose_sample = read.table(firehose_sample.file, sep="\t", header=T, quote='', stringsAsFactors=F)

a = table(firehose_sample[,c("individual_id","cell_type")])

write.table(a, file="Summary_Data_Table.txt", sep="\t", row.names=T, col.names=T, quote=FALSE)
#############
#CTC: Upload bam files 

#argument 1: the full file path for the downloaded bam files
#argument 2: the template (from LIMS) for filling in with full bam file paths. 
#argument 3: sample set for upload to firehose

ls $1/*bam | grep -v Solexa > all_bams.txt

Rscript ../scripts/firehose_template_annotate.R all_bams.txt $2 $3

echo "importing samples to workspace: An_CTC_SIGMA_QC..."
fiss sample_import An_CTC_SIGMA_QC upload_lowpass_bams.txt

echo "importing sample set to workspace: An_CTC_SIGMA_QC..."
echo "sample set is $3"
fiss sset_import An_CTC_SIGMA_QC upload_sset.txt

echo "starting QC workflow..."
fiss flow_start An_CTC_SIGMA_QC $3 QC_AutoCorrelation_Workflow

ctc_firehose
============

Directions and upload scripts for connecting LIMS to Firehose (CTC QC and analysis)

1. Get Firehose upload file from LIMS 
2. IF walk up seq: Download BAM files to long term storage (/cga/golub/ctc/raw_data/)
3. RUN: cd /ctc_firehose/data/
4. RUN: sh ../../scripts/firehose_template_annotate_wrapper.sh [path_to_raw_bams] [FirehoseUploadTable[ lowpass_[date]


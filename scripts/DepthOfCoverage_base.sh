#!/bin/bash

# GATK toolkit

job_prefix='DepthOfCoverage';
hg19_ref=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta 
ucsc_hg19_ref='/humgen/gsa-hpprojects/GATK/data/ucsc.hg19/ucsc.hg19.fasta';
RefSeq='/xchip/singtex/references/GATK_data/refGene_b37.sorted.txt'
targets_interval_list='/xchip/singtex/references/refGene_sorted_b37.interval_list'
dbsnp_reference=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dbsnp.vcf
run_GATK="java -Xmx4000m -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar";
chrom_intervals='/xchip/singtex/references/chrom_intervals.list'

input_bam=$1
output_dir=$2
output_prefix=$3;
queue=$4

io_source_rd=`/broad/tools/scripts/io_resource_for_file /xchip/singtex` ;
job_group='/chengz/GATK';
job_project='GATK';
job_name=${output_prefix}.DOC
job_log=${output_dir}/${output_prefix}.log
bsub -q $queue -R "rusage[mem=4, ${io_source_rd}=4]" -g ${job_group} -P ${job_project} -J ${job_name} -o ${job_log} "${run_GATK} -T DepthOfCoverage -R ${hg19_ref} --partitionType sample --minMappingQuality 5 --summaryCoverageThreshold 1 --summaryCoverageThreshold 2 --summaryCoverageThreshold 5 --summaryCoverageThreshold 10 --summaryCoverageThreshold 50 --out ${output_dir}/${output_prefix}.WGS.DepthOfCoverage ${input_bam} -L 1:1-249250621"

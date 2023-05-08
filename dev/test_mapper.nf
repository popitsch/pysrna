#!/usr/bin/env nextflow

//params.config_file = "${workflow.launchDir}/srna-pipeline.config.json"

log.info "====================================="
log.info "Config file 			: ${params.config_file}"
log.info "Input data files      : ${params.data}"
log.info "Input transcriptome   : ${params.transcriptome}"
log.info "====================================="
log.info "\n"

/*
 * Create a channel for input files: either (unaligned) BAM files or fastq.gz/fq.gz files
 */
preprocessed_fq = Channel
    .fromPath( params.data)
    .ifEmpty { exit 1, "Cannot find any BAM/FASTQs matching: ${params.data}" }
    .map { file -> tuple(file.simpleName, file) }

transcriptome1 = Channel
    .fromPath( params.transcriptome)

/*
 * Prefix mapping with Tailor.
 * BAM files are fixed so they can be viewed in IGV.
 * Extra params passed to tailor via 'params.mapping_param.extra_param', e.g., '-v' to allow mismatches
 */
process map_reads {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false
	module 'samtools/1.10-foss-2018b:fastqc/0.11.8-java-1.8:python/3.7.2-gcccore-8.2.0'
	publishDir "results/mapped_reads", mode: 'copy'
	input:
		set name, file(fqz) from preprocessed_fq
		path(p1) from transcriptome1
	output:
		set name, file("${name}.bam"), file("${name}.bam.bai") into mapped_bam, mapped_bam2, mapped_bam3
		set file("*_fastqc*") into mapped_bam_fqc
	script:
    """
		gunzip -c ${fqz} > ${name}.fq
		${params.cmd.tailor_cmd} map -i ${name}.fq ${params.mapping_param.extra_param} -l ${params.mapping_param.min_prefix_match} -p ${params.dataset_name}.index.tailor -o ${name}.tailor.sam
		samtools sort -o ${name}.tailor.bam ${name}.tailor.sam
    	samtools index ${name}.tailor.bam
		${params.cmd.main_cmd} fix_tailor_bam --bam ${name}.tailor.bam --outdir .
		mv ${name}.tailor_fixed.bam ${name}.bam
		mv ${name}.tailor_fixed.bam.bai ${name}.bam.bai
		fastqc ${name}.bam
    """
}
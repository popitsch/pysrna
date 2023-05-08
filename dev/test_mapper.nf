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
Channel
    .fromPath( params.data)
    .ifEmpty { exit 1, "Cannot find any BAM/FASTQs matching: ${params.data}" }
    .map { file -> tuple(file.simpleName, file) }
    .into{ preprocessed_fq1; preprocessed_fq2; preprocessed_fq3 }

Channel
    .fromPath( params.transcriptome).collect().into{ transcriptome1; transcriptome2; transcriptome3 }

/*
 * Prefix mapping with Tailor.
 * BAM files are fixed so they can be viewed in IGV.
 * Extra params passed to tailor via 'params.mapping_param.extra_param', e.g., '-v' to allow mismatches
 */
process map_reads_tailor {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false
	module 'samtools/1.10-foss-2018b:fastqc/0.11.8-java-1.8:python/3.7.2-gcccore-8.2.0'
	publishDir "results/mapped_reads_tailor", mode: 'copy'
	input:
		set name, file(fqz) from preprocessed_fq1
		path(p1) from transcriptome1
	output:
		set name, file("${name}.bam"), file("${name}.bam.bai") into mapped_bam1
		set file("*_fastqc*") into mapped_bam_fqc1
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


process map_reads_srnaMapper {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false
	module 'bwa/0.7.17-foss-2018b'
	publishDir "results/mapped_reads_srnaMapper:fastqc/0.11.8-java-1.8", mode: 'copy'
	input:
		set name, file(fqz) from preprocessed_fq2
		path(p1) from transcriptome2
	output:
		set name, file("${name}.bam"), file("${name}.bam.bai") into mapped_bam2
		set file("*_fastqc*") into mapped_bam_fqc2
	script:
    """
		bwa index ${params.dataset_name}.fa
		gunzip -c ${fqz} > ${name}.fq
		${params.cmd.srnaMapper_cmd} -r ${name}.fq -g ${params.dataset_name} -o ${name}.sam -e 3
		samtools sort -o ${name}.bam ${name}.sam
    	samtools index ${name}.bam
		fastqc ${name}.bam
    """
}
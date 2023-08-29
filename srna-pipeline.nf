#!/usr/bin/env nextflow
log.info "====================================="
log.info "Config file 			: ${params.config_file}"
log.info "Dataset 				: ${params.dataset_name}"
log.info "Input data files      : ${params.data}"
log.info "====================================="
log.info "\n"

/*
 * Create a channel for input files: either (unaligned) BAM files or fastq.gz/fq.gz files 
 */
data_input = Channel
    .fromPath( params.data)
    .ifEmpty { exit 1, "Cannot find any BAM/FASTQs matching: ${params.data}" }
    .map { file -> tuple(file.simpleName, file) }

/*
 * Preprocessing: parse reads (align anchor seq, extract umis and sRBC sequences) and filter
 * if anchor not found. 
 */
process parse_reads {
	cpus 1
	memory '64 GB'
	time 4.h
	//cache false 
	publishDir "results/parsed_reads", mode: 'copy'
	input: 
		set name, file(dat) from data_input
	output: 
		set name, file("${name}.pass.fq.gz") into parsed_fq
		set name, file("${name}.filtered.fq.gz") into parsed_fq_filtered
		set name, file("${name}.stats.tsv.gz") into parsed_fq_stats
		set name, file("${name}.srbc_stats.tsv.gz") into parsed_fq_stats_srbc
		set file("*_fastqc*") into parsed_fqc
	script:
    """
		${params.cmd.main_cmd} parse_reads --config ${params.config_file} --config_prefix demux_param --dat ${dat} --out .
		gzip ${name}.pass.fq
		gzip ${name}.filtered.fq
		gzip ${name}.stats.tsv
		gzip ${name}.srbc_stats.tsv
		fastqc ${name}.fq.gz
		fastqc ${name}.filtered.fq.gz
    """
} 	

/**
* Remove remaining adapters and hard trim bases from 5'-end
*/
process preprocess_reads {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	publishDir "results/preprocessed_reads", mode: 'copy'
	input: 
		set name, file(fqz) from parsed_fq
	output: 
		set name, file("${name}_trimmed.fq.gz") into preprocessed_fq1
	script:
   """
		fastp -i ${fqz} -o ${name}_trimmed.fq --length_required ${params.demux_param.min_read_len} --adapter_sequence=${params.demux_param.anchor_seq} --trim_front1=${params.demux_param.fptrim}
		gzip ${name}_trimmed.fq
   """
} 	


/**
 * Count and filter reads stemming from spikein sequences
 */
process count_spikein_reads {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	publishDir "results/spikein_counts", mode: 'copy'
	input: 
		set name, file(fqz) from preprocessed_fq1
	output: 
		set name, file("${name}_trimmed.fq.pass.fq.gz") into preprocessed_fq, preprocessed_fq2
		set file("${name}_trimmed.fq.spikein.fq.gz"), file("${name}_trimmed.fq.counts_spikein.tsv") into preprocessed_fq_spikeins
	script:
   """
		${params.cmd.main_cmd} count_spikein_reads --config ${params.config_file} --config_prefix spikein_param --fq ${fqz} --out .
		gzip ${name}_trimmed.fq.spikein.fq
		gzip ${name}_trimmed.fq.pass.fq
   """
} 	

/**
 * Build transcriptome
 */
process build_transcriptome {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	publishDir "results/transcriptome", mode: 'copy'
	output: 
		tuple path("${params.dataset_name}.*") into transcriptome1, transcriptome2, transcriptome3, transcriptome4
	script:
   """	
		# create <name>.fa, <name>.gtf.gz and other transcriptome-related files
		${params.cmd.main_cmd} build_transcriptome --config ${params.config_file} --config_prefix "transcriptome_param" --out . --name ${params.dataset_name}
		# index with tailor
		${params.cmd.tailor_cmd} build -i ${params.dataset_name}.fa -p ${params.dataset_name}.index.tailor
   """
}

/**
 * Optional:Calc mappability of transcriptome
 * Only if transcriptome_param.calc_mappability=true
 */
process calc_transcriptome_mappability {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	publishDir "results/transcriptome", mode: 'copy'
	when:
		params.transcriptome_param.calc_mappability
	input:
    	path(p1) from transcriptome2
	output: 
		tuple path("${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz*") into mappability
	script:
   """	
		# calc mappability with genmap
		genmap index -F ${params.dataset_name}.fa -I ${params.dataset_name}.genmap_index
		genmap map -K ${params.transcriptome_param.mappability_k} -E ${params.transcriptome_param.mappability_e} -I ${params.dataset_name}.genmap_index -O ${params.dataset_name}.genmap_mappability -bg
		bedtools sort -i ${params.dataset_name}.genmap_mappability.bedgraph | bgzip > ${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz && tabix -p bed ${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz 
   """
}


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

/**
 * Optional:downsample to ensure IGV is not crashing. for QC only
 * Only if mapping_param.downsample_reads=true
 */
process downsample_reads {
	tag "${name}"
	cpus 1
	time 1.h
	publishDir "results/mapped_reads_downsampled", mode: 'copy'
	when:
		params.mapping_param.downsample_reads
	input:
    	set name, file(bam), file(bai) from mapped_bam2
    	path(p1) from transcriptome4
	output:
		file "*" into downsampled_reads
	script:
	"""
		${params.cmd.main_cmd} downsample_per_chrom --bam ${bam} -m 10000 --out .
	"""
}


/*
 * Optional: extract a sample of unmapped reads
 * Only if mapping_param.extract_unmapped_sample=true
 */
process extract_unmapped_reads {
	tag "${name}"
	cpus 1
	time 1.h
	publishDir "results/unmapped_reads_downsampled", mode: 'copy'
	when:
		params.mapping_param.extract_unmapped_sample
	input:
    	set name, file(bam), file(bai), file(fqz) from mapped_bam3.join(preprocessed_fq2)
	output:
		tuple path("${name}_unmapped_sample.fq.gz"), path("${name}_unmapped_sample.stats.txt") into unmapped_read_samples
	script:
	"""
		# print out mapped headers
		samtools view -F 4 ${bam} | awk '{ print "@"\$1 }' | sort -u -k 1b,1 > mapped.txt
		# grap only the fq headers.
		zcat ${fqz} | awk '(NR % 4 == 1)' |  sort -u -k 1b,1 > reads.txt
		# join and print only un-paired (-v)
		join --nocheck-order -v 1 reads.txt mapped.txt > unmapped.txt
		# print stats
		wc -l mapped.txt unmapped.txt > ${name}_unmapped_sample.stats.txt
		# create fastq
		head -n 1000 unmapped.txt > unmapped.head.txt
		LC_ALL=C zgrep --no-group-separator -x -A 3 -f unmapped.head.txt ${fqz} | gzip > ${name}_unmapped_sample.fq.gz
	"""
}

/*
 * Count reads
 */
process count_reads {
	tag "${name}"
	cpus 1
	time 1.h
	publishDir "results/counts", mode: 'copy'
	input:
    	set name, file(bam), file(bai) from mapped_bam
    	path(p1) from transcriptome3
	output:
		file "*" into read_counts
	script:
	"""
		${params.cmd.main_cmd} count_srna_reads --bam ${bam} --anno ${params.dataset_name}.gff3.gz --config ${params.config_file} --config_prefix "counting_param" --out . --name ${name}
	"""
}



/*
 * Optional: Create qc results
 * Only if calc_qc=true
 */
process qc_results {
    cpus 1
    publishDir "results/", mode: 'copy'
	when:
		params.calc_qc
    input:
    	file ("*") from read_counts.collect()
    output:
    	file("data.rds") into results
    	file("qc_plots/*pdf") into qc_plots
    script:
    """
        ${params.cmd.qc_cmd} ${params.config_file} .
        mkdir qc_plots && mv *.pdf qc_plots
    """
}

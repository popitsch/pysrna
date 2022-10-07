#!/usr/bin/env nextflow

//params.config_file = "${workflow.launchDir}/srna-pipeline.config.json"

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
    .ifEmpty { exit 1, "Cannot find any BAMs matching: ${params.bams}" }
    .map { file -> tuple(file.baseName, file) }

/*
 * Preprocessing: parse reads (align anchor seq, extract umis and sRBC sequences) and filter
 * if anchor not found. 
 */
process parse_reads {
	cpus 1
	memory '64 GB'
	time 4.h
	//cache false 
	module 'python/3.7.2-gcccore-8.2.0:fastqc/0.11.8-java-1.8'
	publishDir "results/parsed_reads", mode: 'copy'
	input: 
		set name, file(dat) from data_input
	output: 
		set name, file("${name}.fq.gz") into parsed_fq
		set name, file("${name}.filtered.fq.gz") into parsed_fq_filtered
		set name, file("${name}.stats.tsv.gz") into parsed_fq_stats
		set file("*_fastqc*") into parsed_fqc
	script:
    """
		${params.cmd.main_cmd} parse_reads --config ${params.config_file} --config_prefix demux_param --dat ${dat} --out .
		gzip ${name}.fq
		gzip ${name}.filtered.fq
		gzip ${name}.stats.tsv
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
	module 'fastp/0.20.1-gcc-8.2.0-2.31.1'
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
	module 'python/3.7.2-gcccore-8.2.0'
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
	module 'python/3.7.2-gcccore-8.2.0:bedtools/2.27.1-foss-2018b'
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
 * Calc mappability of transcriptome
 */
process calc_transcriptome_mappability {
	cpus 1
	memory '64 GB'
	time 2.h
	//cache false 
	module 'htslib/1.9-foss-2018b:bedtools/2.27.1-foss-2018b'
	publishDir "results/transcriptome", mode: 'copy'
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

/**
 * Downsample to ensure IGV is not crashing:: for QC only
 */
process downsample_reads {
	tag "${name}"
	cpus 1
	time 1.h
	module 'python/3.7.2-gcccore-8.2.0'
	publishDir "results/mapped_reads_downsampled", mode: 'copy'
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
 */
process extract_unmapped_reads {
	tag "${name}"
	cpus 1
	time 1.h
	module 'samtools/1.10-foss-2018b'
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
		samtools view -F 4 ${bam} | awk '{ print "@"\$1 }' | sort -u > mapped.txt
		# grap only the fq headers.
		zcat ${fqz} | awk '(NR % 4 == 1)' |  sort -u > reads.txt
		# join and print only un-paired (-v)
		join -v 1 reads.txt mapped.txt > unmapped.txt
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
	module 'python/3.7.2-gcccore-8.2.0'
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
 * Count with feature_counts
 */
/**process run_featureCounts {
	tag "${params.dataset_name}"
	cpus 1
	time 1.h
	module 'subread/2.0.1-gcc-7.3.0-2.30'
	publishDir "results/feature_counts", mode: 'copy'
	input:
    	set name, file(bam), file(bai) from mapped_bam
    	path(p1) from transcriptome3
	output:
		file "*" into feature_counts
	script:
		// -f                  Perform read counting at feature level
		// -O                  Assign reads to all their overlapping meta-features
		// -t <string>         Specify feature type(s) in a GTF annotation. 'exon' by default.
		// -g <string>         Specify attribute type in GTF annotation. 'gene_id' by default.
		// -M                  Multi-mapping reads will also be counted. 
		// -d <int>            Minimum fragment/template length, 50 by default.
		// -s <int or string>  Perform strand-specific read counting.
		// --fraction          Assign fractional counts to features. 
	"""
	featureCounts -T ${task.cpus} \
		-a ${params.dataset_name}.gtf.gz \
		-t ${params.counting_param.features} \
		-g ${params.counting_param.gene_id} \
		-o ${name}.featureCounts.txt \
		-O \
		-d ${params.mapping_param.min_prefix_match} \
		--primary \
		-f \
		-s 1 \
		--extraAttributes ${params.counting_param.extra_attributes} \
		${bam}
	"""
}
*/

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
import groovy.json.JsonOutput

if( !params.sample_sheet ) {
    error "Missing required parameter: --sample_sheet"
}

if( !file(params.sample_sheet).exists() ) {
    error "Sample sheet not found: ${params.sample_sheet}"
}

sample_sheet = file(params.sample_sheet)

log.info "====================================="
log.info "Config file           : ${params.config_file}"
log.info "Dataset               : ${params.dataset_name}"
log.info "Sample_sheet          : ${params.sample_sheet}"
log.info "Input data files      : ${params.data}"
log.info "=====================================\n"

////////////////////////////////////////////////////////////
// Input channel
////////////////////////////////////////////////////////////

data_input = Channel
    .fromPath(params.data)
    .ifEmpty { error "Cannot find any BAM/FASTQs matching: ${params.data}" }
    .map { f -> tuple(f.simpleName, f) }

////////////////////////////////////////////////////////////
// parse_reads
////////////////////////////////////////////////////////////

process PARSE_READS {
    cpus 1
    memory '64 GB'
    time 4.h
    publishDir "results/parsed_reads", mode: 'copy'

    input:
    tuple val(name), path(dat)
    path sample_sheet

    output:
    tuple val(name), path("${name}.pass.fq.gz"), emit: pass
    tuple val(name), path("${name}.filtered.fq.gz"), emit: filtered
    tuple path("${name}.stats.tsv.gz"), path("${name}.srbc_stats.tsv.gz"), emit: stats
    path "*_fastqc*", emit: qc

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    ${params.cmd.main_cmd} parse_reads \
        --config ${params.config_file} \
        --config_prefix demux_param \
        --dat ${dat} --out .

    gzip ${name}.pass.fq
    gzip ${name}.filtered.fq
    gzip ${name}.stats.tsv
    gzip ${name}.srbc_stats.tsv

    fastqc ${name}.pass.fq.gz
    fastqc ${name}.filtered.fq.gz
    """
}

////////////////////////////////////////////////////////////
// preprocess_reads
////////////////////////////////////////////////////////////

process PREPROCESS_READS {
    cpus 1
    memory '64 GB'
    time 2.h
    publishDir "results/preprocessed_reads", mode: 'copy'

    input:
    tuple val(name), path(fqz)

    output:
    tuple val(name), path("${name}_trimmed.fq.gz"), emit: trimmed

    script:
    """
    fastp \
        -i ${fqz} \
        -o ${name}_trimmed.fq \
        --length_required ${params.demux_param.min_read_len} \
        --adapter_sequence ${params.demux_param.anchor_seq} \
        --trim_front1 ${params.demux_param.fptrim}

    gzip ${name}_trimmed.fq
    """
}

////////////////////////////////////////////////////////////
// count_spikein_reads
////////////////////////////////////////////////////////////

process COUNT_SPIKEIN_READS {
    cpus 1
    memory '64 GB'
    time 2.h
    publishDir "results/spikein_counts", mode: 'copy'

    input:
    tuple val(name), path(fqz)

    output:
    tuple val(name), path("${name}_trimmed.fq.pass.fq.gz"), emit: pass
    tuple path("${name}_trimmed.fq.spikein.fq.gz"),
                     path("${name}_trimmed.fq.counts_spikein.tsv"),
                     emit: spikeins

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    ${params.cmd.main_cmd} count_spikein_reads \
        --config ${params.config_file} \
        --config_prefix spikein_param \
        --fq ${fqz} --out .

    gzip ${name}_trimmed.fq.pass.fq
    gzip ${name}_trimmed.fq.spikein.fq
    """
}

////////////////////////////////////////////////////////////
// build_transcriptome
////////////////////////////////////////////////////////////

process BUILD_TRANSCRIPTOME {
    cpus 1
    memory '64 GB'
    time 2.h
    publishDir "results/transcriptome", mode: 'copy'

    output:
    path("${params.dataset_name}.*")

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    ${params.cmd.main_cmd} build_transcriptome \
        --config ${params.config_file} \
        --config_prefix transcriptome_param \
        --out . --name ${params.dataset_name}

    ${params.cmd.tailor_cmd} build \
        -i ${params.dataset_name}.fa \
        -p ${params.dataset_name}.index.tailor
    """
}

////////////////////////////////////////////////////////////
// calc_transcriptome_mappability (optional)
////////////////////////////////////////////////////////////

process CALC_TRANSCRIPTOME_MAPPABILITY {
    cpus 1
    memory '64 GB'
    time 2.h
    publishDir "results/transcriptome", mode: 'copy'

    when:
    params.transcriptome_param.calc_mappability

    input:
    path(files)

    output:
    path "${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz*", emit: mappability

    script:
    """
    genmap index -F ${params.dataset_name}.fa -I ${params.dataset_name}.genmap_index
    genmap map \
        -K ${params.transcriptome_param.mappability_k} \
        -E ${params.transcriptome_param.mappability_e} \
        -I ${params.dataset_name}.genmap_index \
        -O ${params.dataset_name}.genmap_mappability -bg

    bedtools sort -i ${params.dataset_name}.genmap_mappability.bedgraph |
        bgzip > ${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz
    tabix -p bed ${params.dataset_name}.genmap_mappability.sorted.bedgraph.gz
    """
}

////////////////////////////////////////////////////////////
// map_reads
////////////////////////////////////////////////////////////

process MAP_READS {
    cpus 1
    memory '64 GB'
    time 2.h
    publishDir "results/mapped_reads", mode: 'copy'

    input:
    tuple val(name), path(fqz)
    path(txfiles)

    output:
    tuple val(name), path("${name}.bam"), path("${name}.bam.bai"), emit: bam
    path "*_fastqc*", emit: qc

    script:
    """
    gunzip -c ${fqz} > ${name}.fq

    ${params.cmd.tailor_cmd} map \
        -i ${name}.fq \
        ${params.mapping_param.extra_param} \
        -l ${params.mapping_param.min_prefix_match} \
        -p ${params.dataset_name}.index.tailor \
        -o ${name}.tailor.sam

    samtools sort -o ${name}.tailor.bam ${name}.tailor.sam
    samtools index ${name}.tailor.bam

    ${params.cmd.main_cmd} fix_tailor_bam \
        --bam ${name}.tailor.bam --outdir .

    mv ${name}.tailor_fixed.bam ${name}.bam
    mv ${name}.tailor_fixed.bam.bai ${name}.bam.bai

    fastqc ${name}.bam
    """
}

////////////////////////////////////////////////////////////
// downsample_reads (optional)
////////////////////////////////////////////////////////////

process DOWNSAMPLE_READS {
    cpus 1
    time 1.h
    publishDir "results/mapped_reads_downsampled", mode: 'copy'

    when:
    params.mapping_param.downsample_reads

    input:
    tuple val(name), path(bam), path(bai)
    path(txfiles)

    output:
    path "*", emit: downsampled

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    ${params.cmd.main_cmd} downsample_per_chrom \
        --bam ${bam} -m 10000 --out .
    """
}

////////////////////////////////////////////////////////////
// extract_unmapped_reads (optional)
////////////////////////////////////////////////////////////

process EXTRACT_UNMAPPED_READS {
    cpus 1
    time 1.h
    publishDir "results/unmapped_reads_downsampled", mode: 'copy'

    when:
    params.mapping_param.extract_unmapped_sample

    input:
    tuple val(name), path(bam), path(bai), path(fqz)

    output:
    tuple path("${name}_unmapped_sample.fq.gz"),
          path("${name}_unmapped_sample.stats.txt"),
          emit: unmapped

    script:
    """
    samtools view -F 4 ${bam} | awk '{print "@"\$1}' | sort -u > mapped.txt
    zcat ${fqz} | awk '(NR%4==1)' | sort -u > reads.txt
    join -v1 reads.txt mapped.txt > unmapped.txt

    wc -l mapped.txt unmapped.txt > ${name}_unmapped_sample.stats.txt

    head -n 1000 unmapped.txt > unmapped.head.txt
    zgrep -A3 -x -f unmapped.head.txt ${fqz} |
        gzip > ${name}_unmapped_sample.fq.gz
    """
}

////////////////////////////////////////////////////////////
// count_reads
////////////////////////////////////////////////////////////

process COUNT_READS {
    cpus 1
    time 1.h
    publishDir "results/counts", mode: 'copy'

    input:
    tuple val(name), path(bam), path(bai)
    path(txfiles)

    output:
    path "*", emit: counts

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    ${params.cmd.main_cmd} count_srna_reads \
        --bam ${bam} \
        --anno ${params.dataset_name}.gff3.gz \
        --config ${params.config_file} \
        --config_prefix counting_param \
        --out . --name ${name}
    rm config.json
    """
}

////////////////////////////////////////////////////////////
// qc_results (optional)
////////////////////////////////////////////////////////////

process QC_RESULTS {
    cpus 1
    publishDir "results", mode: 'copy'

    when:
    params.calc_qc

    input:
    val parsed_files
    val count_files
    val spikein_files
    path sample_sheet

    output:
    path "data.rds", emit: data
    path "qc_plots/*pdf", emit: plots

    script:
    """
    echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > config.json
    
    mkdir -p parsed_reads counts spikein_counts
    # Stage PARSE_READS outputs
    cp ${parsed_files.join(' ')} parsed_reads/
    cp ${spikein_files.join(' ')} spikein_counts/
    cp ${count_files.join(' ')} counts/
    
    ${params.cmd.qc_cmd} ${params.config_file} . .
    mkdir -p qc_plots
    mv *.pdf qc_plots
    """
}

////////////////////////////////////////////////////////////
// Workflow
////////////////////////////////////////////////////////////

workflow {

    parsed = PARSE_READS(data_input, sample_sheet)

    trimmed = PREPROCESS_READS(parsed.pass)

    spikein = COUNT_SPIKEIN_READS(trimmed.trimmed)

    tx = BUILD_TRANSCRIPTOME()

    CALC_TRANSCRIPTOME_MAPPABILITY(tx)

    mapped = MAP_READS(spikein.pass, tx)

    DOWNSAMPLE_READS(mapped.bam, tx)

    EXTRACT_UNMAPPED_READS(mapped.bam.join(spikein.pass))

    counted = COUNT_READS(mapped.bam, tx)

    QC_RESULTS(parsed.stats.collect(), counted.counts.collect(), spikein.spikeins.collect(), sample_sheet)
}


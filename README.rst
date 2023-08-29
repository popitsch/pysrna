# pysrna


`pysrna` is a python based analysis pipeline for small RNA-seq data.
It consists of a configurable [Nextflow](https://www.nextflow.io/) pipeline that orchestrates various
bioinformatics third-party tools as well as custom python scripts maintained in this project.

Contents
========

 * [Features](#features)
 * [Installation](#installation)
 * [Usage](#usage)
 * [Configuration](#configuration)
 * [Output Structure](#output-structure)

Features
========

TODO

+ A
+ B
+C

Installation
============

TODO

Third party tools
============

The following 3rd party tools are required to run the analysis pipeline. The shown version numbers \
correspond to our development/testing environment.

* nextflow/21.04.1
* python/3.7.2
* fastp/0.20.1
* fastqc/0.11.8
* bedtools/2.27.1
* samtools/1.10
* htslib/1.9
* tailor/1.1 (apptainer provided)
* R/4.0.2 (only for qc)

Execute the following command to check whether tailor runs as expected:
`apptainer run third_party_tools/tailor_1.1_linux_static.img tailor_v1.1_linux_static`

Configuration file
============

The following is a commented configuration file used for one of our experiments.

.. code-block:: JSON
    {
        "dataset_name": "Hen1_mouse_titr",  # Name of the dataset (used for file naming)
        "config_file": "config.json",       # This configuration file
        "sample_sheet": "sample_sheet.tsv", # sample sheet (see below)
        "data": "01_ngs_raw_mouse/*.fastq.gz", # glob pattern linking the input FASTQ files
        "cmd": {
                "main_cmd": "python srna-pipelines/python/main.py", # command for executing the main python script
                "tailor_cmd": "singularity run third_party_tools/tailor_1.1_linux_static.img tailor_v1.1_linux_static",
                "qc_cmd": "Rscript --vanilla srna-pipelines/R/srna_qc.R" # command for executing the R qc script
        },
        "demux_param": {
                "anchor_seq": "AGATCGGAAGAGCACACGTCT",  # expected adapter sequence
                "umi_len": 6,                           # UMI length
                "srbc_len": 5,                          # sRBC length
                "min_aln_score": 0.9,                   # minimum (length-normalized) alignment score
                "min_read_len": 18,                     # minimum accepted read length
                "fptrim": 4,                            # Number of 5'-end hard-trimmed nucleotides
                "filter_wrong_srbc": true               # If true, reads with unexpected sRBC will be filtered
        },
        "spikein_param": {
                "spikein_fa": "ref/spikeins.fa",        # FASTA file containing the spike-in sequences
                "spikein_meta":  "ref/spikeins_meta.tsv"# TSV file containing spike-in metadata (see below)
        },
        "transcriptome_param": {
                "genome_fa": "ref/genomes/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa", # reference genome
                "gene_anno": "ref/mirgenedb/mmu_nochr.sorted.gff3.gz", # mirgenedb annotation file
                "main_feature": "pre_miRNA", # ID of the main feature to be considered
                "gene_id": "ID",    # Name of the GFF attribute used to store gene ids
                "gene_name": "ID",  # Name of the GFF attribute used to store gene names
                "mappability_k": 18, # k-parameter for calculating the mappability with genmap
                "mappability_e": 2,  # e-parameter for calculating the mappability with genmap
                "amp_extension": 25  # Number of bases added before/after the created transcriptome sequences
        },
        "mapping_param": {
            "min_prefix_match": 18, # prefix length as passed to tailor via the -l parameter
            "extra_param": ""      # possible extra parameters passd to tailor
        },
        "counting_param": {
                "features": {
                        "pre_miRNA": "general_profile", # mapping of 'pre_miRNA' GFF entries to the general counting profile
                        "miRNA": "miRNA_profile"        # mapping of 'miRNA' GFF entries to the mature miRNA counting profile
                        },
                "gene_id": "ID",                        # Name of the GFF attribute used to store gene ids
                "extra_attributes": "Name",             # Extra GFF attributes that will be copied to the output file
                "write_bam": true                       # If true, then debugging BAM files will be written.
        },
        "calc_qc": false    # If true, QC R script will be called.
    }


Sample sheet
============

The sample_sheet is a simple TSV file that maps samples (identified by their FASTQ filename prefix) to
the expected sRBC sequence. For running the QC pipeline it additionaly requires a column containing the number
of raw (sequenced) reads per sample (column: raw_reads). If a genotype column is provided, then it will
be used by the QC Rmd to group/colour some of the analysis plots.

Additionally, the sample_sheet.tsv may contain arbitrary optional meta-data columns that are useful for
subsequent data analysis/QC. The following is an example sample sheet used for one of our experiments.

.. code-block::
    sample_num  filename_prefix         NaIO4_oxidized  perc_rna  sample_name   sRBC_adaptor  sRBC   organism  raw_reads
    1           224823_S18_L004_R1_001  yes             0.01      hen2_ox_001   1             CAGTG  mouse     14061279
    2           224824_S19_L004_R1_001  yes             0.1       hen2_ox_01    2             AGCAA  mouse     10534236
    3           224825_S20_L004_R1_001  yes             1         hen2_ox_1     3             GGTAT  mouse     9161046

- raw_reads: number of raw reads, used for plotting filtering statistics
- genotype: genotype of sample, used for plotting in srna_qc.R


Spikein meta data file
======================

This is a simple TSV file containing the following 3 columns:

* si_name : name of the spike-in as provided in the spike-in FASTA file
* si_len : length of the spike-in seqeunce
* si_conc : expected concentration of the spike-in.
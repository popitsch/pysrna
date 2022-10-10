#!/usr/bin/env Rscript
require(data.table)
require(rjson)
require(dplyr)
require(tidyr)
require(ggplot2)
require(stringr)
require(cowplot)
require(readr)
require(testthat)
require(tidylog)
require(glue)
require(rtracklayer)
require(tximport)
require(DESeq2)
require(tibble)
library(broom)
require(ggpmisc)

count=dplyr::count
select=dplyr::select

# to ensure stripping '\0' (nul) from character vector
options(arrow.skip_nul = TRUE)

# load data from TSV. Supports gzipped files
load_table = function(dataF, append="", header=T, nrows=Inf, skip=0, colClasses=NULL) {
  dataF = as.character(paste0(dataF, append))
  print(paste("Loading", dataF))
  if ( endsWith(dataF, ".gz") ) {
    return(fread(cmd=paste('gunzip -c', dataF), header=header, sep="\t", skip=skip, na.strings=c("na","NA",".", "None"), nrows=nrows, colClasses=colClasses))
  } else {
    return(fread(dataF, header=header, sep="\t", skip=skip, na.strings=c("na","NA",".", "None"), nrows=nrows,
                 colClasses=colClasses))
  }
}
# multiple plots with single title
my_plot_grid = function(..., main=NULL) {
  plot_row=plot_grid(...)
  title <- ggdraw() + draw_label( main, fontface = 'bold', x = 0, hjust = 0 ) + theme(plot.margin = margin(0, 0, 0, 7))
  return (plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1)))
}


# load data
build_dataset = function(sample_sheet, config, results_dir) {
  parsed_read_stats=list()
  counted_read_stats=list()
  srbc_stats=list()
  cnt=list()
  iso=list()
  tai=list()
  sc=list()
  for ( s in sample_sheet %>% pull(sample_name)) {
    fn=sample_sheet %>% filter(sample_name==s) %>% pull(filename_prefix)
    parsed_read_stats[[s]] = load_table(paste0(results_dir, '/parsed_reads/',fn, '.stats.tsv.gz'))
    srbc_stats[[s]] = load_table(paste0(results_dir, '/parsed_reads/',fn, '.srbc_stats.tsv.gz'))
    sc[[s]] = load_table(paste0(results_dir, '/spikein_counts/', fn, '_trimmed.fq.counts_spikein.tsv'))
    cnt[[s]] = load_table(paste0(results_dir, '/counts/', fn, '.counts.tsv')) %>% 
      left_join( load_table(paste0(results_dir, '/counts/', fn, '.meta.tsv')), by=c('gene_id', 'type') )
    iso[[s]] = load_table(paste0(results_dir, '/counts/', fn, '.iso.tsv')) %>% 
      left_join( load_table(paste0(results_dir, '/counts/', fn, '.meta.tsv')), by=c('gene_id', 'type') )
    tai[[s]] = load_table(paste0(results_dir, '/counts/', fn, '.tails.tsv'),
                          colClasses=c('gene_id'='character','type'='character','tail'='character')) %>%  # set colClasses for empty tables!
      left_join( load_table(paste0(results_dir, '/counts/', fn, '.meta.tsv')), by=c('gene_id', 'type') )
    counted_read_stats[[s]] = load_table(paste0(results_dir, '/counts/',fn, '.stats.tsv'))
  }
  parsed_read_stats = parsed_read_stats %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name)))
  counted_read_stats = counted_read_stats %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name)))
  srbc_stats = srbc_stats%>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name)))
  sc = sc %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name))) %>% 
    left_join(load_table(conf$spikein_param$spikein_meta), by='si_name') %>% 
    mutate(methylated=ifelse(grepl('mX',si_name, fixed = TRUE), 'yes', 'no' )) %>% 
    left_join(sample_sheet, by='sample_name') 
  sc = sc %>% 
    mutate(si_name=factor(si_name, levels=unique(sc$si_name)))
  
  d=list()
  d[['sample_sheet']] = sample_sheet
  d[['cnt']] = cnt %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name))) %>% 
    left_join(sample_sheet, by='sample_name') %>% 
    mutate(frac_tailed=reads_tailed/reads,
           frac_tailed_norm=reads_tailed_norm/reads_norm,
           frac_invalid=reads_invalid/(reads+reads_invalid))
  d[['iso']] = iso %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name))) %>% 
    left_join(sample_sheet, by='sample_name')
  d[['tai']] = tai %>% 
    bind_rows(.id='sample_name') %>% 
    mutate(sample_name=factor(sample_name, levels=levels(sample_sheet$sample_name))) %>% 
    left_join(sample_sheet, by='sample_name') %>% 
    mutate(tail_len=nchar(tail)) 
  d[['sc']] = sc
  d[['parsed_read_stats']]=parsed_read_stats
  d[['counted_read_stats']]=counted_read_stats
  d[['srbc_stats']]=srbc_stats
  return (d)
}

# ===============================================================
# MAIN
# ===============================================================
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1 | length(args)>3) {
  stop("usage: srna_qc.R <config> [<outdir>] [<resultsdir>]", call.=FALSE)
} 
srna_config=args[1]
home_dir=paste0(dirname(srna_config),'/')
if (length(args)==1) {
  resultsdir=paste0(home_dir,'/results/')
  outdir=paste0(resultsdir,'analysis/')
} else if (length(args)==2) {
  resultsdir=paste0(home_dir,'/results/')
  outdir=args[2]
} else {
  resultsdir=args[3]
  outdir=args[2]
}
# create results dir?
if (!dir.exists(resultsdir)) {
  print("Results dir not found")
  exit()
} 
# create results dir?
if (!dir.exists(outdir)) {
  dir.create(outdir)
} 

conf=fromJSON(paste(readLines(srna_config), collapse=""))

sample_sheet = load_table(conf$sample_sheet) %>% 
  mutate(sample_name=factor(sample_name, levels = unique(sample_name)))
d = build_dataset(sample_sheet, conf, resultsdir)

data_file=paste0(outdir,'/data.rds')
saveRDS(d, data_file, compress = FALSE)
print(paste0("Done. Load data via d=readRDS(data_file.rds)"))
# ===============================================================
# Filtering stats
# ===============================================================
p1 = d$parsed_read_stats %>% 
  group_by(sample_name) %>% 
  filter(category=='read_count') %>% 
  left_join(sample_sheet %>% select(sample_name, raw_reads), by='sample_name') %>% 
  mutate(frac=value/raw_reads) %>% 
  ggplot(aes(key, frac, fill=key)) +
  geom_col() +
  scale_fill_manual( values = c( "pass"="darkgreen", "filtered"="red", 
                                 "wrong_srbc"="grey", "too_short"="grey", "no_adapter"="grey"  ), guide = "none" ) +
  facet_wrap(sample_name~.) +
  ylim(0,1) + xlab("") + ylab("") +
  coord_flip() + 
  ggtitle("Filtering stats")

known_srbcs=sample_sheet %>% pull(sRBC)
p2 = d$parsed_read_stats %>% 
  filter(category=='srbc') %>% 
  left_join(sample_sheet %>% select(sample_name, sRBC, raw_reads), by='sample_name') %>% 
  group_by(sample_name) %>% 
  mutate(true_srbc=case_when(key==sRBC~'true',key %in% known_srbcs~'known', TRUE ~'other')) %>% 
  mutate(total=sum(value)) %>% 
  group_by(sample_name, true_srbc, total) %>% 
  summarise(share=sum(value)) %>% 
  ungroup() %>% 
  mutate(frac=share/total) %>% 
  ggplot(aes(true_srbc, frac, fill=true_srbc)) +
  geom_col() +
  scale_fill_manual( values = c( "true"="darkgreen", "other"="red" ), guide = "none" ) +
  facet_wrap(sample_name~.)+
  ylim(0,1) + xlab("") + ylab("") +
  coord_flip() + 
  ggtitle("Fraction of reads with true, known or unknown (other sequence) sRBC 5-mer")

p3 = d$counted_read_stats %>% 
  left_join(sample_sheet, by='sample_name') %>% 
  pivot_wider(names_from='key', values_from='value') %>% 
  mutate(frac=counted_reads/mapped_reads) %>% 
  ggplot(aes(sample_name,frac, fill=genotype)) +
  geom_col(position = 'dodge') +
  #scale_fill_manual( values = c( "NaIO4_oxidized"="darkgreen" ), guide = "none" ) +
  xlab("") + ylab("") +
  coord_flip() + 
  ggtitle("Fraction of mapped-reads assigned to pri_miRNA annotations")

p4 = d$counted_read_stats %>% 
  left_join(sample_sheet, by='sample_name') %>% 
  filter(key %in% c('counted_reads')) %>% 
  ggplot(aes(sample_name,value, fill=genotype)) +
  geom_col() + xlab("") + ylab("") +
  coord_flip() + 
  ggtitle("Raw numbers of pri_miRNA assigned reads")

my_plot_grid(p1,p2,p3,p4, labels=c('A','B','C','D'), main='Filter statistics')
ggsave(paste0(outdir,'/qc_pipeline_filtering_stats.pdf'), width=12, height=9)


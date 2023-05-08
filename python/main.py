'''
@author: niko.popitsch@univie.ac.at
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys, json, logging, random, time, gzip, math, random
from collections import Counter, OrderedDict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import pandas as pd
import numpy as np
from intervaltree import IntervalTree, Interval
import pysam
import tqdm
import hashlib
from pathlib import Path
from itertools import zip_longest
from utils import *

def build_transcriptome_section(anno_file, anno_fmt, genome_fa, main_feature, gene_name, amp_extension, padding, out_file_fa, out_fa, out_fa_chrsize, out_fa_dict, out_gff3):
    """ Builds transcriptome from the passed genome/gtf and writes a GFF file"""
    # create dict feature_chrom: (amplicon_start, amplicon_end, amplicon_name) from main features (to avoid long names)
    genome = pysam.FastaFile(genome_fa) # @UndefinedVariable
    ampcoords={}
    ampchrom_names=set()
    for gtf_line in tqdm.tqdm(pysam.TabixFile(anno_file, mode="r").fetch(parser=pysam.asTuple())): # @UndefinedVariable
        if gtf_line[2] == main_feature:
            fchrom,fstart,fend,info=gtf_line[0], int(gtf_line[3])-1, int(gtf_line[4]),parse_info(gtf_line[8], fmt=anno_fmt)
            if fchrom not in genome.references:
                   print("WARN: skipping feature as target chromosome %s not found in reference FASTA %s!" % (fchrom, genome_fa) )
                   continue         
            ampstart=max(1,fstart-amp_extension)
            ampend=fend+amp_extension
            # new chrom name
            ampchrom=info[gene_name].lower() if gene_name in info else fchrom+"_"+str(fstart+1)+"_"+str(fend)
            while ampchrom in ampchrom_names: # make sure that there are no name collisions!
                ampchrom=ampchrom+'_1'                
            if fchrom not in ampcoords:
                ampcoords[fchrom]=IntervalTree()
            ampcoords[fchrom].addi(ampstart,ampend, ampchrom)  
    # merge overlapping intervals and concat names
    for fchrom in ampcoords:
        ampcoords[fchrom].merge_overlaps(data_reducer=lambda a,b: '+'.join([a,b]), strict=False)
    # create FASTA + GFF 
    written_fa=set()  
    for gtf_line in tqdm.tqdm(pysam.TabixFile(anno_file, mode="r").fetch(parser=pysam.asTuple())): # @UndefinedVariable
        fchrom,fstart,fend,info=gtf_line[0], int(gtf_line[3])-1, int(gtf_line[4]),parse_info(gtf_line[8], fmt=anno_fmt)
        if fchrom not in ampcoords:
            continue # skip as gtf chrom not in genome
        gi=next(iter(ampcoords[fchrom][fstart:fend])) # get amplicon interval (genomic coords)
        ampstart, ampend, ampchrom=gi.begin, gi.end, gi.data
        offset=padding+(fstart-ampstart)+1
        seq='N'*padding+genome.fetch(fchrom, ampstart, ampend)+'N'*padding
        ampchrom_len=len(seq)
        # FASTA                    
        if ampchrom not in written_fa:        
            print('>%s'% ampchrom, file=out_fa)
            print(format_fasta(seq), file=out_fa)
            # chrom.sizes file
            print('%s\t%i'% (ampchrom,ampchrom_len), file=out_fa_chrsize)
            # DICT
            print('@SQ\tSN:%s\tLN:%i\tM5:%s\tUR:file:%s'% (ampchrom,
                                                           ampchrom_len,
                                                           hashlib.md5(seq.encode('utf-8')).hexdigest(), # M5 MD5 checksum of the sequence in the uppercase, excluding spaces but including pads (as ‘*’s)
                                                           os.path.abspath(out_file_fa)), file=out_fa_dict)
            written_fa.add(ampchrom)
        gtf_line=list(gtf_line)
        gtf_line[0]=ampchrom
        gtf_line[3]=offset
        gtf_line[4]=offset+fend-fstart-1
        gtf_line[8]=';'.join(["%s=%s"%(k,v) for k,v in info.items()]) # write in GFF3 fmt
        print("\t".join([str(x) for x in gtf_line]), file=out_gff3)      
    #
    # genome = pysam.FastaFile(genome_fa) # @UndefinedVariable
    # new_chrom_names=set()
    # gcoords={}
    # for gtf_line in tqdm.tqdm(pysam.TabixFile(anno_file, mode="r").fetch(parser=pysam.asTuple())): # @UndefinedVariable
    #     fchrom,fstart,fend=gtf_line[0], int(gtf_line[3])-1, int(gtf_line[4])
    #     ampstart=max(1,fstart-amp_extension)
    #     ampend=fend+amp_extension
    #     offset=padding+(fstart-ampstart)+1
    #     if gtf_line[2] == main_feature:  
    #         if fchrom not in genome.references:
    #            print("WARN: skipping feature as target chromosome %s not found in reference FASTA %s!" % (fchrom, genome_fa) )
    #            continue        
    #         seq='N'*padding+genome.fetch(fchrom, ampstart, ampend)+'N'*padding
    #         info=parse_info(gtf_line[8], fmt=anno_fmt) 
    #         new_chrom=info['Name'] if 'Name' in info else fchrom+"_"+str(fstart+1)+"_"+str(fend)
    #         while new_chrom in new_chrom_names: # make sure that there are no name collisions!
    #             new_chrom=new_chrom+'_1'
    #         new_chrom_names.add(new_chrom)
    #         new_chrom_len=ampend-ampstart+1+2*padding # length of new chrom
    #         # FASTA                            
    #         print('>%s'% new_chrom, file=out_fa)
    #         print(format_fasta(seq), file=out_fa)
    #         # chrom.sizes file
    #         print('%s\t%i'% (new_chrom,new_chrom_len), file=out_fa_chrsize)
    #         # DICT
    #         print('@SQ\tSN:%s\tLN:%i\tM5:%s\tUR:file:%s'% (new_chrom,
    #                                                        new_chrom_len,
    #                                                        hashlib.md5(seq.encode('utf-8')).hexdigest(), # M5 MD5 checksum of the sequence in the uppercase, excluding spaces but including pads (as ‘*’s)
    #                                                        os.path.abspath(out_file_fa)), file=out_fa_dict)
    #         gtf_line=list(gtf_line)
    #         gtf_line[0]=new_chrom
    #         gtf_line[3]=offset
    #         gtf_line[4]=offset+fend-fstart-1
    #         gtf_line[8]=';'.join(["%s=%s"%(k,v) for k,v in info.items()]) # write in GFF3 fmt
    #         print("\t".join([str(x) for x in gtf_line]), file=out_gff3)   
    #         if fchrom not in gcoords:
    #             gcoords[fchrom]=IntervalTree()
    #         gcoords[fchrom].addi(fstart,fend,(new_chrom, offset))  
    #     else:
    #         # add enveloped features that are no 'main features'
    #         if fchrom not in gcoords:
    #             continue
    #         reg=gcoords[fchrom][fstart:fend]
    #         if len(reg)!=1: # there must be exactly one candidate
    #             continue
    #         if not next(iter(reg)).contains_interval(Interval(fstart,fend)):
    #             continue
    #         new_chrom, offset=next(iter(reg)).data  
    #         offset+=fstart-next(iter(reg))[0]
    #         gtf_line=list(gtf_line)
    #         gtf_line[0]=new_chrom
    #         gtf_line[3]=offset
    #         gtf_line[4]=offset+fend-fstart-1
    #         print("\t".join([str(x) for x in gtf_line]), file=out_gff3)   
    #     gtf_line[2]
        
def build_transcriptome(name, config, outdir, config_prefix=[]):
    """ Build a transcriptome file that includes (optional) spike-in sequences """
    genome_fa = get_config(config, config_prefix+['genome_fa'], required=True)
    gene_anno_file=get_config(config, config_prefix+['gene_anno'], required=True)
    gene_anno_fmt='GTF' if '.gtf' in gene_anno_file else 'GFF3'
    spikein_fa=get_config(config, config_prefix+['spikein_fa'], required=False)
    padding=get_config(config, config_prefix+['padding'], default_value=10, required=False)
    main_feature=get_config(config, config_prefix+['main_feature'], default_value='exon', required=False)
    amp_extension=get_config(config, config_prefix+['amp_extension'], default_value=10, required=False) # extend amplicon seq updownstream
    gene_name=get_config(config, config_prefix+['gene_name'], default_value='Name', required=False) # default is miRbase default name
    gene_id=get_config(config, config_prefix+['gene_id'], default_value='ID', required=False) # default is miRbase default ID
    #
    out_file_fa=outdir+'/'+name+'.fa'
    out_file_fa_chrsize=outdir+'/'+name+'.fa.chrom.sizes'
    out_file_fa_dict=outdir+'/'+name+'.fa.dict'
    out_file_feature_meta=outdir+'/'+name+'.feature_meta.tsv'
    out_file_gff3=outdir+'/'+name+'.gff3'
    out_file_spikeins_anno=None
    with open(out_file_feature_meta, 'w') as out_feature_meta:
        with open(out_file_fa, 'w') as out_fa:
            with open(out_file_fa_chrsize, 'w') as out_fa_chrsize:
                with open(out_file_fa_dict, 'w') as out_fa_dict:
                    with open(out_file_gff3, 'w') as out_gff3:
                        # create meta TSV file
                        print('\t'.join(['gene_id', 'feature_type']), file=out_feature_meta)
                        for gtf_line in pysam.TabixFile(gene_anno_file, mode="r").fetch(parser=pysam.asTuple()):
                            info=parse_info(gtf_line[8], fmt=gene_anno_fmt) 
                            if gene_id in info:
                                print('\t'.join([info[gene_id], gtf_line[2]]), file=out_feature_meta)
                        # build transcriptome GTF file
                        build_transcriptome_section(gene_anno_file, gene_anno_fmt, genome_fa, main_feature, gene_name, amp_extension, padding, out_file_fa, out_fa, out_fa_chrsize, out_fa_dict, out_gff3)
                        # add tsv/gtf info for spikeins
                        if spikein_fa is not None:
                            out_file_spikeins_anno=outdir+'/'+Path(spikein_fa).stem+('.gtf' if gene_anno_fmt is 'GTF' else '.gff3')
                            spikein=pysam.FastaFile(spikein_fa) # @UndefinedVariable
                            with open(out_file_spikeins_anno, 'w') as out_spikeins_gtf:
                                for spikein_chr in spikein.references:
                                    # add spikein annotation. They will be 'main' features.
                                    if gene_anno_fmt is 'GTF':
                                        info_line=gene_id+' "'+spikein_chr+'"; '+gene_name+' "'+spikein_chr+'";'
                                    else:   
                                        info_line=gene_id+'='+spikein_chr+';'+gene_name+'='+spikein_chr+';'
                                    gtf_line=[spikein_chr, 'spikein', main_feature, 1, spikein.get_reference_length(spikein_chr),
                                              '.', '+', '.', info_line]
                                    print("\t".join([str(x) for x in gtf_line]), file=out_spikeins_gtf) 
                                    print('\t'.join([spikein_chr, 'spikein']), file=out_feature_meta)
                            sort_bgzip_and_tabix(out_file_spikeins_anno, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)
                            out_file_spikeins_anno=out_file_spikeins_anno+'.gz'
                            build_transcriptome_section(out_file_spikeins_anno, gene_anno_fmt, spikein_fa, main_feature, gene_name, 0, padding, out_file_fa, out_fa, out_fa_chrsize, out_fa_dict, out_gff3)
    # compress + index output files            
    sort_bgzip_and_tabix(out_file_gff3, seq_col=0, start_col=3, end_col=4, line_skip=0, zerobased=False)
    pysam.faidx(out_file_fa)# @UndefinedVariable    
    print("Created resources:")
    print("FASTA file + idx:\t"+out_file_fa)
    print("CHROMSIZE file:\t"+out_file_fa_chrsize)
    print("DICT file:\t"+out_file_fa_dict)
    print("GFF3 file:\t"+out_file_gff3)
    print("TSV file:\t"+out_file_feature_meta)
    
def parse_read(query_name, query_sequence, query_qualities, config, config_prefix=[], expected_srbc=None):
    """ Parse single reads """
    anchor_seq=get_config(config, config_prefix+['anchor_seq'], required=True)
    umi_len=get_config(config, config_prefix+['umi_len'], required=True)
    srbc_len=get_config(config, config_prefix+['srbc_len'], required=True)
    aln = pairwise2.align.globalxs( # globalxs(sequenceA, sequenceB, open, extend) -> alignments
            anchor_seq, 
            query_sequence, 
            -2, 
            -1,
            penalize_end_gaps=(False,False), # penalize starting/ending gaps
            score_only=False,
            one_alignment_only=True)[0]
    _, _, aln_score, _, _ = aln
    aln_score=aln_score/len(anchor_seq)
    if aln_score < get_config(config,config_prefix+['min_aln_score'], 0.9):
        return True, query_name, query_sequence, query_qualities, aln_score, None, None, 'no_adapter'
    # get start pos of anchor in read
    start=len(aln[0])-len(aln[0].lstrip('-'))
    start_srbc=start-srbc_len
    start_umi=start_srbc-umi_len
    #end=len([x for x in aln[1][:len(aln[0].rstrip('-'))] if x != '-'])
    if (start_srbc<0) or (start_umi<0):
        return True, query_name, query_sequence, query_qualities, aln_score, None, None, 'too_short' # filter if remaining read too short
    # extract umi, srbc and clip read
    srbc=query_sequence[start_srbc:start_srbc+srbc_len]
    umi=query_sequence[start_umi:start_umi+umi_len]
    trimmed_query_sequence=query_sequence[:start_umi]
    trimmed_query_qualities=query_qualities[:start_umi]
    if len(trimmed_query_sequence)<get_config(config,config_prefix+['min_read_len'], 22):
        return True, query_name, query_sequence, query_qualities, aln_score, umi, srbc, 'too_short' # filter if remaining read too short
    if (expected_srbc is not None) and (srbc!=expected_srbc):
        return True, query_name, query_sequence, query_qualities, aln_score, umi, srbc, 'wrong_srbc' # filter if srbc is not matching expectation
    trimmed_query_name='_'.join([query_name, srbc, umi])
    return False, trimmed_query_name, trimmed_query_sequence, trimmed_query_qualities, aln_score, umi, srbc, None

def parse_reads(dat_file, config, outdir, config_prefix=[]):
    """ Extract info from reads """
    sample_name=Path(dat_file[:-3]).stem if dat_file.endswith('.gz') else Path(dat_file).stem
    out_prefix=outdir+sample_name
    filter_wrong_srbc= get_config(config, config_prefix+['filter_wrong_srbc'], required=True)
    if filter_wrong_srbc:
        sample_sheet=pd.read_csv(get_config(config, ['sample_sheet'], required=True), sep='\t')
        srbcs={k:v.strip() for k,v in zip(sample_sheet['filename_prefix'],sample_sheet['sRBC']) }
        expected_srbc = srbcs[sample_name]
    else:
        expected_srbc=None
        srbcs=None
    stats=Counter()
    stats['umi', 'pass']=set()
    srbc_stats=Counter()
    mean_aln_score_filtered, mean_aln_score_pass=list(), list()
    with open(out_prefix+'.filtered.fq', 'w') as fq_filtered:
        with open(out_prefix+'.pass.fq', 'w') as fq_pass:
            if dat_file.endswith(".fastq.gz") or dat_file.endswith(".fq.gz"):
                with gzip.open(dat_file, 'r') as in1:
                    it1=grouper(in1, 4, '')
                    for read in tqdm.tqdm(it1):
                        query_name, query_sequence, _, query_qualities = [x.decode("utf-8").strip() for x in read]
                        query_name=query_name[1:] # remove '@' prefix
                        filtered, query_name, query_sequence, query_qualities, aln_score, umi, srbc, filter_str=parse_read(query_name, query_sequence, query_qualities, config, config_prefix=config_prefix, expected_srbc=expected_srbc)
                        if filtered:                    
                            out=fq_filtered 
                            stats['read_count', 'filtered']+=1    
                            stats['read_count', filter_str]+=1                  
                            mean_aln_score_filtered.append(aln_score)                
                        else:
                            out=fq_pass          
                            stats['read_count', 'pass']+=1           
                            stats['umi', 'pass'].add(umi)           
                            mean_aln_score_pass.append(aln_score)                
                        print('@%s\n%s\n+\n%s' % (query_name,
                                                  query_sequence, 
                                                  query_qualities if query_qualities is not None else '???'), 
                                                  file=out)
                        srbc_stats[
                            'NA' if srbc is None else srbc,
                            '1' if filtered else '0',
                            'NA' if expected_srbc is None else '1' if srbc==expected_srbc else '0',
                            'NA' if srbcs is None else '1' if srbc in srbcs.values() else '0'
                            ]+=1
            elif dat_file.endswith(".bam"):# unaligned BAM file inpiut
                samfile = pysam.AlignmentFile(dat_file, "rb", check_sq=False)# @UndefinedVariable
                for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
                    query_name, query_sequence, query_qualities = read.query_name, read.query_sequence, read.query_qualities
                    filtered, query_name, query_sequence, query_qualities, aln_score, umi, srbc, filter_str=parse_read(query_name, query_sequence, query_qualities, config, config_prefix=config_prefix, expected_srbc=expected_srbc)
                    if filtered:                    
                        out=fq_filtered 
                        stats['read_count', 'filtered']+=1    
                        stats['read_count', filter_str]+=1                  
                        mean_aln_score_filtered.append(aln_score)                
                    else:
                        out=fq_pass          
                        stats['read_count', 'pass']+=1           
                        stats['umi', 'pass'].add(umi)           
                        mean_aln_score_pass.append(aln_score)                
                    print('@%s\n%s\n+\n%s' % (query_name,
                                              query_sequence, 
                                              ''.join(map(lambda x: chr( x+33 ), query_qualities)) if query_qualities is not None else '???'), 
                                              file=out)
                    srbc_stats[
                            'NA' if srbc is None else srbc,
                            '1' if filtered else '0',
                            'NA' if expected_srbc is None else '1' if srbc==expected_srbc else '0',
                            'NA' if srbcs is None else '1' if srbc in srbcs.values() else '0'
                            ]+=1
                samfile.close()
            else:
                print("ERR: Unknown input file format for file %s" % dat_file)
                sys.exit(1)
    stats['mean_aln_score', 'filtered']=np.mean(mean_aln_score_filtered)
    stats['mean_aln_score', 'pass']=np.mean(mean_aln_score_pass)
    stats['umi', 'pass']=len(stats['umi', 'pass'])
    with open(out_prefix+'.stats.tsv', 'w') as out:
        print('\t'.join(['category', 'key', 'value']), file=out)
        for a,b in stats:
            print('\t'.join([str(x) for x in [a,b,stats[a,b]]]), file=out)
    with open(out_prefix+'.srbc_stats.tsv', 'w') as out:
        print('\t'.join(['srbc', 'filtered', 'expected', 'known', 'count']), file=out)
        for (a,b,c,d),e in srbc_stats.items():
            print('\t'.join([str(x) for x in [a,b,c,d,e]]), file=out)
    print("all done.")

def fix_tailor_bam(bam_file, outdir):
    """ Extract info from reads """
    stats=Counter()
    samfile = pysam.AlignmentFile(bam_file, "rb")# @UndefinedVariable
    reflen={k:v for k,v in zip(samfile.header.references, samfile.header.lengths)}
    out_file=outdir+'/'+Path(bam_file).stem+'_fixed.bam'
    samout = pysam.AlignmentFile(out_file, "wb", template=samfile ) 
    for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
        max_pos=reflen[read.reference_name]
        if read.reference_end < max_pos:
            samout.write(read)
        else:
            stats['filtered']+=1
    samout.close()
    try:
        pysam.index(out_file) # @UndefinedVariable
    except Exception as e:
        print("error sorting+indexing bam: %s" % e)
    print(stats)


def is_spikein_read(spikein, read, max_mm=0):
    """ tests whether the passed read stems most likely from this spikein """
    if max_mm==0:
        return spikein in read #exact match
    aln = pairwise2.align.globalxs( # globalxs(sequenceA, sequenceB, open, extend) -> alignments
            spikein, 
            read, 
            -2, 
            -1,
            penalize_end_gaps=(False,False), # penalize starting/ending gaps
            score_only=False,
            one_alignment_only=True)[0]
    _, _, aln_score, _, _ = aln
    mismatches=len(spikein)-aln_score
    if mismatches>max_mm:
        return False
    return True
    

def count_spikein_reads(fq_file, config, outdir, config_prefix=[]):
    """ Counts and filters spikeins reads """
    spikein_fa=get_config(config, config_prefix+['spikein_fa'], required=False)
    spikein=pysam.FastaFile(spikein_fa) # @UndefinedVariable
    spikeins={ref:''.join([x for x in spikein.fetch(ref) if x!='N']) for ref in spikein.references}
    spikein_seeds={ref:(seq[:math.floor(len(seq)/2)], seq[math.floor(len(seq)/2):] ) for ref,seq in spikeins.items()}
    max_mm=get_config(config, config_prefix+['max_mm'], default_value=1, required=False) # maximum mm 
    # iterate reads       
    si_counts=Counter()
    si_umi_counts=Counter()
    out_file_pass=outdir+'/'+Path(fq_file).stem+'.pass.fq'
    out_file_spikein=outdir+'/'+Path(fq_file).stem+'.spikein.fq'
    out_file_stats=outdir+'/'+Path(fq_file).stem+'.counts_spikein.tsv'
    with open(out_file_pass, 'w') as out_pass:                 
        with open(out_file_spikein, 'w') as out_spikein:                 
            with open(out_file_stats, 'w') as out_stats:
                print('\t'.join(['si_name', 'si_seq', 'si_counts', 'si_max_umi']), file=out_stats)  
                with gzip.open(fq_file, 'r') as in1:
                    it1=grouper(in1, 4, '')
                    for read1 in it1:
                        rn1, rs1, rc1, rq1 = [x.decode("utf-8").strip() for x in read1]
                        found=False
                        for si_name,(si_seed1, si_seed2) in spikein_seeds.items():
                            if (si_seed1 in rs1) or (si_seed2 in rs1):
                                if is_spikein_read(spikeins[si_name], rs1, max_mm=max_mm):
                                    si_counts[si_name]+=1
                                    si_umi_counts[si_name,rs1]+=1
                                    found=True
                        if not found:
                            si_counts['NA']+=1 # count 'pass' reads
                        print('\n'.join([rn1, rs1, rc1, rq1]), file=out_pass if not found else out_spikein)
                    for si_name, si_seq in spikeins.items():
                        si_umis=[ value for key, value in si_umi_counts.items() if key[0]==si_name]
                        max_umi=max(si_umis) if len(si_umis)>0 else 'NA'
                        print('\t'.join([str(x) for x in [si_name, si_seq, si_counts[si_name], max_umi]]), file=out_stats)
                    print('\t'.join([str(x) for x in ['NA', 'NA', si_counts['NA'], 'NA']]), file=out_stats)

# Supported counting profiles
supported_profiles={
    'miRNA_profile':{
        'fp_pre_tolerance': 5,
        'fp_ext_tolerance': 5,
        'tp_ext_tolerance': 5
        },
    'general_profile': {
        }
    }
# for debug output
read_colors={
    'wrong_strand':     '200,200,200',
    'fp_pre_tolerance': '125,0,0',
    'fp_ext_tolerance': '0,125,0',
    'tp_pre_tolerance': '0,0,125',
    'tp_ext_tolerance': '125,0,125'   
    }

def count_srna_reads(name, bam_file, anno_file, config, outdir, config_prefix=[]):
    """ Counts and filters srna reads """
    features={f:profile for f,profile in get_config(config, config_prefix+['features'], required=True).items()}
    gene_id_name=get_config(config, config_prefix+['gene_id'], required=True)
    extra_attributes=get_config(config, config_prefix+['extra_attributes'], default_value=[]).split(',')
    anno_fmt='GTF' if '.gtf' in anno_file else 'GFF3'

    out_file_cnt=outdir+'/'+name+'.counts.tsv'
    out_file_iso=outdir+'/'+name+'.iso.tsv'
    out_file_tai=outdir+'/'+name+'.tails.tsv'
    out_file_met=outdir+'/'+name+'.meta.tsv'
    out_file_stat=outdir+'/'+name+'.stats.tsv'
    out_file_bam_prefix=outdir+'/'+name+'.debug.bam' if get_config(config, config_prefix+['write_bam']) is not None else None

    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)# @UndefinedVariable    
    stats=Counter()
    mapped_read_names=set()
    for read in samfile.fetch():
        stats['alignments']+=1
        mapped_read_names.add(read.query_name)
    stats['mapped_reads']=len(mapped_read_names)
    samout = { f:pysam.AlignmentFile(out_file_bam_prefix+'.'+f+'.debug.bam', "wb", template=samfile ) for f in features.keys() } if out_file_bam_prefix is not None else None
    counted_read_names=set()
    last_chrom=None
    with open(out_file_cnt, 'w') as out_cnt:
        print('\t'.join(['gene_id', 'type', 'isoforms', 'reads', 'reads_norm', 'reads_tailed', 'reads_tailed_norm', 'reads_invalid']), file=out_cnt)
        with open(out_file_iso, 'w') as out_iso:
            print('\t'.join(['gene_id', 'type', 'iso_fpend', 'reads', 'reads_norm', 'reads_tailed', 'reads_tailed_norm']), file=out_iso)
            with open(out_file_tai, 'w') as out_tai:
                print('\t'.join(['gene_id', 'type', 'pos_fpend', 'tail', 'read_count']), file=out_tai)
                with open(out_file_met, 'w') as out_met:
                    print('\t'.join(['gene_id', 'type', 'chromosome', 'start', 'end', 'strand'] + extra_attributes), file=out_met)
                    for gtf_line in tqdm.tqdm(pysam.TabixFile(anno_file, mode="r").fetch(parser=pysam.asTuple())): # @UndefinedVariable
                        # get feature data
                        achrom,atype,astart,aend,astrand,ainfo=gtf_line[0], gtf_line[2], int(gtf_line[3]), int(gtf_line[4]),gtf_line[6], parse_info(gtf_line[8], fmt=anno_fmt)
                        gene_id=ainfo[gene_id_name]
                        if atype not in features.keys():
                            stats['skipped_uncounted_feature']+=1
                            continue
                        if features[atype] not in supported_profiles.keys():
                            stats['skipped_unsupported_profile']+=1
                            continue
                        profile=supported_profiles[features[atype]]
                        fp_pre_tolerance=get_config(profile, ['fp_pre_tolerance']) # None: will not be checked
                        fp_ext_tolerance=get_config(profile, ['fp_ext_tolerance'])
                        tp_pre_tolerance=get_config(profile, ['tp_pre_tolerance'])
                        tp_ext_tolerance=get_config(profile, ['tp_ext_tolerance'])                                    
                        # add counts to overall stats  
                        if achrom!=last_chrom :        
                            stats['counted_reads']+=len(counted_read_names)
                            counted_read_names=set()
                            last_chrom=achrom
                        feature_stats=Counter()
                        tails=Counter()
                        iso=Counter()
                        for read in samfile.fetch(achrom,astart,aend,until_eof=True):
                            feature_stats['mapped_reads']+=1
                            # check read validity
                            rstart,rend,rstrand=read.reference_start,read.reference_end,'-' if read.is_reverse else '+'
                            if rstrand != astrand: # check read strand
                                # if (samout is not None) and (samout[atype] is not None):
                                #     read.set_tag(tag='YC', value=read_colors['wrong_strand'], value_type="Z")
                                #     read.set_tag(tag='yc', value='wrong_strand', value_type="Z")
                                #     samout[atype].write(read)
                                feature_stats['wrong_strand']+=1
                                continue
                            # get 3'/5' ends
                            rfpend,rtpend = (rstart+1,rend) if rstrand=='+' else (rend,rstart+1)
                            afpend,atpend = (astart,aend) if astrand=='+' else (aend,astart)
                            fp_dist=(afpend-rfpend) if rstrand=='-' else (rfpend-afpend) # pos value: read starts *after* anntation, neg value: read starts before annotation
                            tp_dist=(atpend-rtpend) if rstrand=='-' else (rtpend-atpend) # pos value: read starts *after* annotation, neg value: read starts before annotation
                            # check 5' end of read
                            if fp_pre_tolerance and fp_dist<0 and -fp_dist>fp_pre_tolerance:
                                feature_stats['fp_pre_tolerance']+=1 # read starts to far before 5'-end of annotatione
                                if (samout is not None) and (samout[atype] is not None):
                                    read.set_tag(tag='YC', value=read_colors['fp_pre_tolerance'], value_type="Z")
                                    read.set_tag(tag='yc', value='fp_pre_tolerance', value_type="Z")
                                    samout[atype].write(read)
                                continue
                            if fp_ext_tolerance and fp_dist>0 and fp_dist>fp_ext_tolerance:
                                feature_stats['fp_ext_tolerance']+=1 # read starts to far after 5'-end of annotation
                                if (samout is not None) and (samout[atype] is not None):
                                    read.set_tag(tag='YC', value=read_colors['fp_ext_tolerance'], value_type="Z")
                                    read.set_tag(tag='yc', value='fp_ext_tolerance', value_type="Z")
                                    samout[atype].write(read)
                                continue
                            # check 3' end of read
                            if tp_pre_tolerance and tp_dist<0 and -tp_dist>tp_pre_tolerance:
                                feature_stats['tp_pre_tolerance']+=1 # read ends to far before 3'-end of annotatione
                                if (samout is not None) and (samout[atype] is not None):
                                    read.set_tag(tag='YC', value=read_colors['tp_pre_tolerance'], value_type="Z")
                                    read.set_tag(tag='yc', value='tp_pre_tolerance', value_type="Z")
                                    samout[atype].write(read)
                                continue
                            if tp_ext_tolerance and tp_dist>0 and tp_dist>tp_ext_tolerance:
                                feature_stats['tp_ext_tolerance']+=1 # read ends to far after 35'-end of annotation
                                if (samout is not None) and (samout[atype] is not None):
                                    read.set_tag(tag='YC', value=read_colors['tp_ext_tolerance'], value_type="Z")
                                    read.set_tag(tag='yc', value='tp_ext_tolerance', value_type="Z")
                                    samout[atype].write(read)
                                continue
                            # calculate multimapper normalized counts
                            nh=read.get_tag('NH') if read.has_tag('NH') else 1
                            # count read per 5' isoform 
                            iso[rfpend,'reads']+=1
                            iso[rfpend, 'reads_norm']+=1/nh
                            feature_stats['reads']+=1
                            feature_stats['reads_norm']+=1/nh
                            counted_read_names.add(read.query_name)
                            # count tails
                            softclipped=get_softclip_seq(read)[0 if astrand=='-' else 1]
                            if softclipped is not None:
                                tailpos, tailseq=softclipped
                                tails[tailseq]+=1
                                iso[rfpend, 'reads_tailed']+=1
                                iso[rfpend, 'reads_tailed_norm']+=1/nh
                                feature_stats['reads_tailed']+=1
                                feature_stats['reads_tailed_norm']+=1/nh
                            # write read ?
                            #if (samout is not None) and (samout[atype] is not None):
                            #    samout[atype].write(read)                               
                        isoform_5pends=sorted(set([x for x,_ in iso.keys()]))               
                        # write results per anno
                        print('\t'.join([str(x) for x in [gene_id, atype, len(isoform_5pends), 
                                                 feature_stats['reads'], feature_stats['reads_norm'], # FIXME: these are alignment, not unique read stats 
                                                 feature_stats['reads_tailed'],feature_stats['reads_tailed_norm'],
                                                 feature_stats['mapped_reads']-feature_stats['reads']]]), file=out_cnt)                    
                        # isoforms: one line per found 5'-end
                        for fpend in isoform_5pends:
                            print('\t'.join([str(x) for x in [gene_id, atype, fpend, 
                                             iso[fpend, 'reads'], iso[fpend, 'reads_norm'], 
                                             iso[fpend, 'reads_tailed'],iso[fpend, 'reads_tailed_norm']]]), file=out_iso)
                        # write tailing stats
                        for tail,tailcount in tails.items():
                            print('\t'.join([str(x) for x in [gene_id, atype, fpend, tail, tailcount]]), file=out_tai)
                        # write meta data
                        print('\t'.join([str(x) for x in [gene_id, atype, achrom, astart, aend, astrand] + [ainfo[k] if k in ainfo else 'NA' for k in extra_attributes ] ] ), file=out_met)
    if samout is not None:
        for f in features.keys():
            samout[f].close()
            out_file_bam=out_file_bam_prefix+'.'+f+'.debug.bam'
            try:
                pysam.sort("-o", out_file_bam+'.tmp.bam', out_file_bam) # @UndefinedVariable
                os.replace(out_file_bam+'.tmp.bam', out_file_bam)
                pysam.index(out_file_bam) # @UndefinedVariable
            except Exception as e:
                print("error sorting+indexing bam: %s" % e)
    stats['counted_reads']+=len(counted_read_names)   
    with open(out_file_stat, 'w') as out_stat:
        print('\t'.join(['key', 'value']), file=out_stat)
        for k,v in stats.items():
            print(str(k)+'\t'+str(v), file=out_stat)


def downsample_per_chrom(bam_file, max_reads, outdir):
    """ Random subsampling to guarantee max_reads per chromosome """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)# @UndefinedVariable
    out_file_bam = bam_file+'.subsampled_max%i.bam' % max_reads
    samout =  pysam.AlignmentFile(out_file_bam, "wb", template=samfile )
    read_indices={i.contig:random.sample(range(i.mapped),min(max_reads, i.mapped)) for i in samfile.get_index_statistics()}
    for c in samfile.references:
        idx=0
        ri=set(read_indices[c])
        for read in samfile.fetch(c):
            if idx in ri:
                samout.write(read)    
            idx+=1
    samout.close()
    try:
        pysam.sort("-o", out_file_bam+'.tmp.bam', out_file_bam) # @UndefinedVariable
        os.replace(out_file_bam+'.tmp.bam', out_file_bam)
        pysam.index(out_file_bam) # @UndefinedVariable
    except Exception as e:
        print("error sorting+indexing bam: %s" % e)
        
        
        
def analyze_filtered_reads(dat_file, config, max_reads=1000):
    sample_name=Path(dat_file[:-len('.filtered.fq.gz')]).stem
    sample_sheet=pd.read_csv(get_config(config, ['sample_sheet'], required=True), sep='\t')
    srbcs={k:v.strip() for k,v in zip(sample_sheet['filename_prefix'],sample_sheet['sRBC']) }
    expected_srbc = srbcs[sample_name]
    anchor_seq=get_config(config, config_prefix+['anchor_seq'], required=True)
    umi_len=get_config(config, config_prefix+['umi_len'], required=True)
    srbc_len=get_config(config, config_prefix+['srbc_len'], required=True)
    print(f"expected_srbc: {expected_srbc}")
    cnt=Counter()
    with gzip.open(dat_file, 'r') as in1:
        it1=grouper(in1, 4, '')
        for read in tqdm.tqdm(it1):
            cnt['reads']+=1
            query_name, query_sequence, _, query_qualities = [x.decode("utf-8").strip() for x in read]
            # max alignmnet score of anchor
            aln = pairwise2.align.globalxs( # globalxs(sequenceA, sequenceB, open, extend) -> alignments
                    anchor_seq, 
                    query_sequence, 
                    -2, 
                    -1,
                    penalize_end_gaps=(False,False), # penalize starting/ending gaps
                    score_only=False,
                    one_alignment_only=True)[0]
            _, _, aln_score, _, _ = aln
            aln_score=aln_score/len(anchor_seq)
            if aln_score<=0.9:
                cnt['bad_anchor']+=1
                continue
            start=len(aln[0])-len(aln[0].lstrip('-'))
            query_sequence=query_sequence[:start]
            result=(expected_srbc in query_sequence, len(query_sequence)-umi_len-srbc_len>18)
            cnt[result]+=1
            if result==(True,True):
                print(query_name, query_sequence)
            if cnt['reads']>max_reads:
                break
    print(cnt)
            # 


usage = '''

  Copyright (C) 2022 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':

    parser = {}

    parser["build_transcriptome"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["build_transcriptome"].add_argument("-n", "--name", type=str, required=True, dest="name", metavar="name", help="file name prefix of the created transcriptome")
    parser["build_transcriptome"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["build_transcriptome"].add_argument("-p", "--config_prefix", type=str, required=False, default=None, dest="config_prefix", metavar="config_prefix", help="optional config file category")
    parser["build_transcriptome"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["parse_reads"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["parse_reads"].add_argument("-d", "--dat", type=str, required=True, dest="dat_file", metavar="dat_file", help="Input file, either unaligned BAM or FASTQ file")
    parser["parse_reads"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["parse_reads"].add_argument("-p", "--config_prefix", type=str, required=False, default=None, dest="config_prefix", metavar="config_prefix", help="optional config file category")
    parser["parse_reads"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["fix_tailor_bam"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["fix_tailor_bam"].add_argument("-b", "--bam", type=str, required=True, dest="bam_file", metavar="bam_file", help="Unaligned BAM input file")
    parser["fix_tailor_bam"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["count_spikein_reads"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["count_spikein_reads"].add_argument("-f", "--fq", type=str, required=True, dest="fq_file", metavar="fq_file", help="Preprocessed FQ file")
    parser["count_spikein_reads"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["count_spikein_reads"].add_argument("-p", "--config_prefix", type=str, required=False, default=None, dest="config_prefix", metavar="config_prefix", help="optional config file category")
    parser["count_spikein_reads"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["count_srna_reads"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["count_srna_reads"].add_argument("-n", "--name", type=str, required=True, dest="name", metavar="name", help="output file name prefix")
    parser["count_srna_reads"].add_argument("-b", "--bam", type=str, required=True, dest="bam_file", metavar="bam_file", help="Tailor BAM file")
    parser["count_srna_reads"].add_argument("-a", "--anno", type=str, required=True, dest="anno_file", metavar="anno_file", help="Annotation GTF or GFF file")
    parser["count_srna_reads"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["count_srna_reads"].add_argument("-p", "--config_prefix", type=str, required=False, default=None, dest="config_prefix", metavar="config_prefix", help="optional config file category")
    parser["count_srna_reads"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    parser["downsample_per_chrom"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["downsample_per_chrom"].add_argument("-b", "--bam", type=str, required=True, dest="bam_file", metavar="bam_file", help="Tailor BAM file")
    parser["downsample_per_chrom"].add_argument("-m", "--max_reads", type=int, required=True, dest="max_reads", metavar="max_reads", help="Maximum reads per chromosome")
    parser["downsample_per_chrom"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

    #============================================================================
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print("usage: annotate_shRNA.py [-h] " + ",".join(parser.keys()))
        sys.exit(1)
    mod = sys.argv[1]
    if mod not in parser.keys():
        print("Invalid module '%s' selected. Please use one of %s" % (mod, ",".join(parser.keys())))
        sys.exit(1)

    args = parser[mod].parse_args(sys.argv[2:])
    #============================================================================

    # output dir (current dir if none provided)
    outdir = os.path.abspath(args.outdir if args.outdir else os.getcwd())+'/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if mod == "build_transcriptome":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        build_transcriptome(args.name, config, outdir, config_prefix=args.config_prefix.split(',') if args.config_prefix is not None else [])

    if mod == "parse_reads":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        parse_reads(args.dat_file, config, outdir, config_prefix=args.config_prefix.split(',') if args.config_prefix is not None else [])
        
    if mod == "fix_tailor_bam":
        fix_tailor_bam(args.bam_file, outdir)
        
    if mod == "count_spikein_reads":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        count_spikein_reads(args.fq_file, config, outdir, config_prefix=args.config_prefix.split(',') if args.config_prefix is not None else [])

    if mod == "count_srna_reads":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        count_srna_reads(args.name, args.bam_file, args.anno_file, config, outdir, config_prefix=args.config_prefix.split(',') if args.config_prefix is not None else [])

    if mod =="downsample_per_chrom":
        downsample_per_chrom(args.bam_file, args.max_reads, outdir)
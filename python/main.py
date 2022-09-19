'''
@author: niko.popitsch@univie.ac.at
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os, sys, json, logging, random, time
from collections import Counter, OrderedDict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import pandas as pd
import numpy as np
import pysam
import tqdm
from pathlib import Path
from utils import *

def parse_read(read, config):
    read.query_sequence
    anchor_seq=get_config('anchor_seq', config)
    umi_len=get_config('umi_len', config)
    srbc_len=get_config('srbc_len', config)
    aln = pairwise2.align.globalxs( # globalxs(sequenceA, sequenceB, open, extend) -> alignments
            anchor_seq, 
            read.query_sequence, 
            -2, 
            -1,
            penalize_end_gaps=(False,False), # penalize starting/ending gaps
            score_only=False,
            one_alignment_only=True)[0]
    _, _, aln_score, _, _ = aln
    aln_score=aln_score/len(anchor_seq)
    if aln_score < get_config('min_aln_score', config, 0.9):
        return True, read, aln_score, None, None
    # get pos in read
    start=len(aln[0])-len(aln[0].lstrip('-'))
    end=len([x for x in aln[1][:len(aln[0].rstrip('-'))] if x != '-'])
    # extract umi, srbc and clip read
    srbc=read.query_sequence[start-srbc_len:start]
    umi=read.query_sequence[start-srbc_len-umi_len:start-srbc_len]
    query_sequence=read.query_sequence[end:]
    query_qualities=read.query_qualities[end:]
    if len(query_sequence)<get_config('min_read_len', config, 22):
        return True, read, aln_score, None, None # filter if remaining read too short
    read.query_sequence = query_sequence
    read.query_qualities = query_qualities
    read.query_name='_'.join([read.query_name, srbc, umi])
    return False, read, aln_score, umi, srbc

def parse_reads(bam_file, config, outdir):
    """ Extract info from reads """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)# @UndefinedVariable
    out_prefix=outdir+Path(bam_file).stem
    stats=Counter()
    stats['umi', 'pass']=set()
    mean_aln_score_filtered, mean_aln_score_pass=list(), list()
    with open(out_prefix+'.filtered.fq', 'w') as fq_filtered:
        with open(out_prefix+'.fq', 'w') as fq_pass:
            for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
                filtered, read, aln_score, umi, srbc=parse_read(read, config)
                if filtered:                    
                    out=fq_filtered 
                    stats['read_count', 'filtered']+=1                  
                    mean_aln_score_filtered.append(aln_score)                
                else:
                    out=fq_pass          
                    stats['read_count', 'pass']+=1           
                    stats['umi', 'pass'].add(umi)           
                    stats['srbc', srbc]+=1
                    mean_aln_score_pass.append(aln_score)                
                print('@%s\n%s\n+\n%s' % (read.query_name, 
                                               read.query_sequence, 
                                               ''.join(map(lambda x: chr( x+33 ), read.query_qualities)) if read.query_qualities is not None else '???'), 
                                               file=out)
                # if stats['read_count', 'filtered']>10:
                #     break
    samfile.close()
    stats['mean_aln_score', 'filtered']=np.mean(mean_aln_score_filtered)
    stats['mean_aln_score', 'pass']=np.mean(mean_aln_score_pass)
    stats['umi', 'pass']=len(stats['umi', 'pass'])
    with open(out_prefix+'.stats.tsv', 'w') as out:
        print('\t'.join(['category', 'key', 'value']), file=out)
        for a,b in stats:
            print('\t'.join([str(x) for x in [a,b,stats[a,b]]]), file=out)
    print("all done.")

usage = '''

  Copyright (C) 2022 XXX.  All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
'''
if __name__ == '__main__':

    parser = {}

    parser["parse_reads"] = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser["parse_reads"].add_argument("-b", "--bam", type=str, required=True, dest="bam_file", metavar="bam_file", help="Unaligned BAM input file")
    parser["parse_reads"].add_argument("-c", "--config", type=str, required=True, dest="config_file", metavar="config_file", help="JSON config file")
    parser["parse_reads"].add_argument("-o", "--outdir", type=str, required=False, dest="outdir", metavar="outdir", help="output directory (default is current dir)")

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

    if mod == "parse_reads":
        # load and check onfig
        config=json.load(open(args.config_file), object_pairs_hook=OrderedDict)
        parse_reads(args.bam_file, config, outdir)

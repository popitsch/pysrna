from itertools import zip_longest
import pysam
import os, sys
import pybedtools

rcmap=bytes.maketrans(b'ATCGatcgNn', b'TAGCTAGCNN')
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    return seq[::-1].translate(rcmap)

def get_config(config, keys, default_value=None, required=False): 
    """ gets value from dict of dicts or either returns default_value if path does not exist (required=False) or throws an error.
    """
    d = config
    for k in keys:
        if k not in d:
            assert (not required), 'Mandatory config path "%s" missing' % ' > '.join(keys)
            return default_value
        d = d[k]
    return d  

def parse_info(info, fmt='gff3'):
    """ parse GFF3/GTF info section """
    if fmt=='GTF':
        return {k:v.translate({ord(c): None for c in '"'}) for k,v in [a.strip().split(' ') for a in info.split(';') if ' ' in a]}
    return {k:v for k,v in [a.split('=') for a in info.split(';') if '=' in a]}

def grouper(iterable, n, fillvalue=None):
    """ Groups n lines into a list """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def sort_bgzip_and_tabix(in_file, out_file=None, create_index=True, del_uncompressed=True, seq_col=0, start_col=1, end_col=1, line_skip=0, zerobased=False):
    """ Will BGZIP the passed file and creates a tabix index with the given params if create_index is True"""
    if out_file is None:
        out_file=in_file+'.gz'
    pre,post=os.path.splitext(in_file)
    sorted_in_file=pre+'.sorted'+post
    pybedtools.BedTool(in_file).sort().saveas(sorted_in_file)    
    pysam.tabix_compress(sorted_in_file, out_file, force=True) # @UndefinedVariable
    if create_index:
        pysam.tabix_index(out_file, force=True, seq_col=seq_col, start_col=start_col, end_col=end_col, meta_char='#', line_skip=line_skip, zerobased=zerobased) # @UndefinedVariable
    if del_uncompressed:
        os.remove(in_file)
    os.remove(sorted_in_file) # removed sorted version
        
def format_fasta(string, ncol=80):
    return '\n'.join(string[i:i+ncol] for i in range(0, len(string), ncol))

def get_softclip_seq(r):
    left,right=None,None
    pos=0
    for i,(op,l) in enumerate(r.cigartuples):
        if (i==0) & (op==4):
            left=(r.reference_start+pos, r.query_sequence[pos:pos+l])
        if (i==len(r.cigartuples)-1)& (op==4):
            right=(r.reference_start+pos, r.query_sequence[pos:pos+l])
        pos+=l
    return left,right
from itertools import zip_longest

rcmap=bytes.maketrans(b'ATCGatcgNn', b'TAGCTAGCNN')
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    return seq[::-1].translate(rcmap)

def get_config(config, keys, default_value, required=False): 
    """ gets value from dict of dicts or either returns default_value if path does not exist (required=False) or throws an error.
    """
    d = config
    for k in keys:
        if k not in d:
            assert (not required), 'Mandatory config path "%s" missing' % ' > '.join(keys)
            return default_value
        d = d[k]
    return d  

def parse_info(info):
    """ parse GFF3 info section """
    return {k:v for k,v in [a.split('=') for a in info.split(';') if '=' in a]}

def grouper(iterable, n, fillvalue=None):
    """ Groups n lines into a list """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
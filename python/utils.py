from itertools import zip_longest

rcmap=bytes.maketrans(b'ATCGatcgNn', b'TAGCTAGCNN')
def reverse_complement(seq):
    """ Calculate reverse complement DNA sequence """
    return seq[::-1].translate(rcmap)

def get_config(prop, config, default=None):
    if prop not in config:
        assert (default is not None), 'Mandatory config property %s missing' % prop
    return config.get(prop, default)

def parse_info(info):
    """ parse GFF3 info section """
    return {k:v for k,v in [a.split('=') for a in info.split(';') if '=' in a]}

def grouper(iterable, n, fillvalue=None):
    """ Groups n lines into a list """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
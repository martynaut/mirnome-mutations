import requests
from xml.etree import ElementTree
import regex

rev_comp = {
    'G': 'C',
    'C': 'G',
    'A': 'U',
    'T': 'A',
    'U': 'A'
}

iupac_dict = {
    'R': '[GA]{1}',
    'Y': '[UC]{1}',
    'M': '[AC]{1}',
    'K': '[GU]{1}',
    'S': '[GC]{1}',
    'W': '[AU]{1}',
    'H': '[ACU]{1}',
    'B': '[GUC]{1}',
    'V': '[GCA]{1}',
    'D': '[GUA]{1}',
    'N': '[GCUA]{1}',
}

motifs = {
    'UAGGGAW': 'loop',
    'AUUUUUAUUUU': 'loop',
    'GGGU': 'loop',
    'GGAG': 'loop',
    'YGCY': 'loop',
    'UGC': 'loop',
    'UGU': 'loop',
    'AUCUU': 'loop',
    'SMUANY': 'loop',
    'CAUC': 'loop',
    'UAUAA': 'loop',
    'UUUUUCCNUCUUU': 'loop',
    'VCAUCH': 'all',
    'GCAUG': 'all',
    'CAGAC': 'all',
    'UGUNNNNNNNUGU': 'all',
    'GCAGCGC': 'all'
}


def motif_check(row_values, regex_motif, raw_motif):
    if motifs[raw_motif] == 'loop':
        if row_values['type'] != 'loop':
            return '-'
    all_before = regex.finditer(regex_motif, row_values['whole_seq'], overlapped=True)
    for i in all_before:
        motif_start, motif_stop = i.span()
        if row_values['orientation'] == '+':
            if (row_values['pos'] < motif_stop + row_values[1]) \
                    and (row_values['pos'] > motif_start + row_values[1]):
                return 'loss'
        else:
            if (row_values['pos'] < row_values[2] - motif_start) \
                    and (row_values['pos'] > row_values[2] - motif_stop):
                return 'loss'
    return '-'


def change_in_motifs(row_value, motifs_columns):
    for col in motifs_columns:
        if row_value[col] != '-':
            return 1
    return 0


def get_whole_sequence(row):
    URL = "http://genome.ucsc.edu/cgi-bin/das/hg38/dna"
    PARAMS = {'segment': str(row[0]) + ':' + str(row[1]+1) + ',' + str(row[2])}
    r = requests.get(url=URL, params=PARAMS)
    tree = ElementTree.fromstring(r.content)
    seq = []
    for child in tree:
        for cc in child:
            seq.append(cc.text.upper())
    if len(seq) == 1:
        return seq[0].replace(' ', '').replace('\n', '').strip()
    else:
        raise ValueError

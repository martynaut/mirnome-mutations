import hgvs.parser
import hgvs.exceptions


def var_c_p_prep(row):
    hp = hgvs.parser.Parser()

    d = {'chr1': '.11', 'chr2': '.12', 'chr3': '.12', 'chr4': '.12', 'chr5': '.10', 'chr6': '.12', 'chr7': '.14',
         'chr8': '.11', 'chr9': '.12', 'chr10': '.11', 'chr11': '.10', 'chr12': '.12', 'chr13': '.11', 'chr14': '.9',
         'chr15': '.10', 'chr16': '.10', 'chr17': '.11', 'chr18': '.10', 'chr19': '.10', 'chr20': '.11', 'chr21': '.9',
         'chr22': '.11', 'chr23': '.11', 'chr24': '.10'}

    if row['chrom'] == 'chrX':
        hgvs_g = 'NC_000023' + d['chr23'] + ':' + row['to_hgvs_g']
    elif row['chrom'] == 'chrY':
        hgvs_g = 'NC_000024' + d['chr24'] + ':' + row['to_hgvs_g']
    elif len(row['chrom'].replace('chr', '')) == 1:
        hgvs_g = 'NC_00000' + row['chrom'].replace('chr', '') + d[row['chrom']] + ':' + row['to_hgvs_g']
    else:
        hgvs_g = 'NC_0000' + row['chrom'].replace('chr', '') + d[row['chrom']] + ':' + row['to_hgvs_g']

    try:
        var_g = hp.parse_hgvs_variant(hgvs_g)
    except hgvs.exceptions.HGVSParseError:
        return ""

    return str(var_g)


def rev_comp(letter):
    nt_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'U': 'A'
    }
    if len(letter) == 1:
        return nt_dict[letter]
    elif ',' in letter:
        words = letter.split(',')
        new_words = []
        for word in words:
            new_words.append(''.join([nt_dict[x] for x in word]))
        return ','.join(new_words)
    else:
        letter = [nt_dict[x] for x in letter]
        return ''.join(letter)


def tryconvert(value):
    try:
        return str(int(value))
    except ValueError:
        return str(value)

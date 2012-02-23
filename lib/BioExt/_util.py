
from collections import UserList
from copy import deepcopy
from itertools import product
from random import randint, uniform

try:
    from scipy.stats import norm
except:
    norm = None

from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord


__version__ = '0.0.9'


__all__ = [
    '_GAP', '_STOP',
    'randgene',
    'untranslate',
    'pyro_errors',
    'homosplit',
    'intersperse',
    'by_codon',
    'enumerate_by_codon',
    'AmbigList',
    'translate_ambiguous'
] + ['errorize'] if norm else []


_GAP = '-'
_STOP = '*'


_NUC_AMBIGS = {
    'W': 'AT',
    'S': 'CG',
    'M': 'AC',
    'K': 'GT',
    'R': 'AG',
    'Y': 'CT',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}


_AMINO_AMBIGS = {
    'B': 'DN',
    'J': 'IL',
    'Z': 'EQ',
    'X': 'ACDEFGHIKLMNPQRSTVWY'
}


def randgene(length, ppf):
    s = []
    avoid = c = n1 = n2 = -1
    l = 0
    while l < length:
        # get the length
        lp = ppf(uniform(0, 1)) + 1
        # avoid stop codons
        if l % 3 == 1 and n1 == 3 and lp > 1:
            avoid = 0 # TAA (3, 0, 0)
        elif l % 3 == 2 and n2 == 3:
            if n1 == 0:
                avoid = 2 # TAG (3, 0, 2)
            elif n1 == 2:
                avoid = 0 # TGA (3, 2, 0)
        # avoid the last character and stop codons
        while c == n1 or c == avoid:
            c = randint(0, 3)
        # setup state
        n2 = n1 # negative 1 to negative 2
        n1 = c # negative 1 to char
        avoid = -1 # reset avoid
        # grow seq
        s.append('ACGT'[c] * lp)
        l += lp
    return ''.join(s)[:length]


def untranslate(seq):
    d = {}
    for a in 'ACGT':
        for b in 'ACGT':
            for c in 'ACGT':
                cdn = a + b + c
                aa = translate(cdn)
                if aa not in d:
                    d[aa] = []
                d[aa].append(cdn)
    r = []
    for aa in seq:
        o = d[aa]
        r.append(o[randint(0, len(o) - 1)])
    return ''.join(r)


# take from:
# Characteristics of 454 pyrosequencing data—enabling realistic simulation with flowsim
# Susanne Balzer, Ketil Malde, Anders Lanzén, Animesh Sharma, and Inge Jonassen
# Bioinformatics Vol. 26 ECCB 2010, pages i420–i425
# doi:10.1093/bioinformatics/btq365
def pyro_errors(x):
    if x < 0:
        raise ValueError('x must be >= 0')
    elif x < 6:
        return [(0.1230, 0.0737),
                (1.0193, 0.1227),
                (2.0006, 0.1585),
                (2.9934, 0.2188),
                (3.9962, 0.3168),
                (4.9550, 0.3863)][x]
    else:
        return (x, 0.03494 + x * 0.0685)


if norm:
    def errorize(sequence, randlen=None):
        if not isinstance(sequence, str):
            raise ValueError('sequence must be of type str')
        if randlen is None:
            randlen = lambda l: round(norm.ppf(uniform(0, 1), *pyro_errors(l)))
        alph = tuple(set(sequence))
        l0 = len(alph) - 1
        r = [alph[randint(0, l0)] * randlen(0)]
        for x in homosplit(sequence):
            r.append(x[0] * randlen(len(x)))
            r.append(alph[randint(0, l0)] * randlen(0))
        return ''.join(r)


def homosplit(iterable):
    it = iter(iterable)
    a = next(it)
    i = 1
    for b in it:
        if b != a:
            yield a * i
            a = b
            i = 0
        i += 1
    yield a * i


def intersperse(iterable, delimiter, n=1):
    if n < 1:
        msg = "cannot intersperse every n = '%d' < 1 elements" % n
        raise ValueError(msg)
    n -= 1 # we yield a value manually, so yield in n - 1 groups
    it = iter(iterable)
    yield next(it) # yield manually
    while True:
        for _ in range(n):
            yield next(it)
        # grab a value manually before
        # yielding a delimiter to avoid
        # terminal delimiters
        saved = next(it)
        yield delimiter
        yield saved # yield manually


def by_codon(seq, gap_char=_GAP):
    for _, cdn in enumerate_by_codon(seq, gap_char):
        yield cdn


def enumerate_by_codon(seq, gap_char=_GAP):
    if isinstance(seq, SeqRecord):
        seq = seq.seq.tostring()
    elif isinstance(seq, Seq):
        seq = seq.tostring()
    elif not isinstance(seq, str):
        raise ValueError('can only enumerate codons of a SeqRecord, Seq, or str')

    seqlen = len(seq)
    num_cdns = seqlen // 3
    for i in range(num_cdns):
        pos = 3 * i
        cdn = seq[pos:min(seqlen, pos + 3)]
        cdn += gap_char * (3 - len(cdn))
        yield (pos, cdn)


class AmbigList(UserList):
    def __init__(self, aminos):
        super(AmbigList, self).__init__(aminos)


def translate_ambiguous(seq, gap_char=_GAP, trim_gaps=True):
    if isinstance(seq, SeqRecord):
        seqstr = seq.seq.tostring()
    elif isinstance(seq, Seq):
        seqstr = seq.tostring()
    elif not isinstance(seq, str):
        raise ValueError('can only enumerate codons of a SeqRecord, Seq, or str')

    if trim_gaps:
        seqstr = seqstr.replace(gap_char, '')
    seqstr = seqstr.upper()

    aminos = []
    gap_cdn = 3 * gap_char
    for _, cdn in enumerate_by_codon(seqstr, gap_char):
        # if we're not trimming gaps,
        # convert gap codons into single codons
        if cdn == gap_cdn:
            aminos.append(set('-'))
            continue
        # otherwise, combinatorial fun
        nucs = []
        for nuc in cdn:
            if nuc in _NUC_AMBIGS:
                nucs.append(_NUC_AMBIGS[nuc])
            else:
                nucs.append(nuc)
        aminos.append(set(translate(''.join(p)) for p in product(*nucs)))

    return AmbigList(aminos)

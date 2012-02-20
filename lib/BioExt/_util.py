
from collections import UserList
from copy import deepcopy
from itertools import product

from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord


__version__ = '0.0.9'


__all__ = [
    '_GAP', '_STOP',
    'by_codon',
    'enumerate_by_codon',
    'AmbigList',
    'translate_ambiguous'
]


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

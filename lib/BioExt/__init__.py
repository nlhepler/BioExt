
from copy import deepcopy
from itertools import product

from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord


__version__ = '0.0.9'


__all__ = [
    'by_codon',
    'enumerate_by_codon',
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


def by_codon(seq, gap_char=_GAP):
    for _, cdn in enumerate_by_codon(seq, gap_char):
        yield cdn


def enumerate_by_codon(seq, gap_char=_GAP):
    if isinstance(seq, SeqRecord):
        seq = seq.seq.data
    elif isinstance(seq, Seq):
        seq = seq.data
    elif not isinstance(seq, str):
        raise ValueError('can only enumerate codons of a SeqRecord, Seq, or str')

    seqlen = len(seq)
    num_cdns = seqlen // 3
    for i in range(num_cdns):
        pos = 3 * i
        cdn = seq[pos:min(seqlen, pos + 3)]
        cdn += gap_char * (3 - len(cdn))
        yield (pos, cdn)


def translate_ambiguous(seq):
    if isinstance(seq, SeqRecord):
        seqstr = seq.seq.data
    elif isinstance(seq, Seq):
        seqstr = seq.data
    elif not isinstance(seq, str):
        raise ValueError('can only enumerate codons of a SeqRecord, Seq, or str')

    seqstr = seqstr.replace(_GAP, '')

    aminos = []
    for _, cdn in enumerate_by_codon(seqstr):
        nucs = []
        for nuc in cdn:
            if nuc in _NUC_AMBIGS:
                nucs.append(_NUC_AMBIGS[nuc])
            else:
                nucs.append(nuc)
        aminos.append([translate(''.join(p)) for p in product(*nucs)])

    return aminos

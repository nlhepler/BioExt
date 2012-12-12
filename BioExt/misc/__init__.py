# -*- coding: utf-8 -*-

from __future__ import division, print_function

try:
    from UserList import UserList
except:
    from collections import UserList

from itertools import product
from random import randint, random

from Bio.Seq import Seq, translate as _translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature


__all__ = [
    '_GAP', '_STOP',
    'randgene',
    'homosplit',
    'intersperse',
    'by_codon',
    'enumerate_by_codon',
    'AmbigList',
    'translate_ambiguous',
    'translate',
    'compute_cigars'
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


def randgene(length, ppf):
    s = []
    avoid = c = n1 = n2 = -1
    l = 0
    while l < length:
        # get the length
        lp = round(ppf(random())) + 1
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


def _translate_gapped(seq, *args, **kwds):
    if isinstance(seq, SeqRecord):
        s = str(seq.seq)
    elif isinstance(seq, Seq):
        s = str(seq)
    elif isinstance(seq, str):
        s = seq
    else:
        raise ValueError("can only translate sequences of type SeqRecord, Seq, or str")
    gaps = 0
    lwr = 0
    protein = ''
    for i in range(0, len(s), 3):
        j = min(i + 3, len(s))
        if s[i:j] == '---'[:j - i]:
            if not gaps:
                protein += _translate(s[lwr:i].replace('-', 'N'))
            gaps += 1
        elif gaps:
            protein += '-' * gaps
            gaps = 0
            lwr = i
    if gaps:
        protein += '-' * gaps
    else:
        protein += _translate(s[lwr:len(seq)].replace('-', 'N'))
    return protein


def translate(seq, *args, **kwds):
    if isinstance(seq, SeqRecord):
        features = []
        for feature in seq.features:
            # ignore features on the opposite strand
            strand = feature.location.strand
            if strand is not None and strand < 0:
                continue
            features.append(
                SeqFeature(
                    FeatureLocation(
                        feature.location.start.position // 3,
                        feature.location.end.position // 3
                        # for protein sequences,
                        # having a strand doesn't make sense
                    ),
                    type=feature.type
                )
            )
        return SeqRecord(
            Seq(_translate_gapped(seq.seq, *args, **kwds)),
            id=seq.id,
            name=seq.name,
            description=seq.description,
            features=features
        )
    elif isinstance(seq, Seq):
        return Seq(_translate_gapped(seq, *args, **kwds))
    elif isinstance(seq, str):
        return _translate_gapped(seq, *args, **kwds)
    else:
        raise ValueError('can only translate sequences of type SeqRecord, Seq, or str')


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
        aminos.append(set(_translate(''.join(p)) for p in product(*nucs)))

    return AmbigList(aminos)


def compute_cigars(alignment, reference):
    ncol = alignment.get_alignment_length()
    if len(reference) != ncol:
        raise ValueError('Reference must be the same length as the alignment')

    GAPS = ('-', '.')
    for record in alignment:
        # find start, end of record in the ref
        start, end = 0, ncol - 1
        for i in range(ncol):
            if record[i] not in GAPS:
                start = i
                break
        for i in range(ncol, 0, -1):
            if record[i - 1] not in GAPS:
                end = i
                break

        # compute cigar and edit distance
        cigar = ''
        edit_distance = 0
        mode = ''
        count = 0
        for i in range(start, end):
            ref = reference[i].upper()
            query = record[i].upper()
            if ref in GAPS and query in GAPS:
                # if both are gaps, skip
                pass
            elif ref in GAPS:
                m = 'I' # insertion
                edit_distance += 1
            elif query in GAPS:
                m = 'D' # deletion
                edit_distance += 1
            elif ref == query:
                m = '=' # match
            else:
                m = 'X' # mismatch
                edit_distance += 1
            # cigar handling
            if not mode:
                mode = m
                count = 1
            elif m == mode:
                count += 1
            else:
                cigar += '%d%s' % (count, mode)
                mode = m
                count = 1
        # don't forget the last bit
        cigar += '%d%s' % (count, mode)

        # inject the annotations and yield
        record.annotations['CIGAR'] = cigar
        record.annotations['edit_distance'] = edit_distance
        record.annotations['position'] = start + 1 # 1-indexed
        record.annotations['length'] = end - start
        record.annotations['reference_name'] = reference.name

        yield record

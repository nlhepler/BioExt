# -*- coding: utf-8 -*-

from __future__ import division, print_function

try:
    from UserList import UserList
except:
    from collections import UserList

from copy import copy
from itertools import groupby, product
from random import randint, random
from re import compile as re_compile

from Bio.Seq import Seq, translate as _translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature


__all__ = [
    '_GAP', '_STOP',
    'randgene',
    'homocount',
    'homosplit',
    'intersperse',
    'by_codon',
    'enumerate_by_codon',
    'AmbigList',
    'translate_ambiguous',
    'translate',
    'compute_cigar',
    'gapless',
    'gapful',
    'labelstream',
    'clip'
    ]


_GAP = '-'
_GAPS = (_GAP, '.')
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
            avoid = 0  # TAA (3, 0, 0)
        elif l % 3 == 2 and n2 == 3:
            if n1 == 0:
                avoid = 2  # TAG (3, 0, 2)
            elif n1 == 2:
                avoid = 0  # TGA (3, 2, 0)
        # avoid the last character and stop codons
        while c == n1 or c == avoid:
            c = randint(0, 3)
        # setup state
        n2 = n1  # negative 1 to negative 2
        n1 = c  # negative 1 to char
        avoid = -1  # reset avoid
        # grow seq
        s.append('ACGT'[c] * lp)
        l += lp
    return ''.join(s)[:length]


def homocount(iterable):
    it = iter(iterable)
    a = next(it)
    i = 1
    for b in it:
        if b != a:
            yield (a, i)
            a = b
            i = 0
        i += 1
    yield (a, i)


def homosplit(iterable):
    for a, i in homocount(iterable):
        yield a * i


def intersperse(iterable, delimiter, n=1):
    if n < 1:
        msg = "cannot intersperse every n = '%d' < 1 elements" % n
        raise ValueError(msg)
    n -= 1  # we yield a value manually, so yield in n - 1 groups
    it = iter(iterable)
    yield next(it)  # yield manually
    while True:
        for _ in range(n):
            yield next(it)
        # grab a value manually before
        # yielding a delimiter to avoid
        # terminal delimiters
        saved = next(it)
        yield delimiter
        yield saved  # yield manually


def by_codon(seq, gap_char=_GAP):
    for _, cdn in enumerate_by_codon(seq, gap_char):
        yield cdn


def enumerate_by_codon(seq, gap_char=_GAP):
    if isinstance(seq, SeqRecord):
        seq = seq.seq.tostring()
    elif isinstance(seq, Seq):
        seq = seq.tostring()
    elif not isinstance(seq, str):
        msg = 'can only enumerate codons of a SeqRecord, Seq, or str'
        raise ValueError(msg)

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
        msg = "can only translate sequences of type SeqRecord, Seq, or str"
        raise ValueError(msg)
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
        msg = 'can only translate sequences of type SeqRecord, Seq, or str'
        raise ValueError(msg)


def translate_ambiguous(seq, gap_char=_GAP, trim_gaps=True):
    if isinstance(seq, SeqRecord):
        seqstr = seq.seq.tostring()
    elif isinstance(seq, Seq):
        seqstr = seq.tostring()
    elif not isinstance(seq, str):
        msg = 'can only enumerate codons of a SeqRecord, Seq, or str'
        raise ValueError(msg)

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


def compute_cigar(reference, record, reference_name=None, new_style=False):

    if reference_name is None:
        reference_name = reference.name

    ncol = len(record)

    # find start, end of record in the ref
    start, end = 0, ncol - 1
    for i in range(ncol):
        if record[i] not in _GAPS:
            start = i
            break
    for i in range(ncol, 0, -1):
        if record[i - 1] not in _GAPS:
            end = i
            break

    MATCH = '=' if new_style else 'M'
    MISMATCH = 'X' if new_style else 'M'

    # compute cigar and edit distance
    cigar = ''
    edit_distance = 0
    mode = ''
    count = 0
    for i in range(start, end):
        ref = reference[i].upper()
        query = record[i].upper()
        if ref in _GAPS and query in _GAPS:
            # if both are gaps, skip
            pass
        elif ref in _GAPS:
            m = 'I'  # insertion
            edit_distance += 1
        elif query in _GAPS:
            m = 'D'  # deletion
            edit_distance += 1
        elif ref == query:
            m = MATCH  # match
        else:
            m = MISMATCH  # mismatch
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
    record.annotations['position'] = start
    record.annotations['length'] = end - start
    record.annotations['reference_name'] = reference_name

    return record


def gapless(seq):
    if not any(char in _GAPS for char in seq):
        return seq
    regexp = re_compile('[{0}]+'.format(''.join(_GAPS)))
    if isinstance(seq, str):
        return regexp.sub('', seq)
    elif isinstance(seq, Seq):
        return Seq(regexp.sub('', str(seq)), seq.alphabet)
    elif isinstance(seq, SeqRecord):
        # TODO: support features and letter_annotations here
        return SeqRecord(
            Seq(regexp.sub('', str(seq.seq)), seq.seq.alphabet),
            id=seq.id,
            name=seq.name,
            dbxrefs=copy(seq.dbxrefs),
            # features=seq.features,
            description=seq.description,
            annotations=copy(seq.annotations)
            # letter_annotations=seq.letter_annotations
            )
    else:
        raise ValueError('seq must have type SeqRecord, Seq, or str')


_cigar_regexp = re_compile(r'([0-9]+)([M=XID])')


def gapful(record, insertions=True):
    p = 0
    modes = 'M=XI' if insertions else 'M=X'
    cigparts = []
    seqparts = ['-' * (record.annotations['position'])]
    for m in _cigar_regexp.finditer(record.annotations['CIGAR']):
        num, mode = int(m.group(1)), m.group(2)
        if mode in modes:
            cigparts.append(m.group(0))
            seqparts.append(str(record.seq[p:(p + num)]))
        elif mode == 'D':
            seqparts.append('-' * num)
        if mode != 'D':
            p += num
    return SeqRecord(
        Seq(''.join(seqparts), record.seq.alphabet),
        id=record.id,
        name=record.name,
        dbxrefs=copy(record.dbxrefs),
        # features = seq.features,
        description=record.description,
        annotations=copy(record.annotations),
        # letter_annotations=record.letter_annotations
        )


def rle_encode(iterable):
    rle = []
    for name, group in groupby(iterable):
        length = 0
        for _ in group:
            length += 1
        rle.append('{0:d}{1:s}'.format(length, name))
    return ''.join(rle)


def clip(record, start, end, insertions=True, span=False):
    if span and record.annotations['position'] > start:
        return None

    bases = []
    poses = []
    modes = []

    pos = record.annotations['position']
    seq = iter(record.seq)
    for m in _cigar_regexp.finditer(record.annotations['CIGAR']):
        num, mode = int(m.group(1)), m.group(2).upper()
        if not insertions and mode == 'I':
            # if we're skipping insertions, consume those bases
            for _ in range(num):
                next(seq)
        elif mode in 'DMIX=':
            for _ in range(num):
                # consume a base even if we're not going to keep it
                base = None if mode == 'D' else next(seq)
                if pos >= start and pos < end:
                    bases.append(base)
                    poses.append(pos)
                    modes.append(mode)
                pos += 0 if mode == 'I' else 1

    if not len(poses):
        return None

    if span and (min(poses) > start or max(poses) < (end - 1)):
        return None

    annotations = copy(record.annotations)
    annotations['CIGAR'] = rle_encode(modes)
    annotations['position'] = min(poses)

    return SeqRecord(
        Seq(''.join(base for base in bases if base), record.seq.alphabet),
        id=record.id,
        name=record.name,
        dbxrefs=copy(record.dbxrefs),
        # features = seq.features,
        description=record.description,
        annotations=annotations,
        # letter_annotations=record.letter_annotations
        )

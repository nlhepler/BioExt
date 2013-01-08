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
    'compute_cigar'
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


def compute_cigar(reference, record, new_style=False):
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
            m = 'I' # insertion
            edit_distance += 1
        elif query in _GAPS:
            m = 'D' # deletion
            edit_distance += 1
        elif ref == query:
            m = MATCH # match
        else:
            m = MISMATCH # mismatch
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

    return record

def _align(aln, ref, seq):
    return aln(ref, seq)

def align_to_refseq(
    reference,
    records,
    score_matrix=None,
    do_codon=True,
    reverse_complement=True,
    expected_identity=None,
    keep_insertions=False,
    quiet=None,
    **kwargs
    ):

    if keep_insertions:
        raise ValueError('keeping insertions is unsupported at this time')

    if score_matrix is None:
        from BioExt.scorematrix import BLOSUM62
        score_matrix = BLOSUM62.load()

    # drop-in compatibility with hy454
    do_codon = kwargs.get('codon', do_codon)
    reverse_complement = kwargs.get('revcomp', reverse_complement)

    # import necessary components
    from os import getenv
    from Bio.Align import MultipleSeqAlignment
    from BioExt.aligner import Aligner
    from sklearn.externals.joblib import Parallel, delayed

    aln = Aligner(
        score_matrix,
        do_codon=do_codon,
        expected_identity=expected_identity
        )

    if reverse_complement:
        def rc(records):
            for record in records:
                yield record
                yield record.reverse_complement()
        records_ = rc(records)
    else:
        records_ = records

    results = Parallel(
        n_jobs=int(getenv('NCPU', -1)),
        verbose=0,
        pre_dispatch='2 * n_jobs'
        )(
            delayed(aln)(reference, record)
            for record in records_
            )

    def finalize(results):
        # grab the best orientation, if we doing that
        if reverse_complement:
            if len(results) % 2 != 0:
                raise ValueError('len(results) should be a multiple of 2')
            results_ = (
                _1 if _1[0] >= _2[0] else _2
                for _1, _2 in zip(results[::2], results[1::2])
                )
        else:
            results_ = results
        # make translation table
        if do_codon:
            tbl = ''.join(
                'N' if c == _GAP else c
                for c in (chr(i) for i in range(256))
                )
        # iterate through the results
        for score, ref, seq in results_:
            # strip insertions
            seq_ = ''.join(s for r, s in zip(ref, seq) if r != _GAP)
            # convert partial-codon deletions to Ns (frameshift deletion)
            if do_codon:
                gap_cdn = 3 * _GAP
                seq_ = ''.join(
                    cdn if (_GAP not in cdn) or cdn == gap_cdn
                    else cdn.translate(tbl)
                    for cdn in by_codon(seq_)
                    )
            # create return value of the proper type (probably overkill)
            if isinstance(seq, SeqRecord):
                rv = SeqRecord(
                    Seq(seq_),
                    id=seq.id,
                    name=seq.name,
                    description=seq.description,
                    features=seq.features
                    )
            elif isinstance(seq, Seq):
                rv = Seq(seq_)
            else:
                rv = seq_
            yield score, rv

    alignment, discard = MultipleSeqAlignment([]), []
    for score, record in finalize(results):
        if aln.expected(score):
            alignment.append(record)
        else:
            discard.append(record)

    return alignment, discard

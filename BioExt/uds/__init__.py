
from __future__ import division, print_function

import os

from ctypes import c_char
from multiprocessing import Array
from operator import itemgetter
from sys import stderr

from Bio.Align import MultipleSeqAlignment

from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord

from BioExt.aligner import Aligner
from BioExt.joblib import Parallel, delayed
from BioExt.misc import compute_cigar, gapful, gapless


__all__ = [
    'align_to_refseq'
    ]


def _rc(record):
    if isinstance(record, str):
        return rc(record)
    elif isinstance(record, Seq):
        return record.reverse_complement()
    elif isinstance(record, SeqRecord):
        return SeqRecord(
            record.seq.reverse_complement(),
            id=record.id,
            name=record.name,
            description=record.description
            )
    else:
        raise ValueError(
            'record must be one of str, Bio.Seq, or Bio.SeqRecord'
            )


# aln, ref, ref_name, and do_revcomp are set by set_globals below
def _align(record):
    records = (record, _rc(record)) if do_revcomp else (record,)
    score, ref_, record = max(
        (aln(ref.value.decode('utf-8'), record) for record in records),
        key=itemgetter(0)
        )
    record_ = compute_cigar(ref_, record, ref_name)
    return score, record_


def _set_globals(*args):
    for key, value in args:
        globals()[key] = value


def _align_par(
        reference,
        records,
        score_matrix,
        do_codon,
        reverse_complement,
        expected_identity,
        discard,
        output,
        quiet=True
        ):

    try:
        n_jobs = int(os.environ.get('NCPU', -1))
    except ValueError:
        n_jobs = -1

    if n_jobs == 0:
        n_jobs = 1

    aln = Aligner(
        score_matrix,
        do_codon=do_codon,
        expected_identity=expected_identity
        )

    if isinstance(reference, str):
        refstr = reference
    elif isinstance(reference, Seq):
        refstr = str(reference)
    elif isinstance(reference, SeqRecord):
        refstr = str(reference.seq)
    else:
        raise ValueError(
            'reference must be one of str, Bio.Seq, Bio.SeqRecord'
            )

    reference_ = Array(c_char, refstr.encode('utf-8'))

    def keep(score, record):
        if aln.expected(score):
            return True
        elif discard:
            discard(record)
        return False

    if quiet:
        def delayed_(i, fn):
            return delayed(fn)
    else:
        def delayed_(i, fn):
            print('\rdispatched: {0:9d} reads'.format(i), end='', file=stderr)
            stderr.flush()
            return delayed(fn)

    rv = output(
        record
        for score, record in Parallel(
            n_jobs=n_jobs,
            verbose=0,
            pre_dispatch='3 * n_jobs',  # triple-buffering
            initializer=_set_globals,
            initargs=[
                ('aln', aln),
                ('ref', reference_),
                ('ref_name', reference.name),
                ('do_revcomp', reverse_complement)
                ]
            ).lazy(
                delayed_(i, _align)(record)
                for i, record in enumerate(records)
                )
        if keep(score, record)
        )

    if not quiet:
        print('', file=stderr)

    return rv


def align_to_refseq(
        reference,
        records,
        score_matrix=None,
        do_codon=True,
        reverse_complement=True,
        expected_identity=None,
        keep_insertions=False,
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

    discards = []

    def discard(record):
        discards.append(record)

    alignment = MultipleSeqAlignment([])

    def output(records):
        for record in records:
            alignment.append(gapful(gapless(record), insertions=False))

    _align_par(
        reference,
        records,
        score_matrix,
        do_codon,
        reverse_complement,
        expected_identity,
        discard,
        output
        )

    return alignment, discards

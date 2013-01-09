
from os import getenv

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt.aligner import Aligner
from BioExt.misc import _GAP, by_codon

from sklearn.externals.joblib import Parallel, delayed


__all__ = [
    'align_to_refseq'
    ]


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

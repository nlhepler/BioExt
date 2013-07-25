#!/usr/bin/env python3

import os
import signal
import sys

from random import shuffle
from re import compile as re_compile

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt.io._SamBamIO import _to_seqrecord
from BioExt.joblib import Parallel, delayed
from BioExt.misc import clip
from BioExt.quiver import NoQVsModelParams, extractFeatures, refineConsensus

import ConsensusCore as cc

import pysam


# global alignment seems to work best ... I don't know why
USE_LOCAL_ALIGNMENT = False
QV_NITER = 20
SHUFFLE_WINDOW_RATIO = 0.2
SHUFFLE_STRIDE_RATIO = 0.5

ndispatched = 0
realph = re_compile(r'[^ACGT]')


def left_pos(seq):
    return seq.annotations['position']


def right_pos(seq):
    for i in range(len(seq) - 1, -1, -1):
        # skip all trailing ambigs
        if seq[i].upper() in 'ACGT':
            return i + 1
    return i


def proper(string):
    return realph.sub('', string.upper())


def _consensus(reads, quiver):

    poa_config = cc.PoaConfig(USE_LOCAL_ALIGNMENT)

    consensus = cc.PoaConsensus.FindConsensus(
        reads,
        poa_config
        ).Sequence()

    if quiver:
        qv_config = cc.QuiverConfig(
            cc.QvModelParams(*NoQVsModelParams.values()),
            cc.ALL_MOVES,
            cc.BandingOptions(4, 5),
            -12.5
            )

        mms = cc.SparseSseQvMultiReadMutationScorer(qv_config, consensus)

        for read in reads:
            mms.AddRead(extractFeatures(read), cc.FORWARD_STRAND)

        consensus, _ = refineConsensus(mms, QV_NITER)

    return consensus


def main(
        bam_handle,
        output_handle,
        max_coverage,
        min_coverage,
        window_size,
        window_stride,
        quiver,
        quiet
        ):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    try:
        n_jobs = int(os.environ.get('NCPU', -1))
    except ValueError:
        n_jobs = -1

    if n_jobs == 0:
        n_jobs = None

    samfile = None
    try:
        pysam.index(bam_handle)
        samfile = pysam.Samfile(bam_handle, 'rb')
        reference = samfile.header['SQ'][0]['SN']
        alignment_length = samfile.header['SQ'][0]['LN']

        windows = []

        for start in range(0, alignment_length, window_stride):
            end = min(alignment_length, start + window_size)
            windows.append((start, end))

        if quiet:
            def delayed_(fn):
                return delayed(fn)
        else:
            def delayed_(fn):
                global ndispatched
                print('\rdispatched: {0:d} jobs'.format(ndispatched), end='', file=sys.stderr)
                sys.stderr.flush()
                ndispatched += 1
                return delayed(fn)

        def read_groups():
            for start, end in windows:
                reads = [
                    str(read_.seq)
                    for read_ in (
                        clip(_to_seqrecord(samfile, read), start, end, span=True)
                        for read in samfile.fetch(reference, start, end)
                        )
                    if read_
                    ]
                if len(reads) < min_coverage:
                    continue
                shuffle(reads)
                yield reads[:max_coverage]

        parts = list(Parallel(
            n_jobs=n_jobs,
            verbose=0,
            pre_dispatch='3 * n_jobs'
            )(
                delayed_(_consensus)(reads, quiver)
                for reads in read_groups()
                )
            )

        if not quiet:
            print('', file=sys.stderr)

        poa_config = cc.PoaConfig(True)

        consensus = cc.PoaConsensus.FindConsensus(
            parts,
            poa_config
            ).Sequence()

        id_ = 'consensus'

        output_handle.write(
            SeqRecord(
                Seq(consensus, single_letter_alphabet),
                id=id_,
                name=id_,
                description=id_
                ).format('fasta')
            )

    finally:
        if samfile:
            samfile.close()

    return 0


if __name__ == '__main__':
    from argparse import ArgumentParser, ArgumentTypeError, FileType

    parser = ArgumentParser(description=(
        "Use PacBio's ConsensusCore PoaConsensus algorithm to generate a consensus"
        ))

    def posint(string):
        try:
            v = int(string)
            assert v > 0
            return v
        except (AssertionError, ValueError):
            raise ArgumentTypeError(
                'must be a positive integer greater than 0'
                )

    parser.add_argument(
        'bam',
        metavar='BAM',
        type=FileType('rb'),
        help='the input sequences, in BAM format'
        )
    parser.add_argument(
        '-c', '--max-coverage',
        metavar='MAX_COVERAGE',
        type=posint,
        default=250,
        help='maximum number of reads to use when estimating consensus'
        )
    parser.add_argument(
        '-m', '--min-coverage',
        metavar='MIN_COVERAGE',
        type=posint,
        default=100,
        help='minimum number of reads to use for estimating consensus'
        )
    parser.add_argument(
        '-o', '--output',
        metavar='FASTA',
        type=FileType('w'),
        default=sys.stdout,
        help='output FASTA file'
        )
    parser.add_argument(
        '-w', '--window-size',
        metavar='WINDOW_SIZE',
        type=posint,
        default=90,
        help='minimum size of sliding windows (use a multiple of 3)'
        )
    parser.add_argument(
        '-s', '--window-stride',
        metavar='WINDOW_STRIDE',
        type=posint,
        default=30,
        help='stride for sliding window consensus estimation'
        )
    parser.add_argument(
        '-Q', '--quiver',
        action='store_true',
        help='use Quiver algorithm to refine consensus'
        )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='do not report the number of sequences dispatched'
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        bam_file = args.bam.name
        args.bam.close()
        retcode = main(
            bam_file,
            args.output,
            args.max_coverage,
            args.min_coverage,
            args.window_size,
            args.window_stride,
            args.quiver,
            args.quiet
            )
    finally:
        if args:
            if args.output != sys.stdout:
                args.output.close()

    sys.exit(retcode)

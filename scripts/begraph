#!/usr/bin/env python3

from __future__ import division, print_function

import sys

from argparse import ArgumentParser, ArgumentTypeError, FileType

import matplotlib as mpl

from Bio import AlignIO, SeqIO


def csvtype(string):
    return string.split(',')


def spectype(string):
    try:
        colspec = []
        if ':' in string:
            xyz = [int(v) for v in string.split(':', 2) if v]
            xyz[0] -= 1
            colspec.extend(v for v in range(*xyz))
        else:
            colspec.extend(int(v) - 1 for v in string.split(',') if v)
        return colspec
    except ValueError:
        msg = (
            'COLSPEC must be a comma-delimited list of 1-indices or '
            'a colon-delimited range specification start:stop[:skip]'
            )
        raise ArgumentTypeError(msg)


def main(
        infile, mode, outfile, format, refidx,
        bins=None, mean_median=True,
        colspec=None, labels=None,
        nseq=None
        ):
    # set these before loading graphing from BioExt, otherwise warnings everywhere
    if format == 'pdf':
        mpl.use('pdf')
    elif format == 'png':
        mpl.use('agg')
    elif format == 'svg':
        mpl.use('svg')
    else:
        raise ValueError('invalid format: {0}'.format(format))

    if mode == 'histogram':
        from BioExt.graphing import graph_readlength_histogram

        lengths = [
            len(seq)
            for i, seq
            in enumerate(SeqIO.parse(infile, 'fasta'))
            if i != refidx
            ]

        graph_readlength_histogram(
            lengths,
            bins=bins,
            mean_median=mean_median,
            filename=outfile,
            format=format
            )
    elif mode == 'seqlogo':
        # this has to be done after setting use() above
        from BioExt.graphing import graph_logo

        msa = AlignIO.read(infile, 'fasta')

        graph_logo(
            msa,
            colspec,
            outfile,
            format=format,
            labels=labels,
            refidx=refidx
            )
    elif mode in ('coverage', 'majority', 'both'):
        # this import has to be done after setting use() above
        from BioExt.graphing import graph_coverage_majority

        msa = AlignIO.read(infile, 'fasta')

        graph_coverage_majority(
            msa,
            mode,
            nseq,
            outfile,
            format=format,
            refidx=refidx
            )
    else:
        raise ValueError('unknown mode {0}'.format(mode))

    return 0


if __name__ == '__main__':
    parser = ArgumentParser(
        description='graph various aspects of sequencing data'
        )

    subparsers = parser.add_subparsers(dest='cmd')

    parser_msa = subparsers.add_parser(
        'msa',
        description='draw a coverage and/or majority graph of a multiple sequence alignment'
        )
    parser_msa.add_argument(
        'input',
        metavar='MSA',
        type=FileType('r'),
        help='aligned FASTA file'
        )
    parser_msa.add_argument(
        'output',
        metavar='OUTPUT',
        help='output file'
        )
    parser_msa.add_argument(
        '-f', '--format',
        choices=('pdf', 'png'), default='pdf',
        help='output format'
        )
    parser_msa.add_argument(
        '-m', '--mode',
        choices=('coverage', 'majority', 'both'), default='both',
        help='graph coverage, majority, or both'
        )
    parser_msa.add_argument(
        '-n', '--nseq',
        type=int, default=None,
        help='number of input sequences'
        )
    parser_msa.add_argument(
        '-r', '--refidx',
        type=int, default=None,
        help='omit reference sequence which is REFIDX-th in the file (1-indexed)'
        )

    parser_hist = subparsers.add_parser(
        'hist',
        description='draw a histogram of read lengths'
        )
    parser_hist.add_argument(
        'input',
        metavar='FASTA',
        type=FileType('r'),
        help='FASTA file'
        )
    parser_hist.add_argument(
        'output',
        metavar='OUTPUT',
        help='OUTPUT file'
        )
    parser_hist.add_argument(
        '-b', '--bins',
        type=int, default=50,
        help='number of bins in the histogram'
        )
    parser_hist.add_argument(
        '-f', '--format',
        choices=('pdf', 'png'), default='pdf',
        help='output format'
        )
    parser_hist.add_argument(
        '-m', '--mean-median',
        dest='mean_median',
        action='store_true',
        help='colorize mean and median in the histogram'
        )
    parser_hist.add_argument(
        '-r', '--refidx',
        type=int, default=None,
        help='omit reference sequence which is REFIDX-th in the file (1-indexed)'
        )

    parser_seqlogo = subparsers.add_parser(
        'seqlogo',
        description='generate a sequence logo from an MSA'
        )
    parser_seqlogo.add_argument(
        'input',
        metavar='MSA',
        type=FileType('r'),
        default=sys.stdin,
        help='aligned FASTA file'
    )
    parser_seqlogo.add_argument(
        'colspec',
        metavar='COLSPEC',
        type=spectype,
        help=(
            'comma-delimited list of 1-indices or '
            'a colon-delimited range specification start:stop[:skip]'
            )
    )
    parser_seqlogo.add_argument(
        'output',
        metavar='OUTPUT',
        help='output file'
    )
    parser_seqlogo.add_argument(
        '-l', '--labels',
        metavar='LABELS',
        type=csvtype,
        help='comma-delimited list of column labels'
    )
    parser_seqlogo.add_argument(
        '-f', '--format',
        choices=('pdf', 'png'), default='pdf',
        help='output format'
        )
    parser_seqlogo.add_argument(
        '-r', '--refidx',
        type=int, default=None,
        help='omit reference sequence which is REFIDX-th in the file (1-indexed)'
    )

    ns = None
    retcode = -1
    try:
        ns = parser.parse_args()
        if ns.cmd == 'msa':
            retcode = main(
                ns.input, ns.mode, ns.output, ns.format, ns.refidx,
                nseq=ns.nseq)
        elif ns.cmd == 'hist':
            retcode = main(
                ns.input, 'histogram', ns.output, ns.format, ns.refidx,
                bins=ns.bins, mean_median=ns.mean_median)
        elif ns.cmd == 'seqlogo':
            retcode = main(
                ns.input, 'seqlogo', ns.output, ns.format, ns.refidx,
                colspec=ns.colspec, labels=ns.labels)
        else:
            parser.print_help()
    finally:
        if ns is not None:
            if 'input' in ns and not ns.input in (None, sys.stdin):
                ns.input.close()
    sys.exit(retcode)

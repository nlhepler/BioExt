#!/usr/bin/env python3

import signal
import argparse

from Bio import SeqIO
from BioExt.io._SamBamIO import _to_seqrecord
from BioExt.misc import gapful

import pysam


def regiontype(string):
    try:
        start, end = string.split(':')
        start = int(start) - 1 if start else None
        end = int(end) if end else None

        if start:
            assert start >= 0
        if end:
            assert end >= 0
        if start and end:
            assert end > start

        return start, end
    except:
        raise argparse.ArgumentTypeError("'{0:s}' is not in the form of START:END".format(string))


def main(bam_file, out_handle, start=None, end=None):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    try:
        # Index bam file in order to samfile.fetch
        pysam.index(bam_file)
        samfile = pysam.Samfile(bam_file, 'rb')
        length = samfile.header['SQ'][0]['LN']
        fetch_args = []

        if start is not None or end is not None:
            fetch_args.extend([samfile.getrname(0), start, end])

        SeqIO.write(
            (
                (seq + ('-' * (length - len(seq))))[start:end]
                for seq in (
                    gapful(_to_seqrecord(samfile, record), insertions=False)
                    for record in samfile.fetch(*fetch_args)
                    )
                ),
            out_handle,
            'fasta'
            )
    finally:
        if samfile is not None:
            samfile.close()

    return 0


if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(
        description='convert a BAM file to a MSA, removing insertions'
        )

    parser.add_argument(
        '-r', '--region',
        metavar='REGION',
        type=regiontype,
        default=(None, None),
        help='only include sequences in a certain REGION'
        )
    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('rb'),
        help='input BAM file'
        )
    parser.add_argument(
        'output',
        default=sys.stdout,
        metavar='OUTPUT',
        type=argparse.FileType('w'),
        help='output FASTA MSA file'
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        bam_file = args.input.name
        args.input.close()
        retcode = main(bam_file, args.output, *args.region)
    finally:
        if args is not None:
            if args.output != sys.stdout:
                args.output.close()

    sys.exit(retcode)

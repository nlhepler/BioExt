#!/usr/bin/env python3

import signal
import argparse

from itertools import chain

from Bio import SeqIO

from BioExt.args import add_reference
from BioExt.io import BamIO


def main(in_handle, bam_file, reference=None):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    records = SeqIO.parse(in_handle, 'fasta')

    if reference is None:
        try:
            reference = next(records)
            records = chain.from_iterable(((reference,), records))
        except StopIteration:
            raise ValueError('provide a non-empty MSA')

    BamIO.write(
        records,
        bam_file,
        reference
        )

    return 0


if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(
        description='convert a BAM file to a MSA, removing insertions'
        )

    add_reference(parser, '-r', '--reference')
    parser.add_argument(
        'input',
        default=sys.stdout,
        metavar='INPUT',
        type=argparse.FileType('r'),
        help='input FASTA MSA file'
        )
    parser.add_argument(
        'output',
        metavar='OUTPUT',
        type=argparse.FileType('wb'),
        help='output BAM file'
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        bam_file = args.output.name
        args.output.close()
        retcode = main(args.input, bam_file, args.reference)
    finally:
        if args is not None:
            if args.input != sys.stdout:
                args.input.close()

    sys.exit(retcode)

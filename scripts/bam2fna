#!/usr/bin/env python3

import signal

from Bio import SeqIO
from BioExt.io import BamIO


def main(bam_file, out_handle):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    SeqIO.write(BamIO.parse(bam_file), out_handle, 'fasta')

    return 0


if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description='convert a BAM file to a FASTA file'
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
        help='output FASTA file'
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        bam_file = args.input.name
        args.input.close()
        retcode = main(bam_file, args.output)
    finally:
        if args is not None:
            if args.output != sys.stdout:
                args.output.close()

    sys.exit(retcode)

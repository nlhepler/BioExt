#!/usr/bin/env python3

import signal
import argparse

from Bio import SeqIO

from BioExt.collections import Counter


def transform(record):
    return hash(str(record.seq).upper())


def newrecord(counts, record, sep='_N:'):
    r = transform(record)
    count = counts[r]
    if count == 0:
        return None
    s = '{0:s}{1:d}'.format(sep, count)
    record.id += s
    record.description += s
    del counts[r]
    return record


def main(in_handle, out_handle, sep='_N:'):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    counts = Counter(transform(record) for record in SeqIO.parse(in_handle, 'fasta'))

    in_handle.seek(0)

    SeqIO.write(
        (
            record_
            for record_ in (
                newrecord(counts, record, sep)
                for record in SeqIO.parse(in_handle, 'fasta')
                )
            if record_),
        out_handle,
        'fasta'
        )

    return 0


if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(
        description='convert a BAM file to a MSA, removing insertions'
        )

    parser.add_argument(
        '-s', '--separator',
        metavar='SEP',
        type=str,
        default='_N:',
        help="sequence count will be appended as '(SEP)COUNT'"
        )
    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('r'),
        help='input FASTA MSA file'
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
        retcode = main(args.input, args.output, args.separator)
    finally:
        if args is not None:
            if args.input != sys.stdin:
                args.input.close()
            if args.output != sys.stdout:
                args.output.close()

    sys.exit(retcode)

#!/usr/bin/env python3

import signal, sys

from argparse import ArgumentParser, FileType

from Bio import SeqIO

from BioExt.misc import translate


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    parser = ArgumentParser(description='translate a FASTA nucleotide file')
    parser.add_argument('-f', '--frame', type=int, choices=range(3), default=0)
    parser.add_argument('input',  nargs='?', type=FileType('r'), default=sys.stdin)
    parser.add_argument('output', nargs='?', type=FileType('w'), default=sys.stdout)

    ns = parser.parse_args(args)

    SeqIO.write((translate(s[ns.frame:]) for s in SeqIO.parse(ns.input, 'fasta')), ns.output, 'fasta')

    return 0


if __name__ == '__main__':
    sys.exit(main())

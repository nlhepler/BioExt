#!/usr/bin/env python3

import signal
import argparse

from itertools import chain

from Bio.Seq import Seq

from BioExt.args import add_reference
from BioExt.misc import gapful
from BioExt.io import BamIO


GAP = '-'


def trim_edges(record, reflen):
    record_ = gapful(record, insertions=False)

    # remove these so BamIO.write regenerates the CIGAR string
    for annotation in ('CIGAR', 'position'):
        if annotation in record_.annotations:
            del record_.annotations[annotation]

    # trim from the left first
    seq = []
    for i, ltr in enumerate(record_.seq):
        if i % 3 == 0 and ltr != GAP:
            # now trim from the right
            rest = record_.seq[i:]
            j = len(rest) - 1
            for ltr_ in reversed(rest):
                if j % 3 == 2 and ltr_ != GAP:
                    # increment because python thinks 0-indexing and [,) is intuitive...
                    upr = j + 1
                    seq.append(str(rest[:upr]))
                    seq.append(GAP * (len(rest) - upr))
                    break
                j -= 1
            break
        else:
            seq.append(GAP)

    # must be reference length
    seq.append(GAP * (reflen - len(record_)))

    record_.seq = Seq(''.join(seq), record_.seq.alphabet)

    return record_


def main(in_file, out_file, reference=None):

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass

    records = BamIO.parse(in_file)

    if reference is None:
        try:
            reference = next(records)
            records = chain.from_iterable(((reference,), records))
        except StopIteration:
            raise ValueError('provide a non-empty MSA')

    reflen = len(reference)

    BamIO.write(
        (trim_edges(record, reflen) for record in records),
        out_file,
        reference
        )

    return 0


if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(
        description='clip non-codon aligned edges in a BAM file'
        )

    add_reference(parser, '-r', '--reference')
    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('rb'),
        help='input BAM file'
        )
    parser.add_argument(
        'output',
        metavar='OUTPUT',
        type=argparse.FileType('wb'),
        help='output BAM file'
        )

    args = parser.parse_args()

    in_file = args.input.name
    out_file = args.output.name
    args.input.close()
    args.output.close()

    retcode = main(in_file, out_file, args.reference)

    sys.exit(retcode)

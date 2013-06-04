#!/usr/bin/env python3

import sys

from re import compile as re_compile

from Bio.Seq import Seq

from BioExt.io import BamIO

from pysam import Samfile


recig = re_compile(r'(?:(\d+)([MIDX=]))')


def clip(reads, reflen):
    for read in reads:
        col = read.annotations['position']
        cig = []
        seq = []
        lwr = 0
        for m in recig.finditer(read.annotations['CIGAR']):
            op = m.group(2)
            nop = int(m.group(1))
            upr = lwr + nop
            if (op == 'I' and col < 1) or (op == 'I' and col >= reflen):
                pass
            else:
                cig.append('{0:d}{1:s}'.format(nop, op))
                seq.append(str(read.seq[lwr:upr]))
            if op != 'I':
                col += nop
            lwr = upr
        read.seq = Seq(''.join(seq))
        read.annotations['CIGAR'] = ''.join(cig)
        yield read


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    f = Samfile(args[0])
    header = f.header
    f.close()

    reflen = header['SQ'][0]['LN']

    BamIO.write(clip(BamIO.parse(args[0]), reflen), args[1], header=header)

    return 0

if __name__ == '__main__':
    sys.exit(main())

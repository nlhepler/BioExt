#!/usr/bin/env python3

import sys

from math import log
from operator import itemgetter

import numpy as np

from Bio import SeqIO

from BioExt.rateclass import RateClass, p_bg


def idx2nuc(idx):
    if idx == 0:
        return 'A'
    elif idx == 1:
        return 'C'
    elif idx == 2:
        return 'G'
    elif idx == 3:
        return 'T'
    else:
        raise ValueError('unknown idx {0}'.format(idx))


def nuc2idx(nuc):
    if nuc == 'A':
        return 0
    elif nuc == 'C':
        return 1
    elif nuc == 'G':
        return 2
    elif nuc == 'T':
        return 3
    else:
        raise ValueError('unknown nuc {0}'.format(nuc))


def weighted_harmonic_mean(params):
    return sum(p[0] for p in params) / sum(p[0] / p[1] for p in params)


def weighted_mean(params):
    return sum(p[0] * p[1] for p in params) / sum(p[0] for p in params)


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    mat = None

    with open(args[0]) as handle:
        msa = SeqIO.parse(handle, 'fasta')

        for nr, read in enumerate(msa, start=1):
            print('\rreads:    {0:7d}'.format(nr + 1), end='', file=sys.stderr)
            if mat is None:
                mat = np.zeros((len(read), 4), dtype=int)
            if len(read) != mat.shape[0]:
                raise RuntimeError('Invalid read length: input is not an MSA')
            for i, nuc in enumerate(read):
                try:
                    j = nuc2idx(nuc)
                except ValueError:
                    continue
                mat[i, j] += 1
            if nr % 100 == 0:
                sys.stderr.flush()

    sys.stderr.flush()

    nrow, ncol = mat.shape
    data = []
    css = []

    for i in range(nrow):
        cov = mat[i, :].sum()
        max_idx, maj = max(enumerate(mat[i, :]), key=itemgetter(1))
        css.append(idx2nuc(max_idx))
        for j in range(4):
            k = mat[i, j]
            if j == max_idx or k == 0:
                continue
            print('{0}\t{1}'.format(cov, cov - k))
            data.append((cov, cov - k))

    return 0

    css_ = ''.join(css)
    lg_L, aicc, params = RateClass(data, 3)()

    bg = params[0][1]
    # bg = weighted_mean(params)
    # bg = weighted_harmonic_mean(params)

    print('', file=sys.stderr)
    print('consensus: ', css_, file=sys.stderr)
    print('lg_L:      {0: .7f}'.format(lg_L), file=sys.stderr)
    print('aicc:      {0: .7f}'.format(aicc), file=sys.stderr)
    print('background:{0: .7f}'.format(bg), file=sys.stderr)
    print('rates:    [', ', '.join('{0:.7f}'.format(p[1]) for p in params), ']', file=sys.stderr)
    print('weights:  [', ', '.join('{0:.7f}'.format(p[0]) for p in params), ']', file=sys.stderr)

    lg_bg = log(bg)
    lg_invbg = log(1.0 - bg)

    for i in range(nrow):
        cov = mat[i, :].sum()
        max_idx, maj = max(enumerate(mat[i, :]), key=itemgetter(1))
        css = idx2nuc(max_idx)
        for j in range(4):
            if j == max_idx:
                continue
            k = mat[i, j]
            p = p_bg(lg_bg, lg_invbg, cov, k)
            if p < 0.01:
                var = idx2nuc(j)
                print('{0}\t{1}\t{2}:{3}\t{4}:{5}\t{6}'.format(
                    i + 1, cov, css, maj, var, k, p))

    return 0


if __name__ == '__main__':
    sys.exit(main())

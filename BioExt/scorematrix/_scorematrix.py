
from __future__ import division, print_function

from itertools import chain
from math import log
from re import compile as re_compile
from warnings import warn

import numpy as np

from BioExt.collections import OrderedDict


__all__ = [
    'DNA65',
    'DNA70',
    'DNA80',
    'DNA88',
    'DNA95',
    'parse_scorematrix',
    'DNAScoreMatrix',
    'DNAExpIdScoreMatrix',
    'ProteinScoreMatrix'
]


pletters = 'ARNDCQEGHILKMFPSTWYVBZX*'
dletters = 'ACGT'


def parse_scorematrix(smpath):
    with open(smpath) as fh:
        ws = re_compile(r'\s+')
        comment = re_compile(r'^\s*#')
        S = fh.read().split('\n')
        T = [s.strip() for s in S if not comment.match(s)]
        U = [ws.sub(' ', t).split(' ') for t in T if len(t) > 0]
        V = [u[1:] for u in U[1:]]
        W = [[int(w) for w in v] for v in V]
        lettersX = ''.join(U[0]).upper()
        lettersY = ''.join([u[0] for u in U[1:]]).upper()
        if len(lettersX) >= 20:
            letters = pletters
            klass = ProteinScoreMatrix
        else:
            letters = dletters
            klass = DNAScoreMatrix
        if not set(letters).issubset(set(lettersX + lettersY)):
            msg = "scoring matrix '%s' is insufficiently descriptive" % smpath
            raise RuntimeError(msg)
        if lettersX != lettersY or lettersX[:len(letters)] != letters:
            cols = [lettersX.index(l) for l in letters]
            rows = [lettersY.index(l) for l in letters]
            return klass([[W[i][j] for j in cols] for i in rows])
        else:
            return klass(W)


def _normalize(matrix):
    max_ = matrix.max()
    min_ = matrix.min()
    if max_ == min_:
        matrix[:] = 0.5
    else:
        matrix -= min_
        matrix /= max_ - min_
    return matrix


class ScoreMatrix(object):

    def __init__(self, matrix, letters):
        self.__matrix = matrix
        self.__letters = ''.join(letters) # make sure it's a string
        self.__format = '%% %ds' % max(len(str(d)) for d in chain(*matrix))

    @property
    def letters(self):
        return self.__letters

    def tondarray(self):
        N = len(self.__letters)
        mat = np.zeros((N, N), dtype=float)
        for i in range(N):
            for j in range(N):
                mat[i, j] = self.__matrix[i][j]
        return mat

    def __getitem__(self, key):
        if not isinstance(key, tuple) or len(key) != 2:
            raise KeyError("indexing into a ScoreMatrix requires 2D indices")
        return self.__matrix[key[0]][key[1]]

    def _tohyphy(self, names):
        if not isinstance(names, tuple) or len(names) != 2:
            raise ValueError("must provide (letter, matrix) variable names in that order")
        lname, mname = names
        mstr = '{\n {%s}\n}' % '}\n {'.join(
            ', '.join(self.__format % str(s) for s in row) for row in self.__matrix
        )
        return '%s = "%s";\n%s = %s;' % (lname, self.__letters[:20], mname, mstr)

    def __call__(self, ref, seq, mismatch=0):
        return ScoreMatrix.score(self, ref, seq, mismatch)

    # so we don't have to worry about calling load() on DNA* matrices
    def load(self):
        return self

    def score(self, ref, seq, mismatch=0):
        if not len(ref) == len(seq):
            raise ValueError("reference and query are not aligned")
        alph = dict((v, i) for i, v in enumerate(self.__letters))
        score = 0
        for i, r in enumerate(ref):
            q = seq[i]
            if '-' in (r, q):
                score -= mismatch
                continue
            try:
                score += self.__matrix[alph[r]][alph[q]]
            except KeyError:
                l = r if r not in alph else q
                warn("unknown letter '%s', ignoring" % l)
        return score

    def freqs(self, bias=None):
        from scipy.linalg import eig

        N = len(self.__letters)
        unif = 1 / N

        if bias is None:
            bias = 10 ** int(np.log10(unif) // 1)
        elif bias < 0 or bias > unif:
            raise ValueError(
                'bias must be a probability between [0, {0}]'.format(unif)
                )

        matrix = self.tondarray()
        # convert to rate matrix
        matrix = np.subtract(matrix.T, matrix[np.eye(N, dtype=bool)]).T
        matrix = np.divide(matrix.T, matrix.sum(axis=1)).T
        # compute left eigenvectors
        w, vl, _ = eig(matrix, left=True)
        # extract the frequencies
        freqs = vl.T[np.abs(w - 1) < 1e-7].reshape((N,))

        # handle bias
        if any(freqs < bias):
            freqs -= freqs.min()
            bias_ = bias * freqs.sum() / (1.0 - N * bias)
            freqs += bias_

        # normalize
        freqs /= freqs.sum()

        if freqs.dtype == complex:
            freqs = (f.real for f in freqs)

        return OrderedDict(zip(self.__letters, freqs))


class DNAScoreMatrix(ScoreMatrix):

    def __init__(self, matrix, letters=dletters):
        super(DNAScoreMatrix, self).__init__(matrix, letters)


class DNAExpIdScoreMatrix(DNAScoreMatrix):
    UNIFORM_FREQS = {
        'A': 0.25,
        'C': 0.25,
        'G': 0.25,
        'T': 0.25
    }

    def __init__(self, expected_identity, freqs=None):
        if freqs is None:
            freqs = DNAExpIdScoreMatrix.UNIFORM_FREQS
        if not set(freqs.keys()).issubset(set(dletters)):
            msg = "frequencies provided to not address each of '%s'" % dletters
            raise ValueError(msg)
        lam = 1. / min(list(freqs.values()))
        N = len(dletters)
        pab = expected_identity / N
        pnab = (1. - expected_identity) / (N ** 2 - N)
        matrix = [[0 for _ in range(N)] for _ in range(N)]
        for i, l in enumerate(dletters):
            for j, k in enumerate(dletters):
                if l == k:
                    matrix[i][j] = int(round(lam * log(pab  / (freqs[l] * freqs[k]))))
                else:
                    matrix[i][j] = int(round(lam * log(pnab / (freqs[l] * freqs[k]))))
        super(DNAExpIdScoreMatrix, self).__init__(matrix, dletters)


class ProteinScoreMatrix(ScoreMatrix):

    def __init__(self, matrix, letters=pletters):
        super(ProteinScoreMatrix, self).__init__(matrix, letters)


DNA65 = DNAExpIdScoreMatrix(0.65)
DNA70 = DNAExpIdScoreMatrix(0.70)
DNA80 = DNAExpIdScoreMatrix(0.80)
DNA88 = DNAExpIdScoreMatrix(0.88)
DNA95 = DNAExpIdScoreMatrix(0.95)

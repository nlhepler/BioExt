
from __future__ import division, print_function

from itertools import chain
from math import log
from re import compile as re_compile
from warnings import warn

import numpy as np

from BioExt.collections import OrderedDict
from BioExt.optimize import minimize_bfgs


__all__ = [
    'FrequenciesError',
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


class FrequenciesError(RuntimeError):
    pass


def parse_scorematrix(name, smpath):
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
            return klass(name, [[W[i][j] for j in cols] for i in rows])
        else:
            return klass(name, W)


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

    def __init__(self, name, matrix, letters):
        self.__name = name
        self.__matrix = matrix
        self.__letters = ''.join(letters)  # make sure it's a string
        self.__format = '%% %ds' % max(len(str(d)) for d in chain(*matrix))

    def __repr__(self):
        return self.__name

    def __str__(self):
        return self.__name

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

        def find(v):
            if isinstance(v, str):
                i = self.__letters.find(v)
                if i < 0:
                    raise ValueError("'{0}' is not in the alphabet".format(v))
                return i
            elif isinstance(v, int):
                return v
            else:
                raise ValueError('key must be an integer index or alphabet character')

        a = find(key[0])
        b = find(key[1])
        return self.__matrix[a][b]

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

    def freqs(self, tol=1e-7):

        def M(lam, matrix):
            return np.exp(matrix * lam)

        def Y(lam, matrix):
            return np.linalg.inv(M(lam, matrix))

        def f(lam, matrix):
            return (Y(lam, matrix).sum() - 1) ** 2

        if isinstance(self, DNAScoreMatrix):
            valid_letters = set(dletters[:4])
        elif isinstance(self, ProteinScoreMatrix):
            valid_letters = set(pletters[:20])
        else:
            valid_letters = set(self.letters)

        indices = np.array([letter in valid_letters for letter in self.letters], dtype=bool)
        matrix = self.tondarray()[indices, :][:, indices]

        res = minimize_bfgs(f, 1.0 / 3, args=(matrix,), tol=tol)
        msg = 'unable to determine frequency distribution'

        if not res.success:
            raise FrequenciesError(msg)

        lam = res.x
        freqs = np.zeros((len(self.letters),), dtype=float)
        freqs[indices] = Y(lam, matrix).sum(axis=0)
        freqs /= freqs.sum()

        if not (freqs >= 0).all():
            raise FrequenciesError(msg)

        return OrderedDict(zip(self.__letters, freqs))

    @property
    def __valid_scores(self):
        if isinstance(self, DNAScoreMatrix):
            valid_letters = set(dletters[:4])
        elif isinstance(self, ProteinScoreMatrix):
            valid_letters = set(pletters[:20])
        else:
            return chain.from_iterable((s for s in r) for r in self.__matrix)
        idxs = [i for i, l in enumerate(self.letters) if l in valid_letters]
        return chain.from_iterable((self.__matrix[i][j] for j in idxs) for i in idxs)

    def max(self):
        return max(self.__valid_scores)

    def min(self):
        return min(self.__valid_scores)


class DNAScoreMatrix(ScoreMatrix):

    def __init__(self, name, matrix, letters=dletters):
        super(DNAScoreMatrix, self).__init__(name, matrix, letters)


class DNAExpIdScoreMatrix(DNAScoreMatrix):
    UNIFORM_FREQS = {
        'A': 0.25,
        'C': 0.25,
        'G': 0.25,
        'T': 0.25
    }

    def __init__(self, expected_identity, freqs=None):
        name = 'DNA{0:d}'.format(round(100 * expected_identity))
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
                    matrix[i][j] = int(round(lam * log(pab / (freqs[l] * freqs[k]))))
                else:
                    matrix[i][j] = int(round(lam * log(pnab / (freqs[l] * freqs[k]))))
        super(DNAExpIdScoreMatrix, self).__init__(name, matrix, dletters)


class ProteinScoreMatrix(ScoreMatrix):

    def __init__(self, name, matrix, letters=pletters):
        super(ProteinScoreMatrix, self).__init__(name, matrix, letters)


DNA65 = DNAExpIdScoreMatrix(0.65)
DNA70 = DNAExpIdScoreMatrix(0.70)
DNA80 = DNAExpIdScoreMatrix(0.80)
DNA88 = DNAExpIdScoreMatrix(0.88)
DNA95 = DNAExpIdScoreMatrix(0.95)

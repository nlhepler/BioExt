
from __future__ import division, print_function

from itertools import chain
from math import log
from re import compile as re_compile
from warnings import warn

from numpy import mean


__all__ = [
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


class ScoreMatrix(object):

    def __init__(self, matrix, letters):
        self.__matrix = matrix
        self.__letters = letters
        self.__format = '%% %ds' % max(len(str(d)) for d in chain(*matrix))

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

    def score(self, ref, seq, mismatch=0):
        if not len(ref) == len(seq):
            raise ValueError("reference and query are not aligned")
        alph = dict((v, i) for i, v in enumerate(self.__letters))
        score = 0
        for i, r in enumerate(ref):
            q = seq[i]
            if '-' in (r, q):
                score -= mismatch
            try:
                score += self.__matrix[alph[r]][alph[q]]
            except KeyError:
                l = r if r not in alph else q
                warn("unknown letter '%s', ignoring" % l)
        return score


class DNAScoreMatrix(ScoreMatrix):

    def __init__(self, matrix, letters=dletters):
        super(DNAScoreMatrix, self).__init__(matrix, letters)


class DNAExpIdScoreMatrix(ScoreMatrix):

    def __init__(self, expected_identity, freqs):
        if not set(freqs.keys()).issubset(set(dletters)):
            msg = "frequencies provided to not address each of '%s'" % dletters
            raise ValueError(msg)
        lam = 1. / mean(list(freqs.values()))
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

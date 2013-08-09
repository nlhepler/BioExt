
from __future__ import division, print_function

from random import randint, random

from scipy.stats import norm

from BioExt.misc import homocount


__all__ = [
    'pyro_errors',
    'errorize'
    ]


# take from:
# Characteristics of 454 pyrosequencing data—enabling realistic simulation with flowsim
# Susanne Balzer, Ketil Malde, Anders Lanzén, Animesh Sharma, and Inge Jonassen
# Bioinformatics Vol. 26 ECCB 2010, pages i420–i425
# doi:10.1093/bioinformatics/btq365
_pyro_errors = [
    (0.1230, 0.0737),
    (1.0193, 0.1227),
    (2.0006, 0.1585),
    (2.9934, 0.2188),
    (3.9962, 0.3168),
    (4.9550, 0.3863)
    ]

def pyro_errors(x):
    if x < 0:
        raise ValueError('x must be >= 0')
    elif x < 6:
        return _pyro_errors[x]
    else:
        return (x, 0.03494 + x * 0.0685)


def seqgen(alph, randlen, sequence):
    l0 = len(alph) - 1
    yield alph[randint(0, l0)] * randlen(0)
    for x, n in homocount(sequence):
        yield x * randlen(n)
        yield alph[randint(0, l0)] * randlen(0)


def errorize(sequence, randlen=None):
    if not isinstance(sequence, str):
        raise ValueError('sequence must be of type str')
    if randlen is None:
        randlen = lambda l: round(norm.ppf(random(), *pyro_errors(l)))
    alph = tuple(set(sequence))
    return ''.join(seqgen(alph, randlen, sequence))

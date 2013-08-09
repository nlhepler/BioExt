
from __future__ import division, print_function

from BioExt.scorematrices._scorematrix import *
from BioExt.scorematrices._lazyscorematrix import LazyScoreMatrix


__all__ = ['LazyScoreMatrix']
__all__ += _scorematrix.__all__


def _init():
    from sys import modules
    from BioExt.scorematrices._factory import _smfactory
    for k, v in _smfactory().items():
        __all__.append(k)
        vars(modules[__name__])[k] = v

_init()

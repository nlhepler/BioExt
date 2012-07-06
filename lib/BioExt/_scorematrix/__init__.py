
from ._scorematrix import *


__all__ = []
__all__ += _scorematrix.__all__


def _init():
    from sys import modules
    from ._factory import _smfactory
    for k, v in _smfactory().items():
        __all__.append(k)
        vars(modules[__name__])[k] = v

_init()

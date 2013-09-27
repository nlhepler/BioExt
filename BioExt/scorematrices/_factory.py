
from __future__ import division, print_function

from glob import iglob
from os.path import abspath, basename, join, split
from re import search as re_search


from BioExt.scorematrices._lazyscorematrix import _lazyscorematrix


__all__ = []


_smdir = join(
    split(
        split(
            abspath(__file__)
        )[0]  # _scorematrix/
    )[0],  # _scorematrix/../
    'data',
    'scorematrices'
)


def _smfactory(smdir=_smdir, smfmt='%s.txt'):
    matrices = {}
    for smpath in iglob(join(smdir, smfmt % '*')):
        m = re_search(smfmt % '(.+)', smpath)
        if m:
            name = basename(m.group(1))
            matrices[name] = _lazyscorematrix(name, smdir, basename(smpath))
    return matrices

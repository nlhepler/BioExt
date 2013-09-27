
from __future__ import division, print_function

from os.path import exists, join

from BioExt.scorematrices._scorematrix import parse_scorematrix


__all__ = []


def _lazyscorematrix(name, smdir, smfile):
    path = join(smdir, smfile)
    if exists(path):
        return LazyScoreMatrix(name, path)
    else:
        raise RuntimeError('cannot load Scoring Matrix, data not found!')


class LazyScoreMatrix(object):

    def __init__(self, name, path):
        self.__name = name
        self.__path = path

    def __repr__(self):
        return self.__name

    def __str__(self):
        return self.__name

    def load(self):
        if exists(self.__path):
            return parse_scorematrix(self.__name, self.__path)
        else:
            msg = "cannot load ScoreMatrix '%s', file missing!" % self.__path
            raise RuntimeError(msg)

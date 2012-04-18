
from __future__ import division, print_function

from os.path import join, exists
from os import path

from Bio import SeqIO


__all__ = []


def _lazyseq(seqdir, seqfile):
    seqpath = join(seqdir, seqfile)
    if exists(seqpath):
        return Lazyseq(seqpath)
    else:
        raise RuntimeError('cannot load reference, gene not found!')


class Lazyseq(object):

    def __init__(self, seqpath):
        self._seqpath = seqpath

    def load(self):
        with open(self._seqpath) as fh:
            record = SeqIO.read(fh, 'fasta')
        return record

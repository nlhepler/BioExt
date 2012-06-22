
from __future__ import division, print_function

from os.path import join, exists, splitext
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
        _, ext = splitext(self._seqpath)
        if ext in ('.gb',):
            filetype = 'genbank'
        elif ext in ('.fa', '.faa', '.fna'):
            filetype = 'fasta'
        else:
            msg = "reference has an unknown file type extension '%s'" % ext
            raise ValueError(msg)
        with open(self._seqpath) as fh:
            record = SeqIO.read(fh, filetype)
        return record

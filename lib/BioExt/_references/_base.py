
from os import path

from Bio import SeqIO


__all__ = []


class Lazyseq(object):

    def __init__(self, seqpath):
        self._seqpath = seqpath

    def load(self):
        with open(self._seqpath) as fh:
            record = SeqIO.read(fh, 'fasta')
        return record


class Refseq(object):

    @staticmethod
    def _lazyseq(rdir, rfile):
        fname = path.join(
            path.split(
                path.split(
                    path.abspath(__file__)
                )[0] # _references/
            )[0], # _references/../
            'data',
            'references',
            rdir,
            rfile
        )
        if path.exists(fname):
            return Lazyseq(fname)
        else:
            raise RuntimeError('cannot load reference, gene not found!')


from ._base import Refseq


__all__ = []


class _nl4_3(Refseq):

    @staticmethod
    def _lazyseq(gene):
        return Refseq._lazyseq('nl4-3', 'NL4_3_%s.fa' % gene)

    @property
    def prrt(self):
        return _nl4_3._lazyseq('prrt')

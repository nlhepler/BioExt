
from ._base import Refseq


__all__ = []


class _hxb2(Refseq):

    @staticmethod
    def _lazyseq(gene):
        return Refseq._lazyseq('hxb2', 'HXB2_%s.fa' % gene)

    @property
    def env(self):
        return _hxb2._lazyseq('env')

    @property
    def gag(self):
        return _hxb2._lazyseq('gag')

    @property
    def int(self):
        return _hxb2._lazyseq('int')

    @property
    def nef(self):
        return _hxb2._lazyseq('nef')

    @property
    def pol(self):
        return _hxb2._lazyseq('pol')

    @property
    def pr(self):
        return _hxb2._lazyseq('pr')

    @property
    def prrt(self):
        return _hxb2._lazyseq('prrt')

    @property
    def rev(self):
        return _hxb2._lazyseq('rev')

    @property
    def rt(self):
        return _hxb2._lazyseq('rt')

    @property
    def tat(self):
        return _hxb2._lazyseq('tat')

    @property
    def vif(self):
        return _hxb2._lazyseq('vif')

    @property
    def vpr(self):
        return _hxb2._lazyseq('vpr')

    @property
    def vpu(self):
        return _hxb2._lazyseq('vpu')

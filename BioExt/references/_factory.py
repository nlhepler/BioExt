
from __future__ import division, print_function

from collections import namedtuple
from glob import iglob
from os.path import (abspath, basename,
    join, split, splitext)
from re import (compile as re_compile,
    search as re_search, I as re_I)


from BioExt.references._lazyseq import _lazyseq


__all__ = []


_installrefdirs = []


_refdir = join(
    split(
        split(
            abspath(__file__)
        )[0] # _references/
    )[0], # _references/../
    'data',
    'references'
)


def _reffactory(seqdir, seqfmt, name=None, genes=None):
    if name is None:
        name = re_compile(r'[^0-9A-Z]', re_I).sub('_', basename(seqdir))
    datadict = {}
    if genes is not None:
        for gene in genes:
            datadict[gene] = _lazyseq(seqdir, seqfmt % gene)
    else:
        genes = []
        for seqpath in iglob(join(seqdir, seqfmt % '*')):
            m = re_search(seqfmt % '(.+)', seqpath)
            if m:
                gene = m.group(1)
                genes.append(gene)
                datadict[gene] = _lazyseq(seqdir, basename(seqpath))
    # if the seqdir has the refdir in it, then add it to the list
    # of default-installed reference sequence directories
    if _refdir in seqdir:
        globber = '*' + splitext(seqfmt)[1]
        _installrefdirs.append(
            join(
                'data',
                'references',
                basename(seqdir),
                globber
            )
        )
    return namedtuple(name, genes)(**datadict)


from os.path import join

from ._factory import _refdir, _reffactory


__all__ = ['hxb2', 'nl4_3']


hxb2 = _reffactory(join(_refdir, 'hxb2'), 'HXB2_%s.fa')
nl4_3 = _reffactory(join(_refdir, 'nl4-3'), 'NL4_3_%s.fa')

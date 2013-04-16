
from __future__ import division, print_function

from BioExt.io import _SamBamIO


__all__ = [
    'parse',
    'write'
    ]


def parse(path):
    for record in _SamBamIO._parse('r', path):
        yield record


def write(records, path, reference=None, new_style=False):
    return _SamBamIO._write('w', path, reference, new_style)

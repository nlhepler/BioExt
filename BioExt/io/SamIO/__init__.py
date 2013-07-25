
from __future__ import division, print_function

from BioExt.io import _SamBamIO


__all__ = [
    'parse',
    'write'
    ]


def parse(path):
    return _SamBamIO._parse('r', path)


def write(records, path, reference=None, new_style=False, header=None):
    return _SamBamIO._write('w', path, reference, new_style, header)

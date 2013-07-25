
from os import close, remove
from os.path import exists, getsize
from shutil import move
from tempfile import mkstemp

from BioExt.io import _SamBamIO

from pysam import sort as pysam_sort


__all__ = [
    'parse',
    'write',
    'sort'
    ]


def parse(path, index=True):
    return _SamBamIO._parse('rb', path, index)


def write(records, path, reference=None, new_style=False, header=None):
    return _SamBamIO._write('wb', records, path, reference, new_style, header)


def sort(path):
    try:
        fd, tmp_path = mkstemp()

        close(fd)

        if exists(path) and getsize(path):
            pysam_sort(path, tmp_path)
            tmp_path += '.bam'  # sort adds the .bam suffix automatically
            move(tmp_path, path)
    finally:
        if exists(tmp_path):
            remove(tmp_path)


from os import close, remove
from os.path import exists
from shutil import move
from tempfile import mkstemp

from pysam import view as samtools_view, sort as samtools_sort

from BioExt.io.SamIO import parse as sam_parse, write as sam_write


__all__ = [
    'parse',
    'write'
    ]


def parse(handle):
    try:
        fd, path = mkstemp(); close(fd)

        with open(path, 'w') as fh:
            for record in samtools_view(handle.name):
                fh.write(record)

        with open(path, 'r') as fh:
            for record in sam_parse(fh):
                yield record
    finally:
        if exists(path):
            remove(path)


def write(records, handle, reference, new_style=False):
    handle_mode = handle.mode.lower()
    if 'b' not in handle_mode or 'w' not in handle_mode:
        raise ValueError('BamIO.write() takes a file handle opened for binary write')

    try:
        fd, path = mkstemp(); close(fd)

        with open(path, 'w') as fh:
            sam_write(records, fh, reference, new_style)

        handle.write(samtools_view('-bS', path))
    finally:
        if exists(path):
            remove(path)

def sort(path):
    try:
        fd, new_path = mkstemp(); close(fd)

        samtools_sort(path, new_path)
        new_path += '.bam' # sort adds the .bam suffix automatically
        move(new_path, path)
    finally:
        if exists(new_path):
            remove(new_path)

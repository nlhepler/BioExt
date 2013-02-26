
from __future__ import division, print_function

from os import close, remove
from os.path import exists, getsize
from shutil import move
from tempfile import mkstemp

from pysam import (
    view as samtools_view,
    sort as samtools_sort
    )

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
            count = sam_write(records, fh, reference, new_style)

        if count:
            handle.write(samtools_view('-bS', path))
    finally:
        if exists(path):
            remove(path)

    return count


def sort(path):
    try:
        fd, tmp_path = mkstemp(); close(fd)

        if exists(path) and getsize(path):
            samtools_sort(path, tmp_path)
            tmp_path += '.bam'  # sort adds the .bam suffix automatically
            move(tmp_path, path)
    finally:
        if exists(tmp_path):
            remove(tmp_path)

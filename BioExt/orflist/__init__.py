#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import division, print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


__all__ = ['OrfList']


def _findall(subs, string):
    idxs = []
    for sub in subs:
        # because we cut the string to find the next sub,
        # we need to keep track of how much we've cut off
        prefix = 0
        while True:
            idx = string[prefix:].find(sub)
            # no more to find, break out
            if idx < 0:
                break
            idxs.append(prefix + idx)
            # increment to find the next,
            # possibly-overlapping, instance
            idx += 1
            prefix += idx
    return sorted(idxs)


class OrfList(object):

    def __init__(self, seq, include_stops=True):
        # deal with python 3
        try:
            seqtypes = (Seq, str, unicode)
        except NameError:
            seqtypes = (Seq, str)

        if isinstance(seq, SeqRecord):
            seq = seq.seq
        elif not isinstance(seq, seqtypes):
            raise ValueError('must provide either a SeqRecord, Seq, or str')

        self.__seq = seq

        # accumulate all the potential start codon indices
        useq = seq.upper()

        start_codons = ('ATG', 'AUG')
        stop_codons = ('TAG', 'UAG', 'TAA', 'UAA', 'TGA', 'UGA')

        start_idxs = _findall(start_codons, useq)
        stop_idxs = _findall(stop_codons, useq)

        offset = 3 if include_stops else 0

        orfs = []
        seq_end = len(useq) - offset
        for frame in range(3):
            frame_start_idxs = [idx for idx in start_idxs if idx % 3 == frame]
            # no start codons? skip!
            if len(frame_start_idxs) == 0:
                continue
            first_start = frame_start_idxs[0]
            frame_stop_idxs = [idx for idx in stop_idxs if idx % 3 == frame and idx > first_start]
            # allow ourselves to seek to the last position in the sequence,
            # if we've not already flagged it
            if len(frame_stop_idxs) == 0 or frame_stop_idxs[-1] != seq_end:
                frame_stop_idxs.append(seq_end)
            # find all orfs in this frame,
            # using `i` to save our last position in frame_stop_idxs,
            # and use it to shorten our next scan of stop_codons
            i = 0
            for start_idx in frame_start_idxs:
                for j, stop_idx in enumerate(frame_stop_idxs[i:]):
                    # make sure we grab at least two complete codons
                    if stop_idx > (start_idx + 2):
                        orfs.append((start_idx, stop_idx + offset))
                        i = j
                        break

        # sort the orfs from the longest to shortest
        self.__orfs = sorted(orfs, key=lambda orf: orf[1] - orf[0], reverse=True)

    def __getitem__(self, key):
        start, stop = self.__orfs[key]
        return self.__seq[start:stop]

    def __len__(self):
        return len(self.__orfs)

    def __contains__(self, key):
        return True if key >= 0 and key < len(self.__orfs) else False

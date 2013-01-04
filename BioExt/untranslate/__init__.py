# -*- coding: utf-8 -*-

from __future__ import division, print_function

from operator import itemgetter
from random import random

from Bio.Seq import Seq, translate as _translate
from Bio.SeqRecord import SeqRecord

from BioExt.collections import Counter


__all__ = [
    'UntranslationTable',
    'UniformUntranslationTable',
    'untranslate'
]


def _default_table(prior=0):
    table = {}
    for i in 'ACGT':
        for j in 'ACGT':
            for k in 'ACGT':
                cdn = ''.join((i, j, k))
                aa = _translate(cdn)
                if aa not in table:
                    table[aa] = Counter()
                table[aa].update({cdn: prior})
    # default the unkown amino acid to NNN in all cases
    table['X'] = Counter({'NNN': 1})
    return table


class UntranslationTable(object):

    def __init__(self, seq, prior=0):
        if isinstance(seq, SeqRecord):
            seq = str(seq.seq)
        elif isinstance(seq, Seq):
            seq = str(seq)
        elif not isinstance(seq, str):
            raise ValueError('seq must be of type SeqRecord, Seq, or str')
        table = _default_table(0)
        for i in range(0, len(seq), 3):
            j = i + 3
            if j > len(seq):
                continue
            cdn = seq[i:j].upper()
            aa = _translate(cdn)
            # skip unknown codons, they are irrelevant
            if aa == 'X':
                continue
            if cdn not in table[aa]:
                raise ValueError("sequence uses malformed alphabet '%s'" % cdn)
            table[aa][cdn] += 1
        for aa, cdns in table.items():
            total = prior * len(cdns) + sum(cdns.values())
            unif = 1. / len(cdns)
            cdf = []
            acc = 0.
            for cdn, count in sorted(cdns.items(), key=itemgetter(1)):
                if total:
                    pdf = (count + prior) / total
                else:
                    pdf = unif
                acc += pdf
                cdf.append((acc, cdn))
            table[aa] = cdf
        self.__table = table

    def __getitem__(self, key):
        if key not in self.__table:
            raise ValueError("unknown amino acid '%s'" % key)
        v = random()
        for p, cdn in self.__table[key]:
            if v < p:
                return cdn
        # catch numerical issues
        return self.__table[key][-1][1]


UniformUntranslationTable = UntranslationTable('', 1)

def untranslate(seq, table=UniformUntranslationTable):
    if not isinstance(table, UntranslationTable):
        raise ValueError('table must be an UntranslationTable object')
    r = []
    for aa in seq:
        r.append(table[aa])
    return ''.join(r)



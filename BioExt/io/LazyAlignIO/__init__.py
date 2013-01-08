
from os.path import abspath

from Bio import Alphabet
from Bio import SeqIO


__all__ = [
    'read'
    ]


def read(handle, format):
    return LazyMultipleSeqAlignment(handle, format)


class LazyMultipleSeqAlignment:

    def __init__(self, handle, format, alphabet=None):
        if alphabet is not None:
            if (not isinstance(alphabet, Alphabet.Alphabet)
                    or isinstance(alphabet, Alphabet.AlphabetEncoder)):
                raise ValueError('Invalide alphabet argument')
            self._alphabet = alphabet
        else:
            self._alphabet = Alphabet.single_letter_alphabet

        length = None
        nr = 0
        for r in SeqIO.parse(handle, format):
            if length is None:
                length = len(r)
            elif len(r) != length:
                raise ValueError('Sequences must all be the same length')
            nr += 1

        self._filename = abspath(handle.name)
        self._format = format
        self._alignment_length = length
        self._length = nr

    def __len__(self):
        return self._length

    def __iter__(self):
        with open(self._filename) as handle:
            for record in SeqIO.parse(handle, self._format):
                yield record

    def __getitem__(self, index):
        if isinstance(index, int):
            N = len(self)
            # handle negative indices
            if index < 0 and index >= -N:
                index = index % N
            elif index >= N:
                raise IndexError('index out of range')
            seq = None
            for i, seq_ in enumerate(self):
                if i == index:
                    seq = seq_
                    break
            return seq
        else:
            raise ValueError('invalid index')

    def get_alignment_length(self):
        return self._alignment_length

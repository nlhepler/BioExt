
from __future__ import division, print_function

from collections import defaultdict

import numpy as np

from Bio.Seq import Seq, translate as _translate
from Bio.SeqRecord import SeqRecord

from BioExt.aligner._align import _align, _compute_codon_matrices
from BioExt.scorematrix import ProteinScoreMatrix as _ProteinScoreMatrix


__all__ = ['Aligner']


def _protein_to_codon(protein_matrix):
    codon_matrix = np.ones((64, 64), dtype=float) * -1e4
    dletters = 'ACGT'
    pletters = protein_matrix.letters
    mapping = defaultdict(list)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                cdn = ''.join(dletters[l] for l in (i, j, k))
                aa = pletters.index(_translate(cdn))
                mapping[aa].append(16 * i + 4 * j + k)
    protein_matrix_ = protein_matrix.tondarray()
    M, N = protein_matrix_.shape
    for i in range(M):
        for k in mapping[i]:
            for j in range(N):
                for l in mapping[j]:
                    codon_matrix[k, l] = protein_matrix_[i, j]
    return dletters, codon_matrix


class Aligner:

    def __init__(
            self,
            score_matrix,
            open_insertion=None,
            extend_insertion=None,
            open_deletion=None,
            extend_deletion=None,
            miscall_cost=None,
            do_local=True,
            do_affine=True,
            do_codon=True
            ):

        if do_codon and not isinstance(score_matrix, _ProteinScoreMatrix):
            raise ValueError('codon alignment requires a protein score matrix')

        if do_codon:
            letters, score_matrix_ = _protein_to_codon(score_matrix)
        else:
            letters = score_matrix.letters
            score_matrix_ = score_matrix.tondarray()

        # set the default extension cost to 7.5% of the range
        magic = 40 / 3 # magic denominator for 7.5%
        ext_cost = (score_matrix_.max() - score_matrix_.min()) / magic
        # we take the negation of the minimum,
        # because the implementation assumes these values are penalties,
        # and subtracts them in all places
        min_score = -score_matrix_.min()

        # sane defaults for various penalties
        if open_insertion is None:
            open_insertion = 2.5 * min_score
        if extend_insertion is None:
            extend_insertion = ext_cost
        if open_deletion is None:
            open_deletion = 1.5 * min_score
        if extend_deletion is None:
            extend_deletion = ext_cost
        if miscall_cost is None:
            miscall_cost = min_score

        letters = letters.encode('utf8')
        char_map = np.zeros((256,), dtype=int)
        # this ensures that no matter context,
        # computing the codon index (16 * i + 4 * j + k) will always be 0
        # for any i, j, k in [0, 3]
        char_map[:] = -256
        for i, l in enumerate(letters):
            char_map[l] = i

        if do_codon:
            codon3x5, codon3x4, codon3x2, codon3x1 = _compute_codon_matrices(score_matrix_)
        else:
            empty = np.zeros((0, 0), dtype=float)
            codon3x5 = codon3x4 = codon3x2 = codon3x1 = empty

        self.__nchars = len(letters)
        self.__char_map = char_map
        self.__score_matrix = score_matrix_
        self.__open_insertion = open_insertion
        self.__extend_insertion = extend_insertion
        self.__open_deletion = open_deletion
        self.__extend_deletion = extend_deletion
        self.__miscall_cost = miscall_cost
        self.__do_local = do_local
        self.__do_affine = do_affine
        self.__do_codon = do_codon
        self.__codon3x5 = codon3x5
        self.__codon3x4 = codon3x4
        self.__codon3x2 = codon3x2
        self.__codon3x1 = codon3x1

    def __call__(
            self,
            ref,
            query,
            open_insertion=None,
            extend_insertion=None,
            open_deletion=None,
            extend_deletion=None,
            miscall_cost=None,
            do_local=None,
            do_affine=None
            ):

        # populate defaults from initialization
        if open_insertion is None:
            open_insertion = self.__open_insertion
        if extend_insertion is None:
            extend_insertion = self.__extend_insertion
        if open_deletion is None:
            open_deletion = self.__open_deletion
        if extend_deletion is None:
            extend_deletion = self.__extend_deletion
        if miscall_cost is None:
            miscall_cost = self.__miscall_cost
        if do_local is None:
            do_local = self.__do_local
        if do_affine is None:
            do_affine = self.__do_affine

        if isinstance(ref, SeqRecord):
            ref_ = str(ref.seq)
        elif isinstance(ref, Seq):
            ref_ = str(ref)
        else:
            ref_ = ref

        if isinstance(query, SeqRecord):
            query_ = str(query.seq)
        elif isinstance(query, Seq):
            query_ = str(query)
        else:
            query_ = query

        score, ref_aligned, query_aligned = _align(
            ref_,
            query_,
            self.__nchars,
            self.__char_map,
            self.__score_matrix,
            self.__score_matrix.shape[0],
            open_insertion,
            extend_insertion,
            open_deletion,
            extend_deletion,
            miscall_cost,
            do_local,
            do_affine,
            self.__do_codon,
            self.__codon3x5,
            self.__codon3x4,
            self.__codon3x2,
            self.__codon3x1
            )

        if isinstance(ref, SeqRecord):
            ref_aligned_ = SeqRecord(
                Seq(ref_aligned, ref.seq.alphabet),
                id=ref.id,
                name=ref.name,
                description=ref.description,
                dbxrefs=ref.dbxrefs,
                annotations=ref.annotations
                )
        elif isinstance(ref, Seq):
            ref_aligned_ = Seq(ref_aligned, ref.seq.alphabet)
        else:
            ref_aligned_ = ref_aligned

        if isinstance(query, SeqRecord):
            query_aligned_ = SeqRecord(
                Seq(query_aligned, query.seq.alphabet),
                id=query.id,
                name=query.name,
                description=query.description,
                dbxrefs=query.dbxrefs,
                annotations=query.annotations
                )
        elif isinstance(query, Seq):
            query_aligned_ = Seq(query_aligned, query.seq.alphabet)
        else:
            query_aligned_ = query_aligned

        return score, ref_aligned_, query_aligned_

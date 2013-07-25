
from __future__ import division, print_function

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport exp, log
from libc.stdlib cimport free

dtype = np.float64
itype = np.int

ctypedef np.float64_t dtype_t
ctypedef np.int_t itype_t

cdef char GAP
GAP = '-'

cdef dtype_t A_LARGE_NUMBER
A_LARGE_NUMBER = 1.0e100

cdef extern from "alignment.h":
    dtype_t AlignStrings(
        char *,
        char *,
        char **,
        char **,
        itype_t,
        itype_t *,
        dtype_t *,
        itype_t,
        char,
        dtype_t,
        dtype_t,
        dtype_t,
        dtype_t,
        dtype_t,
        itype_t,
        itype_t,
        itype_t,
        dtype_t *,
        dtype_t *,
        dtype_t *,
        dtype_t *,
        dtype_t *,
        dtype_t *,
        dtype_t *
        ) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
def choose(itype_t n, itype_t k):
    cdef itype_t i
    cdef dtype_t r

    r = 0.0
    for i in range(1, k + 1):
        r += log(n - k + i) - log(i)
    return exp(r)


@cython.boundscheck(False)
@cython.wraparound(False)
def _compute_codon_matrices(dtype_t[:, :] cost_matrix):

    cdef itype_t cdn1, cdn2, i, j, k, l
    cdef np.ndarray[dtype_t, ndim=2] codon3x5, codon3x4, codon3x2, codon3x1

    cdef dtype_t max100, max010, max001, max110, max101, max011, score

    cdef dtype_t penalty3x5, penalty3x4, penalty3x2, penalty3x1

    # these should be taken care of in alignment.c
    penalty3x4 = 0.0  # 1.0
    penalty3x5 = 2 * penalty3x4
    penalty3x2 = 0.0  # 1.0
    penalty3x1 = 2 * penalty3x2

    codon3x5 = np.zeros((64, 10 * 64), dtype=dtype)  # 64 codons, 10 possible placements
    codon3x4 = np.zeros((64, 4 * 64),  dtype=dtype)  # 64 codons, 4 possible placements
    codon3x2 = np.zeros((64, 16 * 3),  dtype=dtype)  # 16 dinucs, 3 possible placements
    codon3x1 = np.zeros((64, 4 * 3),   dtype=dtype)  # 4 nucs, 3 possible placements

    for cdn1 in range(64):
        for i in range(4):
            max100 = -A_LARGE_NUMBER
            max010 = -A_LARGE_NUMBER
            max001 = -A_LARGE_NUMBER
            for j in range(4):
                max110 = -A_LARGE_NUMBER
                max101 = -A_LARGE_NUMBER
                max011 = -A_LARGE_NUMBER
                for k in range(4):
                    cdn2 = 16 * i + 4 * j + k
                    score = cost_matrix[cdn1, cdn2]
                    # fill in 3x5 and 3x4 partial scoring matrices
                    for l in range(10):
                        codon3x5[cdn1, 10 * cdn2 + l] = score - penalty3x5
                    for l in range(4):
                        codon3x4[cdn1, 4 * cdn2 + l] = score - penalty3x4
                    # codon3x1 scores, 1 is the ith position, j/k over "wobble"
                    max100 = max(max100, cost_matrix[cdn1, 16 * i + 4 * j + k])
                    max010 = max(max010, cost_matrix[cdn1, 16 * j + 4 * i + k])
                    max001 = max(max001, cost_matrix[cdn1, 16 * j + 4 * k + i])
                    # codon3x2 score, 1s are in the ith and jth positions, k over "wobble"
                    max110 = max(max110, cost_matrix[cdn1, 16 * i + 4 * j + k])
                    max101 = max(max101, cost_matrix[cdn1, 16 * i + 4 * k + j])
                    max011 = max(max011, cost_matrix[cdn1, 16 * k + 4 * i + j])
                # fill codon3x2 partial scoring matrix
                codon3x2[cdn1, 12 * i + 3 * j + 0] = max110 - penalty3x2
                codon3x2[cdn1, 12 * i + 3 * j + 1] = max101 - penalty3x2
                codon3x2[cdn1, 12 * i + 3 * j + 2] = max011 - penalty3x2
            # fill codon3x1 partial scoring matrix
            codon3x1[cdn1, 3 * i + 0] = max100 - penalty3x1
            codon3x1[cdn1, 3 * i + 1] = max010 - penalty3x1
            codon3x1[cdn1, 3 * i + 2] = max001 - penalty3x1

    return codon3x5, codon3x4, codon3x2, codon3x1


@cython.boundscheck(False)
@cython.wraparound(False)
def _align(
        bytes b_ref,
        bytes b_query,
        itype_t char_count,
        np.ndarray[itype_t] char_map,
        np.ndarray[dtype_t, ndim=2, mode='c'] cost_matrix,
        itype_t cost_stride,
        dtype_t open_insertion,
        dtype_t extend_insertion,
        dtype_t open_deletion,
        dtype_t extend_deletion,
        dtype_t miscall_cost,
        itype_t do_local,
        itype_t do_affine,
        itype_t do_codon,
        np.ndarray[dtype_t, ndim=2, mode='c'] codon3x5,
        np.ndarray[dtype_t, ndim=2, mode='c'] codon3x4,
        np.ndarray[dtype_t, ndim=2, mode='c'] codon3x2,
        np.ndarray[dtype_t, ndim=2, mode='c'] codon3x1,
        np.ndarray[dtype_t, ndim=1, mode='c'] score_matrix,
        np.ndarray[dtype_t, ndim=1, mode='c'] deletion_matrix,
        np.ndarray[dtype_t, ndim=1, mode='c'] insertion_matrix):

    # cast from unicode to char *
    cdef char * ref = b_ref
    cdef char * query = b_query

    cdef char * ref_aligned = NULL
    cdef char * query_aligned = NULL
    cdef bytes b_ref_aligned
    cdef bytes b_query_aligned
    cdef dtype_t score

    if do_codon and len(b_ref) % 3 != 0:
        raise ValueError('when do_codon = True, len(ref) must be a multiple of 3')

    try:
        score = AlignStrings(
            ref, query,
            &ref_aligned, &query_aligned,
            char_count,
            <itype_t *> char_map.data,
            <dtype_t *> cost_matrix.data,
            cost_stride,
            GAP,
            open_insertion, extend_insertion,
            open_deletion, extend_deletion,
            miscall_cost,
            do_local, do_affine, do_codon,
            <dtype_t *> codon3x5.data,
            <dtype_t *> codon3x4.data,
            <dtype_t *> codon3x2.data,
            <dtype_t *> codon3x1.data,
            <dtype_t *> score_matrix.data,
            <dtype_t *> deletion_matrix.data,
            <dtype_t *> insertion_matrix.data)

        if ref_aligned == NULL or query_aligned == NULL:
            free(ref_aligned)
            free(query_aligned)
            raise MemoryError('memory allocation error in AlignStrings(...)')

        # cast char * back to unicode
        b_ref_aligned = ref_aligned
        b_query_aligned = query_aligned

    finally:
        free(ref_aligned)
        free(query_aligned)

    return score, b_ref_aligned, b_query_aligned

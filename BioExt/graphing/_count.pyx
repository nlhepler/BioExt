
from __future__ import division, print_function

import numpy as np
cimport numpy as np
cimport cython

dtype = np.float64
itype = np.int

ctypedef np.float64_t dtype_t
ctypedef np.int_t itype_t
ctypedef np.int8_t ctype_t


@cython.boundscheck(False)
@cython.wraparound(False)
def _count(
        object alignment,
        np.ndarray[itype_t, ndim=1, mode='c'] cols,
        object alphabet,
        np.ndarray[dtype_t, ndim=2, mode='c'] values):

    cdef itype_t nval = values.shape[0]
    cdef itype_t nchar = values.shape[1]
    cdef itype_t ncol = 0
    cdef itype_t npos = 0
    cdef itype_t i = 0
    cdef np.ndarray[dtype_t, ndim=2, mode='c'] counts = None

    if cols is not None:
         ncol = cols.shape[0]

    for seq in alignment:
        if counts is None:
            npos = len(seq)
            if cols is None:
                ncol = npos
            counts = np.zeros((ncol, nchar), dtype=dtype)
        if len(seq) != npos:
            raise ValueError("provided sequences are not all the same length")
        if cols is None:
            for j in range(ncol):
                for k in range(nval):
                    if seq[j] == alphabet[k]:
                        np.add(counts[j, :], values[k, :], counts[j, :])
                        break
        else:
            for j in range(ncol):
                if cols[j] < 0 or cols[j] >= npos:
                    raise ValueError("nonexistent columns specified")
                for k in range(nval):
                    if seq[cols[j]] == alphabet[k]:
                        np.add(counts[j, :], values[k, :], counts[j, :])
                        break
        i += 1

    return counts

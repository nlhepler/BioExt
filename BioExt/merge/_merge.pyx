
from __future__ import division, print_function

from libc.stdlib cimport free, malloc, realloc

from re import compile as re_compile
from sys import stderr

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


cdef extern from "merge.h":
    ctypedef struct opts_t:
        int min_overlap
        int min_reads
        int tol_gaps
        int tol_ambigs

    ctypedef struct triple_t:
        int col
        int ins
        char nuc

    ctypedef struct aligned_t:
        triple_t * data
        int len
        int lpos
        int rpos
        int ncontrib

    cdef char nuc2bits(char nuc) nogil

    cdef char bits2nuc(char bits) nogil

    cdef aligned_t * merge__(
        int nreads,
        aligned_t * reads,
        opts_t * opts,
        int * nclusters
        ) nogil

    cdef void aligned_destroy(aligned_t * read) nogil


cdef aligned_t labelstream(record):

    cdef int col = record.annotations['position'] - 1
    cdef int ins = 0
    cdef int rpos = 0
    cdef int idx = 0
    cdef triple_t * data = <triple_t *>malloc(len(record) * sizeof(triple_t))
    cdef aligned_t record_

    if data == NULL:
        raise MemoryError("failed to allocate memory for record data")

    regexp = re_compile(r'([0-9]+)([M=XID])')
    seq_ = iter(record)

    for m in regexp.finditer(record.annotations['CIGAR']):
        num, mode = int(m.group(1)), m.group(2)
        for _ in range(num):
            rpos += 1

            if mode in 'M=X':
                col += 1
                ins = 0
                nuc = next(seq_)
            elif mode == 'I':
                ins += 1
                nuc = next(seq_)
            elif mode == 'D':
                col += 1
                ins = 0
                continue

            data[ idx ].col = col
            data[ idx ].ins = ins
            data[ idx ].nuc = nuc2bits(nuc)

            idx += 1

    record_.data = data
    record_.len = len(record)
    record_.lpos = record.annotations['position']
    record_.rpos = rpos
    record_.ncontrib = 1

    return record_


# for AlignedRead from pysam
cdef aligned_t labelstream_(read):

    cdef int col = read.pos - 1
    cdef int idx = 0
    cdef int ins = 0
    cdef int rpos = 0
    cdef triple_t * data = <triple_t *>malloc(read.qlen * sizeof(triple_t))
    cdef aligned_t read_

    if data == NULL:
        raise MemoryError("failed to allocate memory for read data")

    for mode, num in read.cigar:
        for _ in range(num):

            assert idx < read.qlen, "Invalid access!"

            rpos += 1

            if mode in (0, 7, 8):
                col += 1
                ins = 0
            elif mode == 1:
                ins += 1
            elif mode == 2:
                col += 1
                ins = 0
                continue

            data[ idx ].col = col
            data[ idx ].ins = ins
            data[ idx ].nuc = nuc2bits(read.seq[idx])

            idx += 1

    read_.data = data
    read_.len = read.qlen
    read_.lpos = read.pos
    read_.rpos = rpos
    read_.ncontrib = 1

    return read_


cdef _to_seqrecord(int idx, aligned_t * cluster):
    if not cluster.len:
        raise ValueError(
            'cluster has length of 0, which is very incorrect'
            )

    next_col = cluster.data[0].col + 1
    seq = [<Py_UNICODE>bits2nuc(cluster.data[0].nuc)]
    cigar = []
    mode = 'I' if cluster.data[0].ins else 'M'
    num = 1

    for i in range(1, cluster.len):
        # handle a skip in the column
        if cluster.data[i].col > next_col:
            if mode == 'D':
                num += cluster.data[i].col - next_col
            else:
                cigar.append('{0:d}{1:s}'.format(num, mode))
                num = cluster.data[i].col - next_col
                mode = 'D'
        # keep track of the next column, seq, and mode_
        next_col = cluster.data[i].col + 1
        seq.append(<Py_UNICODE>bits2nuc(cluster.data[i].nuc))
        mode_ = 'I' if cluster.data[i].ins else 'M'

        if mode_ == mode:
            num += 1
        else:
            cigar.append('{0:d}{1:s}'.format(num, mode))
            num = 1
            mode = mode_

    cigar.append('{0:d}{1:s}'.format(num, mode))

    name = 'cluster_{0:d}_{1:d}r'.format(idx, cluster.ncontrib)

    return SeqRecord(
        Seq(''.join(seq), single_letter_alphabet),
        id=name,
        name=name,
        description=name,
        annotations={'CIGAR': ''.join(cigar), 'position': cluster.lpos}
        )


cdef _merge(bam, reads, int min_overlap, int min_reads, int tol_gaps, int tol_ambigs):
    cdef int nreads = 512
    cdef aligned_t * reads_ = <aligned_t *>malloc(nreads * sizeof(aligned_t))
    cdef opts_t * opts = <opts_t *>malloc(sizeof(opts_t))
    cdef aligned_t * clusters = NULL
    cdef int nclusters
    cdef int idx = 0
    cdef int i = 0
    opts.min_overlap = min_overlap
    opts.min_reads = min_reads
    opts.tol_gaps = tol_gaps
    opts.tol_ambigs = tol_ambigs

    if reads_ == NULL:
        raise MemoryError("failed to allocate memory for reads")

    if opts == NULL:
        raise MemoryError("failed to allocate memory for opts")

    try:
        for read in reads:
            if bam:
                if read.qlen < min_overlap:
                    continue
                reads_[idx] = labelstream_(read)
            else:
                if len(read) < min_overlap:
                    continue
                reads_[idx] = labelstream(read)

            idx += 1

            if idx >= nreads:
                nreads *= 2
                reads_ = <aligned_t *>realloc(reads_, nreads * sizeof(aligned_t))

        nreads = idx
        reads_ = <aligned_t *>realloc(reads_, nreads * sizeof(aligned_t))

        clusters = merge__(nreads, reads_, opts, &nclusters)

        if clusters == NULL:
            raise MemoryError(
                "failed to merge clusters, likely a memory allocation failure"
                )

        clusters_ = []

        for i in range(nclusters):
            clusters_.append(_to_seqrecord(i + 1, &clusters[i]))

    finally:
        free(opts)

        for i in range(nreads):
            aligned_destroy(&reads_[i])

        free(reads_)

        for i in range(nclusters):
            if clusters[i].ncontrib > 1:
                aligned_destroy(&clusters[i])

        free(clusters)

    return clusters_


def merge(reads, min_overlap, min_reads=1, tol_gaps=True, tol_ambigs=True):
    return _merge(False, reads, min_overlap, min_reads, tol_gaps, tol_ambigs)


def merge_(reads, min_overlap, min_reads=1, tol_gaps=True, tol_ambigs=True):
    return _merge(True, reads, min_overlap, min_reads, tol_gaps, tol_ambigs)

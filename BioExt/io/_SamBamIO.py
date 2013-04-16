
from __future__ import division, print_function

from hashlib import md5
from os.path import exists
from re import compile as re_compile

from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt.misc import compute_cigar

import pysam


_VALID_CHARS = re_compile(r'[^A-Z]')
_CIGAR_MODES = 'MIDNSHP=X'


def _from_cigarstring(cigarstr):
    regexp = re_compile(r'([0-9]+)([M=XID])')
    cigar = []
    for m in regexp.finditer(cigarstr):
        num, mode = int(m.group(1)), m.group(2).upper()
        mode_ = _CIGAR_MODES.find(mode)
        cigar.append((mode_, num))
    return cigar


def _to_cigarstring(cigar):
    return ''.join(
        '{0:d}{1:s}'.format(num, _CIGAR_MODES[mode])
        for mode, num in cigar
        )


def _to_seqrecord(samfile, read):

    seq = Seq(
        read.seq.decode('ascii'),
        single_letter_alphabet
        )

    qname = read.qname
    rname = samfile.getrname(read.tid)

    annotations = {}
    annotations['sam_flag'] = read.flag
    annotations['reference_name'] = rname
    annotations['position'] = read.pos
    annotations['mapping_quality'] = read.mapq
    annotations['CIGAR'] = _to_cigarstring(read.cigar)
    annotations['reference_next'] = read.rnext
    annotations['position_next'] = read.pnext
    annotations['template_length'] = read.tlen

    if read.qual:
        letter_annotations = {}
        letter_annotations['phred_quality'] = [ord(q) - 33 for q in read.qual]
    else:
        letter_annotations = None

    record = SeqRecord(
        seq,
        id=qname,
        name=qname,
        description=qname,
        annotations=annotations,
        letter_annotations=letter_annotations
        )

    return record


def _find_tid(header, name):

    try:
        tid = next(
            i for i, sq
            in enumerate(header['SQ'])
            if sq['SN'] == name
            )
    except:
        tid = -1

    return tid


def _from_seqrecord(header, record):

    qname = record.name
    tid = _find_tid(
        header,
        record.annotations.get('reference_name', '*')
        )
    seq = _VALID_CHARS.sub('', str(record.seq).upper()).encode('ascii')
    cigar = _from_cigarstring(record.annotations['CIGAR'])
    pos = record.annotations['position']
    flag = record.annotations.get('sam_flag', 0)
    mapq = record.annotations.get('mapping_quality', 255)
    rnext = _find_tid(
        header,
        record.annotations.get('reference_next', '*')
        )
    pnext = record.annotations.get('position_next', -1)
    tlen = record.annotations.get('template_length', 0)
    qual = record.letter_annotations.get('phred_quality', None)

    if 'edit_distance' in record.annotations:
        tags = [('NM', record.annotations['edit_distance'])]
    else:
        tags = []

    read = pysam.AlignedRead()

    read.qname = qname
    read.tid = tid
    read.pos = pos
    read.seq = seq
    read.cigar = cigar
    read.flag = flag
    read.mapq = mapq
    read.rnext = rnext
    read.pnext = pnext
    read.tlen = tlen
    read.qual = qual
    read.tags = tags

    # TODO: enable this when the method exists
    # read.validate()

    return read


def _parse(mode, path, index=True):
    try:
        # make sure the index is up to date
        if 'b' in mode and index:
            pysam.index(path)

        samfile = pysam.Samfile(path, mode)

        for read in samfile.fetch():
            yield _to_seqrecord(samfile, read)
    finally:
        if samfile is not None:
            samfile.close()


def _write(mode, records, path, reference, new_style):
    try:
        count = 0
        header = {
            'HD': {'VN': '1.4', 'SO': 'unknown'},
            }
        if reference is not None:
            hasher = md5()
            hasher.update(
                _VALID_CHARS.sub('', str(reference.seq).upper()).encode('ascii')
                )
            header['SQ'] = [{
                'SN': reference.name,
                'LN': len(reference),
                'M5': hasher.hexdigest()
                }]
        samfile = pysam.Samfile(path, mode, header=header)

        if reference is not None and isinstance(records, MultipleSeqAlignment):
            def iterate(records):
                for record in records:
                    if ('CIGAR' in record.annotations and
                            'position' in record.annotations):
                        yield record
                    else:
                        yield compute_cigar(reference, record, new_style)
        else:
            def iterate(records):
                for record in records:
                    yield record

        for record in iterate(records):
            if samfile.write(_from_seqrecord(header, record)):
                count += 1
    finally:
        samfile.close()

    return count

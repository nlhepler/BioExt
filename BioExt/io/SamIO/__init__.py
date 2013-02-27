
from __future__ import division, print_function

from hashlib import md5
from re import compile as re_compile

from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt.misc import compute_cigar


__all__ = [
    'parse',
    'write'
    ]


def parse(handle):
    hdr = True
    for line in handle:
        fields = line.split('\t')

        # skip header lines
        if fields[0].upper() not in ('@HD', '@SQ', '@RG', '@PG', '@CO'):
            hdr = False
        if hdr:
            continue

        # TODO: properly handle multi-segment SAM files
        name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fields[:11]

        annotations = {}
        annotations['sam_flag'] = int(flag)
        annotations['reference_name'] = rname
        annotations['position'] = int(pos)
        annotations['mapping_quality'] = int(mapq)
        annotations['CIGAR'] = cigar
        annotations['reference_next'] = rnext
        annotations['position_next'] = int(pnext)
        annotations['template_length'] = int(tlen)

        # optional fields
        for field in fields[11:]:
            fname, ftype, fval = field.split(':', 2)
            if fname.upper() == 'NM':
                annotations['edit_distance'] = int(fval)

        if qual != '*':
            letter_annotations = {}
            letter_annotations['phred_quality'] = [ord(q) - 33 for q in qual]
        else:
            letter_annotations = None

        record = SeqRecord(
            Seq(seq, single_letter_alphabet),
            id=name,
            name=name,
            description=name,
            annotations=annotations,
            letter_annotations=letter_annotations
            )

        yield record


def write(records, handle, reference=None, new_style=False):
    re_alph = re_compile(r'[^A-Z]')

    # header lines
    header = ['@HD\tVN:1.4\tSO:unknown\n']
    if reference is not None:
        reflen = len(reference)
        refname = reference.name
        hasher = md5()
        hasher.update(re_alph.sub('', str(reference.seq).upper()).encode('utf8'))
        header.append(
            '@SQ\tSN:%s\tLN:%d\tM5:%s\n' % (
                refname,
                reflen,
                hasher.hexdigest()
                )
            )

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

    count = 0
    for record in iterate(records):
        # print the header before the first record
        if header:
            for line in header:
                handle.write(line)
            header = None

        # retrieve the optional annotations
        flag = record.annotations.get('sam_flag', 0)
        mapq = record.annotations.get('mapping_quality', 255)
        rname = record.annotations.get('reference_name', '*')
        rnext = record.annotations.get('reference_next', '*')
        pnext = record.annotations.get('position_next', 0)
        tlen = record.annotations.get('template_length', 0)
        # 9 (+33) corresponds to '*'
        quals = record.letter_annotations.get('phred_quality', [9])

        fields = [
            record.name,
            '%d' % flag,
            '%s' % rname,
            '%d' % record.annotations['position'],
            '%d' % mapq,
            record.annotations['CIGAR'],
            rnext,
            '%d' % pnext,
            '%d' % tlen,
            re_alph.sub('', str(record.seq).upper()),
            ''.join(chr(q + 33) for q in quals)
            ]
        # if we have edit distance, add that too
        if 'edit_distance' in record.annotations:
            fields.append('NM:i:%d' % record.annotations['edit_distance'])
        # write!
        handle.write('\t'.join(fields))
        handle.write('\n')
        # increment!
        count += 1

    return count

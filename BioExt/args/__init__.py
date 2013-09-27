
from __future__ import division, print_function


__all__ = [
    'add_alphabet',
    'add_reference',
    'add_scorematrix'
    ]


def add_alphabet(parser, shortarg='-a', longarg='--alphabet'):
    parser.add_argument(
        shortarg, longarg,
        metavar='ALPHABET',
        choices=('amino', 'dna', 'codon'),
        default='codon',
        help='perform an alignment using one of {{{0}}} [default=codon]'.format(
            ', '.join(('amino', 'dna', 'codon'))
            )
        )
    return parser


def add_reference(parser, shortarg='-r', longarg='--reference'):
    from argparse import ArgumentTypeError
    from Bio import SeqIO
    from BioExt.references import hxb2, nl4_3

    references = {
        'HXB2_env': hxb2.env,
        'HXB2_gag': hxb2.gag,
        'HXB2_int': hxb2.int,
        'HXB2_nef': hxb2.nef,
        'HXB2_pol': hxb2.pol,
        'HXB2_pr': hxb2.pr,
        'HXB2_prrt': hxb2.prrt,
        'HXB2_rev': hxb2.rev,
        'HXB2_rt': hxb2.rt,
        'HXB2_tat': hxb2.tat,
        'HXB2_vif': hxb2.vif,
        'HXB2_vpr': hxb2.vpr,
        'HXB2_vpu': hxb2.vpu,
        'NL4-3_prrt': nl4_3.prrt
        }

    def reference(string):
        if string in references:
            return references[string].load()
        try:
            with open(string) as handle:
                ref = next(SeqIO.parse(handle, 'fasta'))
            return ref
        except:
            msg = "'{0}' does not exist or is not a valid FASTA file".format(string)
            raise ArgumentTypeError(msg)

    parser.add_argument(
        shortarg, longarg,
        metavar='REFERENCE',
        type=reference,
        help='REFERENCE FASTA file or {{{0}}}'.format(', '.join(references.keys()))
        )
    return parser


def add_scorematrix(parser, shortarg='-m', longarg='--score-matrix'):
    import BioExt.scorematrices
    from BioExt.scorematrices import (
        DNAScoreMatrix,
        LazyScoreMatrix,
        ProteinScoreMatrix
        )

    score_matrices = {}
    for obj in dir(BioExt.scorematrices):
        val = getattr(BioExt.scorematrices, obj)
        if isinstance(val, (DNAScoreMatrix, LazyScoreMatrix, ProteinScoreMatrix)):
            score_matrices[str(val)] = val

    parser.add_argument(
        shortarg, longarg,
        metavar='SCOREMATRIX',
        type=lambda s: score_matrices.get(s, s),
        choices=sorted(score_matrices.values(), key=str),
        default=score_matrices['BLOSUM62'],
        help='parameterize using one of {{{0}}} [default=BLOSUM62]'.format(
            ', '.join(score_matrices.keys())
            )
        )
    return parser

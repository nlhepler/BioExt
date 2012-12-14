
from __future__ import division, print_function

from sklearn.externals.joblib import delayed, Parallel

from Bio import SeqIO

from BioExt.aligner import Aligner
from BioExt.io import BamIO
from BioExt.misc import compute_cigar
from BioExt.scorematrix import HIV_BETWEEN_F

def align(aln, ref, seq):
    score, ref, seq = aln(ref, seq)
    seq_ = compute_cigar(ref, seq)
    return score, seq_

def main(reffile, seqsfile, outfile):
    hbf = HIV_BETWEEN_F.load()
    aln = Aligner(hbf)

    with open(reffile) as fh:
        ref = SeqIO.read(fh, 'fasta')

    with open(seqsfile) as fh:
        results = Parallel(
            n_jobs=-1,
            verbose=0,
            pre_dispatch='10*n_jobs'
            )(
                delayed(align)(aln, ref, seq)
                for seq in SeqIO.parse(fh, 'fasta')
                )

        scores, seqs = zip(*results)

    with open(outfile, 'wb') as fh:
        BamIO.write(seqs, fh, ref)

    BamIO.sort(outfile)

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(*sys.argv[1:]))

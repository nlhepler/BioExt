
# These functions are borrowed from PacBio's GenomicConsensus software,
# but python3-compatible and with minimal dependencies.
# http://www.github.com/PacificBiosciences/GenomicConsensus

from __future__ import division, print_function

from operator import itemgetter

from BioExt.collections import OrderedDict

import ConsensusCore as cc


fst = itemgetter(0)
snd = itemgetter(1)
NoQVsModelParams = OrderedDict(
    Match=0.0,
    Mismatch=-1.21730327606,
    MismatchS=0.0,
    Branch=-0.371355384588,
    BranchS=0.0,
    DeletionN=-0.250208973885,
    DeletionWithTag=0.0,
    DeletionWithTagS=0.0,
    Nce=-0.250370770693,
    NceS=0.0,
    Merge=-0.371355384588,
    MergeS=0.0
)


def uniqueSingleBaseMutations(tpl, positions=None):
    allBases = "ACGT"
    prevTplBase = None
    positions = positions or range(len(tpl))
    for tplStart in positions:
        tplBase = tpl[tplStart]
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation(cc.SUBSTITUTION, tplStart, subsBase)
        # Insertions---only allowing insertions that are not cognate
        # with the previous base.
        for insBase in allBases:
            if insBase != prevTplBase:
                yield cc.Mutation(cc.INSERTION, tplStart, insBase)
        # Deletion--only allowed if refBase does not match previous tpl base
        if tplBase != prevTplBase:
            yield cc.Mutation(cc.DELETION, tplStart, "-")
        prevTplBase = tplBase


def nearbyMutations(mutations, tpl, neighborhoodSize):
    mutationPositions = map(cc.Mutation.Position, mutations)
    nearbyPositions = set()
    for mp in mutationPositions:
        nearbyPositions.update(
            range(
                max(0, mp - neighborhoodSize),
                min(len(tpl), mp + neighborhoodSize)
                )
            )
    return uniqueSingleBaseMutations(tpl, sorted(nearbyPositions))


def bestSubset(mutationsAndScores, separation):
    """
    Given a list of (mutation, score) tuples, this utility method
    greedily chooses the highest scoring well-separated elements.  We
    use this to avoid applying adjacent high scoring mutations, which
    are the rule, not the exception.  We only apply the best scoring one
    in each neighborhood, and then revisit the neighborhoods after
    applying the mutations.
    """
    input = mutationsAndScores[:]
    output = []

    while input:
        best = max(input, key=snd)
        output.append(best)
        nStart = best[0].Position() - separation
        nEnd = best[0].Position() + separation
        for t in input[:]:
            if nStart <= t[0].Position() <= nEnd:
                input.remove(t)

    return output


def refineConsensus(mms, niter=20):
    # some sane defaults
    mutationNeighborhood = 20
    mutationSeparation = 10

    favorableMutationsAndScores = None
    converged = False

    for _ in range(niter):

        if favorableMutationsAndScores is None:
            mutationsToTry = uniqueSingleBaseMutations(mms.Template())
        else:
            favorableMutations = list(map(fst, favorableMutationsAndScores))
            mutationsToTry = nearbyMutations(
                favorableMutations,
                mms.Template(),
                mutationNeighborhood
                )

        favorableMutationsAndScores = [
            (m, mms.Score(m))
            for m
            in filter(mms.FastIsFavorable, mutationsToTry)
            ]

        if favorableMutationsAndScores:
            bestMutations = list(map(
                fst,
                bestSubset(favorableMutationsAndScores, mutationSeparation)
                ))
            mms.ApplyMutations(bestMutations)
        else:
            # If we can't find any favorable mutations, our work is done.
            converged = True
            break

    return mms.Template(), converged


def extractFeatures(read):
    args = [read]
    for _ in range(5):
        args.append(cc.FloatFeature(len(read)))
    return cc.QvSequenceFeatures(*args)


from __future__ import division, print_function

import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import repeat
from operator import itemgetter
from os import close
from os.path import dirname, join
from re import compile as re_compile
from tempfile import mkstemp

from Bio.Alphabet import (
    Gapped,
    HasStopCodon,
    _verify_alphabet
    )
from Bio.Alphabet.IUPAC import (
    ambiguous_dna,
    ambiguous_rna,
    extended_protein
    )

from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    FormatStrFormatter,
    FuncFormatter
    )

from BioExt.graphing._basefont import Basefont
from BioExt.misc import _GAP, _STOP


__all__ = [
    'graph_coverage_majority',
    'graph_logo'
    ]


FONT_PATH = join(dirname(dirname(__file__)), 'data', 'fonts', 'ttf')


# update the fontManager to handle the Roboto font installed with hy454
# fontManager.ttffiles.extend(findSystemFonts([FONT_PATH]))
# fontManager.ttflist = createFontList(fontManager.ttffiles)
ROBOTO_REGULAR = FontProperties(
    fname=join(FONT_PATH, 'Roboto-Regular.ttf')
    )


# BLACK = '#000000'
GREY = '#969cb0'
# GREEN = '#b8ff25'
# ORANGE = '#ff6e27'
# BLUE = '#189cff'
# RED = '#e80c5b' # '#e80c7a'

# android color swatch
LIGHT_GREY = '#F2F2F2'
LIGHT_BLUE = '#33B5E5'
LIGHT_PURPLE = '#AA66CC'
LIGHT_GREEN = '#99CC00'
LIGHT_ORANGE = '#FFBB33'
LIGHT_RED = '#FF4444'
DARK_GREY = '#DDDDDD'
DARK_BLUE = '#0099CC'
DARK_PURPLE = '#9933CC'
DARK_GREEN = '#669900'
DARK_ORANGE = '#FF8800'
DARK_RED = '#CC0000'

# default grey
DNA_COLORS = defaultdict(
    repeat(GREY).__next__, {
        'A': LIGHT_GREEN,
        'C': LIGHT_BLUE,
        'G': LIGHT_ORANGE,
        'T': LIGHT_RED,
        'U': LIGHT_RED
        }
    )
AMINO_COLORS = defaultdict(
    repeat(GREY).__next__,
    [(l, LIGHT_GREEN) for l in 'KRH'] +
    [(l, LIGHT_BLUE) for l in 'DE'] +
    [(l, LIGHT_RED) for l in 'AVLIPWFM']
    )


def _adjust_spines_outward(ax, spines, points):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', points))
            # do not call set_smart_bounds, it removes ticks for some reason
            # spine.set_smart_bounds(True)


def _max_nonzero_min(values, default=0):
    vmax = max(values)
    vmin = vmax
    for v in values:
        if v > 0 and v < vmin:
            vmin = v
    vmin = default if vmin == vmax else vmin
    return vmax, vmin


def graph_coverage_majority(
    alignment,
    mode,
    nseq=None,
    filename=None,
    dpi=None, figsize=None, format='pdf', transparent=True,
    refidx=None
):
    if not mode:
        mode = 'both'

    if mode not in ('coverage', 'majority', 'both'):
        raise ValueError("mode must be one of 'coverage', 'majority', or 'both'")

    if filename is None:
        fd, filename = mkstemp()
        close(fd)

    if figsize is None:
        # actually results in a graph that is almost 4" high
        # and 6.4" wide -- almost the golden ratio
        figsize = (6, 6)

    if refidx is not None:
        msa = alignment
        alignment = msa[:refidx]
        alignment.extend(msa[refidx + 1:])
        del msa

    if nseq is None:
        nseq = len(alignment)

    fst_pos = 0
    lst_pos = alignment.get_alignment_length()

    # if we don't match the whole reference,
    # cut off the head, tail of gaps
    if False:  # strip:
        leading = re_compile(r'^-*')
        trailing = re_compile(r'-*$')
        l = None
        t = None
        for r in alignment:
            seq = str(r.seq)
            # *ing.search should always match
            l_ = len(leading.search(seq).group(0))
            t_ = len(trailing.search(seq).group(0))
            if l is None or l_ < l:
                l = l_
            if t is None or t_ < t:
                t = t_
        # if we're greater than 0, remove leading, trailing gaps
        if l:
            fst_pos += l
        if t:
            lst_pos -= t

    npos = lst_pos - fst_pos
    xs = np.arange(fst_pos, lst_pos) + 1

    # some heuristic to prevent too many ticks,
    # is complicated by the start and end dynamically business
    # 0.477122 jumps an order 10 at 30% of the next order 10
    specialK = 0.477122
    xdiv = 10 ** int(np.log10(lst_pos) - specialK) or 1
    xsep = xdiv * int(npos / 5 / xdiv + 1)
    xstart = (fst_pos // xsep) * xsep
    xend = lst_pos + xsep
    xticks = np.arange(xstart, xend, xsep)
    # prevent the first two and last two ticks from being too close
    if len(xticks) > 1:
        fix = False
        min_delta = xsep / 4
        delta = xticks[1] - fst_pos
        if delta and delta < min_delta:
            xstart += xsep
            fix = True
        delta = lst_pos - xticks[-2]
        if delta and delta < min_delta:
            xend -= xsep
            fix = True
        if fix:
            xticks = np.arange(xstart, xend, xsep)

    xticks[0] = fst_pos + 1
    xticks[-1] = lst_pos

    ydiv = 10 ** int(np.log10(nseq) - specialK) or 1
    ysep = ydiv * int(nseq / 5 / ydiv + 1)
    yticks = np.arange(0, nseq + ysep, ysep)

    if len(yticks) > 1:
        if nseq - yticks[-2] < ysep / 4:
            yticks = np.arange(0, nseq, ysep)

    yticks[-1] = nseq

    alph = set([])
    for seq in alignment:
        alph.update(ltr.upper() for ltr in seq[fst_pos:lst_pos])

    lmap = dict(
        (ltr, idx)
        for idx, ltr
        in enumerate(
            ltr_
            for ltr_
            in sorted(alph)
            if ltr_ != _GAP
            )
        )

    counts = np.zeros((len(lmap), npos), dtype=int)

    for seq in alignment:
        for i, l in enumerate(seq[fst_pos:lst_pos]):
            if l == _GAP:
                continue
            counts[lmap[l.upper()], i] += 1

    fig = plt.figure(figsize=figsize, dpi=dpi)

    # golden rectangle!
    rect = 0.2, 0.2, 1, 0.618
    ax1 = fig.add_axes(rect)

    # move the axis spines off the data
    _adjust_spines_outward(ax1, ('left', 'right'), 18)

    if mode in ('majority', 'both'):
        if mode == 'majority':
            majorities = counts.max(axis=0) / counts.sum(axis=0)
        else:
            majorities = counts.max(axis=0)

        lines = ax1.plot(
            xs,
            majorities,
            color=LIGHT_RED,
            linewidth=1.0,
            zorder=-2
            )

        for l in lines:
            l.set_clip_on(False)

    if mode in ('coverage', 'both'):
        coverages = counts.sum(axis=0)
        lines = ax1.plot(
            xs,
            coverages,
            color=LIGHT_BLUE,
            linewidth=1.0,
            zorder=-1
            )

        for l in lines:
            l.set_clip_on(False)

    extra_artists = []

    if mode == 'both':
        # create a proxy artist for legend, PolyCollections don't work (heights)
        p1 = Line2D([0, 1], [0, 1], color=LIGHT_BLUE, linewidth=1.0)
        # create a proxy artist for legend, [Lines2D] don't work (majorities)
        p2 = Line2D([0, 1], [0, 1], color=LIGHT_RED, linewidth=1.0)
        leg = ax1.legend(
            [p1, p2],
            ['Coverage', 'Majority proportion'],
            bbox_to_anchor=(0.5, -0.15),
            loc=9,
            ncol=2,
            prop=ROBOTO_REGULAR,
            borderpad=0
        )
        leg.legendPatch.set_alpha(0.)
        extra_artists.append(leg)
    # labels
        ax1.set_ylabel('No. of sequences', fontproperties=ROBOTO_REGULAR)
    elif mode == 'coverage':
        ax1.set_ylabel('Coverage', fontproperties=ROBOTO_REGULAR)
    elif mode == 'majority':
        ax1.set_ylabel('Majority proportion', fontproperties=ROBOTO_REGULAR)

    ax1.set_xlabel('Reference sequence position', fontproperties=ROBOTO_REGULAR)
    ax1.set_xlim((fst_pos + 1, lst_pos))

    major_ticks = ax1.xaxis.get_major_ticks(len(xticks))

    # disable the first and last tick on the x-axis,
    # they're redundant and ugly (esp if minH > 0)
    major_ticks[0].tick1On = False
    major_ticks[-1].tick1On = False

    if mode in ('coverage', 'both'):
        max_cov, min_cov = _max_nonzero_min(coverages)

    if mode in ('majority', 'both'):
        max_maj, min_maj = _max_nonzero_min(majorities)

    def format_percent(x, pos=None):
        return '%1.0f%%' % (100 * x)

    ax1.xaxis.set_tick_params(width=1.0)
    ax1.yaxis.set_tick_params(width=1.0)

    # if we're only doing coverage, include the number
    # of sequences, but with majority in the mix
    # it doesn't make sense
    if mode == 'majority':
        # y-ticks
        ax1.yaxis.set_major_formatter(FuncFormatter(format_percent))
        ax1.set_yticks(np.arange(0, 1.1, 0.2))
        ax1.set_ylim((0, 1))
        # alter axes to show max, min value
        ax1.spines['left'].set_bounds(min_maj, max_maj)
        # if we're not showing the # of sequences on the left,
        # remove the ticks and spine
        major_ticks += ax1.yaxis.get_major_ticks()
        ax1.spines['right'].set_visible(False)
        # x-ticks
        ax1.set_xticks(xticks)
        # set font properties
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels()
        )
    else:
        # y-ticks
        ax1.set_yticks(yticks)
        ax1.set_ylim((0, nseq))
        # alter axes to show max, min value
        ax1.spines['left'].set_bounds(min_cov, max_cov)
        # twin axis
        ax2 = ax1.twinx()
        ax2.xaxis.set_tick_params(width=1.0)
        ax2.yaxis.set_tick_params(width=1.0)
        # move the axis spines off the data
        _adjust_spines_outward(ax2, ('right',), 18)
        ax2.yaxis.set_major_formatter(FuncFormatter(format_percent))
        ax2.set_yticks(np.arange(0, 1.1, 0.2))
        # x-ticks
        ax2.set_xticks(xticks)
        # set transparent here otherwise ax2 doesn't exist
        if transparent:
            ax2.patch.set_alpha(0.)
        # get the major ticks so we can disable them later
        major_ticks += ax2.xaxis.get_major_ticks()
        # disable the top spines, like we do later
        ax2.spines['top'].set_visible(False)
        # set font properties
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels() +
            ax2.xaxis.get_ticklabels() +
            ax2.yaxis.get_ticklabels()
        )

    if mode == 'both':
        # use the axes spines to show the maximum value
        ax1.spines['right'].set_bounds(min_maj, max_maj)
        # set the colors of the axis spines to correspond to the data
        ax1.spines['left'].set_color(LIGHT_BLUE)
        ax1.spines['right'].set_color(LIGHT_RED)
    elif mode == 'coverage':
        # use the axes spines to show the maximum value
        ax1.spines['right'].set_bounds(min_cov, max_cov)

    for label in ticklabels:
        label.set_fontproperties(ROBOTO_REGULAR)

    if transparent:
        fig.patch.set_alpha(0.)
        ax1.patch.set_alpha(0.)

    # remove the upper ticks
    for tick in major_ticks:
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].set_visible(False)

    fig.savefig(
        filename,
        format=format,
        transparent=transparent,
        bbox_extra_artists=extra_artists,
        bbox_inches='tight',
        pad_inches=0.25
    )

    return filename


DNA_ALPHABET = Gapped(ambiguous_dna)
RNA_ALPHABET = Gapped(ambiguous_rna)
AMINO_ALPHABET = HasStopCodon(Gapped(extended_protein, gap_char=_GAP), stop_symbol=_STOP)


def graph_logo(
    alignment,
    columns,
    filename=None,
    dpi=None, edgecolor='k', figsize=None, format='pdf',
    labels=None, linewidth=0., transparent=True,
    refidx=None,
    cutoff=1e-2
):
    if filename is None:
        fd, filename = mkstemp()
        close(fd)

    if figsize is None:
        figsize = (3, 3)

    if labels is None:
        labels = ['%d' % (idx + 1) for idx in columns]

    if refidx is not None:
        msa = alignment
        alignment = msa[:refidx]
        alignment.extend(msa[refidx + 1:])
        del msa

    N = len(columns)

    alph = None
    for alph_ in (DNA_ALPHABET, RNA_ALPHABET, AMINO_ALPHABET):
        for r in alignment:
            r.seq.alphabet = alph_
        if all([_verify_alphabet(r.seq.upper()) for r in alignment]):
            alph = alph_
            break

    skips = _GAP
    if alph in (DNA_ALPHABET, RNA_ALPHABET):
        T = 'T' if alph == DNA_ALPHABET else 'U'
        letters = 'ACG' + T
        ambigs = {
            'M': 'AC',
            'R': 'AG',
            'W': 'A' + T,
            'S': 'CG',
            'Y': 'C' + T,
            'K': 'G' + T,
            'V': 'ACG',
            'H': 'AC' + T,
            'D': 'AG' + T,
            'B': 'CG' + T,
        }
        skips += 'N'
        colors = DNA_COLORS
    elif alph == AMINO_ALPHABET:
        letters = 'ACDEFGHIKLMNPQRSTVWY'
        ambigs = {
            'B': 'DN',
            'J': 'IL',
            'Z': 'EQ',
        }
        skips += _STOP + 'X' + 'OU'
        colors = AMINO_COLORS
    else:
        raise RuntimeError("sequences with indeterminable alphabet provided")

    s = len(letters)
    counts = np.zeros((s, N), dtype=float)

    for r in alignment:
        for j, c in enumerate(columns):
            ltr = r[c].upper()
            if ltr in skips:
                continue
            elif ltr in ambigs:
                frac = 1 / len(ambigs[ltr])
                for ltr_ in ambigs[ltr]:
                    i = letters.index(ltr_)
                    counts[i, j] += frac
            elif ltr in letters:
                i = letters.index(ltr)
                counts[i, j] += 1
            else:
                raise ValueError('unknown letter: {0}'.format(ltr))

    pwm = counts / counts.sum(axis=0)

    # compute the information content at each position
    maxbits = np.log2(s)
    e_n = (s - 1) / (2. * np.log(2) * counts.sum(axis=0))
    R = maxbits * np.ones((N,), dtype=float)
    R -= [-sum(v * np.log2(v) for v in pwm[:, j] if v) for j in range(N)]
    R -= e_n

    heights = np.zeros((s, N), dtype=float)
    identities = np.zeros((s, N), dtype=int)

    for j in range(N):
        ordered_ltrs = sorted(zip(range(s), pwm[:, j]), key=itemgetter(1))
        for (i, (k, v)) in enumerate(ordered_ltrs):
            h = R[j] * v
            heights[i, j] = h if h >= cutoff else 0
            identities[i, j] = k

    font = Basefont(join(FONT_PATH, 'Roboto-Black.ttf'))

    fig = plt.figure(figsize=figsize, dpi=dpi)

    # make each column a vertical golden rect
    rect = 0.2, 0.2, 0.382 * N, 0.618
    ax = fig.add_axes(rect)

    _adjust_spines_outward(ax, ('left',), 9)

    ax.set_ylabel('bits', fontproperties=ROBOTO_REGULAR)

    if figsize is None:
        fig.set_figwidth(N)

    if transparent:
        fig.patch.set_alpha(0.)
        ax.patch.set_alpha(0.)

    # remove the top and right ticks
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.tick2On = False

    # remove the bottom ticks
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1On = False

    # rotate the x-axis labels by 45 degrees to enhance packing
    for label in ax.xaxis.get_ticklabels():
        label.set_rotation(45)

    # set font properties
    for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
        label.set_fontproperties(ROBOTO_REGULAR)

    # disable top and right spines, we don't need them
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    def format_xlabel(x, pos=None):
        idx = np.clip(int(x) - 1, 0, N - 1)
        return labels[idx]

    ax.xaxis.set_major_formatter(FuncFormatter(format_xlabel))
    # avoid too much precision
    ax.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))

    # set the ticks
    ysep = 0.5 if s < 20 else 1.0
    yticks = np.arange(0, maxbits, ysep, dtype=float)
    if maxbits - yticks[-1] < ysep:
        yticks[-1] = maxbits
    else:
        yticks = np.append(yticks, maxbits)
    ax.set_yticks(yticks)
    ax.set_xticks(np.arange(1, N + 1, dtype=float) + 0.5)

    # set the axes limits here AFTER the ticks, otherwise borkage
    ax.set_xlim((1, N + 1))
    ax.set_ylim((0, maxbits))

    idxs = np.arange(1, N + 1)
    bottoms = np.zeros((N,), dtype=float)
    for i in range(s):
        bars = ax.bar(idxs, heights[i, :], width=1., bottom=bottoms)
        bottoms += heights[i, :]
        for j, bar in enumerate(bars):
            if heights[i, j]:
                l = letters[identities[i, j]]
                glyph = font[l]
                ax.add_patch(glyph)
                glyph.set_transform(bar.get_transform())
                glyph.set_edgecolor(edgecolor)
                glyph.set_facecolor(colors[l])
                glyph.set_linewidth(linewidth)
                glyph.set_zorder(-1)
            bar.set_visible(False)

    # set the remaining spine to show the maximum value
    ax.spines['left'].set_bounds(0, max(bottoms))

    # make these uniform
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)

    fig.savefig(
        filename, format=format, transparent=transparent,
        bbox_inches='tight', pad_inches=0.25
        )

    return filename

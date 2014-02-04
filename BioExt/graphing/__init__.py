
from __future__ import division, print_function

import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import repeat
from operator import itemgetter
from os import close
from os.path import dirname, join
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
from matplotlib.patches import Patch
from matplotlib.ticker import (
    FormatStrFormatter,
    FuncFormatter
    )

from BioExt.graphing._count import _count
from BioExt.graphing._basefont import Basefont
from BioExt.misc import _GAP, _STOP


__all__ = [
    'graph_coverage_majority',
    'graph_logo',
    'graph_readlength_histogram'
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


def _magic_ticks(lwr, upr, div=5):
    # some heuristic to prevent too many ticks,
    # is complicated by the start and end dynamically business
    # 0.477122 jumps an order 10 at 30% of the next order 10
    magic = 0.477122

#     div2 = 10 ** int(np.log10(upr - lwr) - magic) or 1
#     sep = div2 * int(upr / div / div2 + 1)
#     ticks = np.arange(lwr, upr + sep, sep)
#
#     if len(ticks) > 1:
#         if upr - ticks[-2] < sep / 4:
#             ticks = np.arange(lwr, upr, sep)
#
#     ticks[-1] = upr

    rng = upr - lwr
    div2 = 10 ** int(np.log10(rng) - magic) or 1
    sep = div2 * int(rng / div / div2 + 1)
    start = (lwr // sep) * sep
    end = upr + sep
    ticks = np.arange(start, end, sep)

    # prevent the first two and last two ticks from being too close
    if len(ticks) > 1:
        fix = False
        min_delta = sep / 4
        delta = ticks[1] - lwr
        if delta and delta < min_delta:
            start += sep
            fix = True
        delta = upr - ticks[-2]
        if delta and delta < min_delta:
            end -= sep
            fix = True
        if fix:
            ticks = np.arange(start, end, sep)

    ## first and last tick
    ticks[0] = lwr
    ticks[-1] = upr

    return ticks


def count_alignment(alignment, columns=None, refidx=None, limit=100):

    aln = iter(alignment)
    records = []

    for i, r in enumerate(aln, start=len(records)):
        if len(records) > limit:
            break
        if i == refidx:
            continue
        records.append((i, r))

    alph = None
    for alph_ in (DNA_ALPHABET, RNA_ALPHABET, AMINO_ALPHABET):
        for _, r in records:
            r.seq.alphabet = alph_
        if all(_verify_alphabet(r.seq.upper()) for _, r in records):
            alph = alph_
            break

    if alph is None:
        raise RuntimeError('unknown alphabet')

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

    def b(r):
        return r.upper()

    def allrecords():
        for i, r in records:
            yield b(str(r.seq))
        for i, r in enumerate(aln, start=i):
            if i == refidx:
                continue
            yield b(str(r.seq))

    alphabet = letters + ''.join(sorted(ambigs.keys()))
    nchar = len(alphabet)
    values = np.zeros((nchar, len(letters)), dtype=float)

    for i, c in enumerate(alphabet):
        if i < len(letters):
            values[i, i] = 1.
        else:
            v = 1. / len(ambigs[c])
            for d in ambigs[c]:
                values[i, alphabet.index(d)] = v

    if columns is not None:
        columns = np.array(columns, dtype=int)

    counts = _count(
        allrecords(),
        columns,
        b(alphabet),
        values
        )

    return counts.transpose(), (letters, colors)


def graph_coverage_majority(
    alignment,
    mode,
    nseq=None,
    filename=None,
    dpi=None, figsize=(6, 6), format='pdf', transparent=True,
    refidx=None
):
    if not mode:
        mode = 'both'

    if mode not in ('coverage', 'majority', 'both'):
        raise ValueError("mode must be one of 'coverage', 'majority', or 'both'")

    if filename is None:
        fd, filename = mkstemp()
        close(fd)

    if nseq is None:
        nseq = len(alignment)
        if refidx is not None:
            nseq -= 1

    counts, _ = count_alignment(alignment, refidx=refidx)
    fst_pos, lst_pos = 0, counts.shape[1]
    xs = np.arange(fst_pos, lst_pos) + 1

    # FIGURE
    fig = plt.figure(figsize=figsize, dpi=dpi)

    # golden rectangle!
    rect = 0.2, 0.2, 1, 0.618
    ax1 = fig.add_axes(rect)

    # move the axis spines off the data
    _adjust_spines_outward(ax1, ('left', 'right'), 3 * figsize[0])

    if transparent:
        fig.patch.set_alpha(0.)
        ax1.patch.set_alpha(0.)

    # PLOT
    if mode in ('majority', 'both'):
        if mode == 'majority':
            majorities = counts.max(axis=0) / counts.sum(axis=0)
        else:
            majorities = counts.max(axis=0)

        max_maj, min_maj = _max_nonzero_min(majorities)

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
        max_cov, min_cov = _max_nonzero_min(coverages)

        lines = ax1.plot(
            xs,
            coverages,
            color=LIGHT_BLUE,
            linewidth=1.0,
            zorder=-1
            )

        for l in lines:
            l.set_clip_on(False)

    # AXES
    if mode != 'majority':
        ax2 = ax1.twinx()
        _adjust_spines_outward(ax2, ('right',), 3 * figsize[0])
        if transparent:
            ax2.patch.set_alpha(0.)

    # TICKS
    ## compute the x-ticks
    xticks = _magic_ticks(fst_pos + 1, lst_pos)

    ## set the x-ticks
    ax1.set_xlim((fst_pos + 1, lst_pos))
    ax1.set_xticks(xticks)
    if mode != 'majority':
        ax2.set_xticks(xticks)

    ## get the major ticks
    major_ticks = ax1.xaxis.get_major_ticks(len(xticks))

    ## disable the first and last tick on the x-axis,
    ## they're redundant and ugly (esp if minH > 0)
    major_ticks[0].tick1On = False
    major_ticks[-1].tick1On = False

    ## compute the y-ticks
    yticks = _magic_ticks(0, nseq)

    ## set the y-ticks, limits, range
    if mode == 'majority':
        ### y-ticks
        ax1.set_ylim((0, 1))
        ax1.set_yticks(np.arange(0, 1.1, 0.2))
        # major ticks
        major_ticks += ax1.yaxis.get_major_ticks()
        # tick labels
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels()
        )
    else:
        ### y-ticks
        ax1.set_ylim((0, nseq))
        ax1.set_yticks(yticks)
        ax2.set_yticks(np.arange(0, 1.1, 0.2))
        # major ticks
        major_ticks += ax2.xaxis.get_major_ticks()
        # tick labels
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels() +
            ax2.xaxis.get_ticklabels() +
            ax2.yaxis.get_ticklabels()
        )

    ## remove the major ticks
    for tick in major_ticks:
        tick.tick2On = False

    ## set the tick label font
    for label in ticklabels:
        label.set_fontproperties(ROBOTO_REGULAR)

    ## width
    ax1.xaxis.set_tick_params(width=1.0)
    ax1.yaxis.set_tick_params(width=1.0)
    if mode != 'majority':
        ax2.xaxis.set_tick_params(width=1.0)
        ax2.yaxis.set_tick_params(width=1.0)

    # SPINES
    ## show max, min value, set visibility
    ax1.spines['top'].set_visible(False)
    if mode == 'majority':
        ax1.spines['left'].set_bounds(min_maj, max_maj)
        ax1.spines['right'].set_visible(False)
    else:
        ax1.spines['left'].set_bounds(min_cov, max_cov)
        ax2.spines['top'].set_visible(False)

    ## more bounds
    if mode == 'both':
        ax1.spines['right'].set_bounds(min_maj, max_maj)
    if mode == 'coverage':
        ax1.spines['right'].set_bounds(min_cov, max_cov)

    ## color
    if mode == 'both':
        ax1.spines['left'].set_color(LIGHT_BLUE)
        ax1.spines['right'].set_color(LIGHT_RED)

    # LABELS
    ## formatter
    def format_percent(x, pos=None):
        return '%1.0f%%' % (100 * x)

    if mode == 'majority':
        ax1.yaxis.set_major_formatter(FuncFormatter(format_percent))
    else:
        ax2.yaxis.set_major_formatter(FuncFormatter(format_percent))

    ax1.set_xlabel('Reference sequence position', fontproperties=ROBOTO_REGULAR)

    if mode == 'majority':
        ax1.set_ylabel('Majority proportion', fontproperties=ROBOTO_REGULAR)
    else:
        ax1.set_ylabel('No. of reads', fontproperties=ROBOTO_REGULAR)
        ax2.set_ylabel('Pct. of reads', rotation=270., fontproperties=ROBOTO_REGULAR)

    # LEGEND
    extra_artists = []

    if mode == 'both':
        # create a proxy artist for legend, PolyCollections don't work (heights)
        p1 = Line2D([0, 1], [0, 1], color=LIGHT_BLUE, linewidth=1.0)
        # create a proxy artist for legend, [Lines2D] don't work (majorities)
        p2 = Line2D([0, 1], [0, 1], color=LIGHT_RED, linewidth=1.0)
        leg = ax1.legend(
            [p1, p2],
            ['Coverage', 'Majority'],
            bbox_to_anchor=(0.5, -0.15),
            loc=9,
            ncol=2,
            prop=ROBOTO_REGULAR,
            borderpad=0
        )
        leg.legendPatch.set_alpha(0.)
        extra_artists.append(leg)

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
    dpi=None, edgecolor='k', figsize=(3, 3), format='pdf',
    labels=None, linewidth=0., transparent=True,
    refidx=None,
    cutoff=1e-2
):
    if filename is None:
        fd, filename = mkstemp()
        close(fd)

    if labels is None:
        labels = ['%d' % (idx + 1) for idx in columns]

    N = len(columns)
    counts, (letters, colors) = count_alignment(alignment, columns)
    pwm = counts / counts.sum(axis=0)
    s = len(letters)

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

    # FIGURE
    fig = plt.figure(figsize=figsize, dpi=dpi)

    # make each column a vertical golden rect
    rect = 0.2, 0.2, 0.382 * N, 0.618
    ax = fig.add_axes(rect)

    _adjust_spines_outward(ax, ('left',), 3 * figsize[0])

    if figsize is None:
        fig.set_figwidth(N)

    if transparent:
        fig.patch.set_alpha(0.)
        ax.patch.set_alpha(0.)

    # DRAW
    font = Basefont(join(FONT_PATH, 'Roboto-Black.ttf'))
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

    # TICKS
    ## set the ticks
    ysep = 0.5 if s < 20 else 1.0
    yticks = np.arange(0, maxbits, ysep, dtype=float)

    if maxbits - yticks[-1] < ysep:
        yticks[-1] = maxbits
    else:
        yticks = np.append(yticks, maxbits)

    ax.set_yticks(yticks)
    ax.set_xticks(np.arange(1, N + 1, dtype=float) + 0.5)

    ## set the axes limits here AFTER the ticks, otherwise borkage
    ax.set_xlim((1, N + 1))
    ax.set_ylim((0, maxbits))

    ## remove the top and right ticks
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.tick2On = False

    ## remove the bottom ticks
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1On = False

    ## rotate the x-axis tick labels by 45 degrees to enhance packing
    for label in ax.xaxis.get_ticklabels():
        label.set_rotation(45)

    ## set tick label font
    for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
        label.set_fontproperties(ROBOTO_REGULAR)

    ## make tick width uniform
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)

    # SPINES
    ## disable top and right spines, we don't need them
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ## set the remaining spine to show the maximum value
    ax.spines['left'].set_bounds(0, max(bottoms))

    # LABELS
    ## formatters
    def format_xlabel(x, pos=None):
        idx = np.clip(int(x) - 1, 0, N - 1)
        return labels[idx]

    ax.xaxis.set_major_formatter(FuncFormatter(format_xlabel))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))

    ax.set_ylabel('bits', fontproperties=ROBOTO_REGULAR)

    fig.savefig(
        filename, format=format, transparent=transparent,
        bbox_inches='tight', pad_inches=0.25
        )

    return filename


def graph_readlength_histogram(
        lengths,
        filename=None,
        bins=50, mean_median=False, rwidth=0.9,
        dpi=None, figsize=(6, 6), format='pdf', transparent=True
        ):

    if filename is None:
        fd, filename = mkstemp()
        close(fd)

    # FIGURE
    fig = plt.figure(figsize=figsize, dpi=dpi)

    rect = 0.2, 0.2, 1, 0.618
    ax = fig.add_axes(rect)

    _adjust_spines_outward(ax, ('left',), 3 * figsize[0])

    if transparent:
        fig.patch.set_alpha(0.)
        ax.patch.set_alpha(0.)

    n, bins, patches = ax.hist(
        lengths,
        bins=bins,
        rwidth=rwidth,
        color=LIGHT_BLUE,
        linewidth=0.
        )

    if mean_median:
        mean = np.mean(lengths)
        median = np.median(lengths)

        for i, upr in enumerate(bins[1:]):
            if upr > mean:
                patches[i].set_color(LIGHT_RED)
                mean = float('Inf')
            if upr > median:
                patches[i].set_color(LIGHT_GREEN)
                median = float('Inf')

    # TICKS
    xlwr, xupr = (int(v) for v in ax.get_xlim())
    ylwr, yupr = (int(v) for v in ax.get_ylim())

    xticks = _magic_ticks(xlwr, xupr)
    yticks = _magic_ticks(ylwr, yupr)

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    ax.set_xlim((xticks[0], xticks[-1]))
    ax.set_ylim((yticks[0], yticks[-1]))

    major_ticks = ax.xaxis.get_major_ticks()

    ## disable the first and last tick on the x-axis,
    ## they're redundant and ugly (esp if minH > 0)
    major_ticks[0].tick1On = False
    major_ticks[-1].tick1On = False

    ## remove the top and right ticks
    for tick in major_ticks + ax.yaxis.get_major_ticks():
        tick.tick2On = False

    ## set tick label font
    for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
        label.set_fontproperties(ROBOTO_REGULAR)

    ## make tick width uniform
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)

    # SPINES
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # LABELS
    ax.set_xlabel('Read length', fontproperties=ROBOTO_REGULAR)
    ax.set_ylabel('No. of reads', fontproperties=ROBOTO_REGULAR)

    # LEGEND
    extra_artists = []

    if mean_median:
        p1 = Patch(color=LIGHT_RED, linewidth=0.)
        p2 = Patch(color=LIGHT_GREEN, linewidth=0.)
        leg = ax.legend(
            [p1, p2],
            ['Mean', 'Median'],
            bbox_to_anchor=(0.5, -0.15),
            loc=9,
            ncol=2,
            prop=ROBOTO_REGULAR,
            borderpad=0
        )

        leg.legendPatch.set_alpha(0.)
        extra_artists.append(leg)

    fig.savefig(
        filename,
        format=format,
        transparent=transparent,
        bbox_extra_artists=extra_artists,
        bbox_inches='tight',
        pad_inches=0.25
    )

    return filename

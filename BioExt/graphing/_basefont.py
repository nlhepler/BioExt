
from __future__ import division, print_function

from os.path import exists

import numpy as np

from matplotlib.path import Path
from matplotlib.patches import PathPatch

from freetype import *


class Basefont(object):

    def __init__(self, filepath, charsize=48):
        if not exists(filepath):
            raise ValueError('No valid font file specified')

        try:
            if isinstance(filepath, unicode):
                filepath = filepath.encode('utf-8')
        except NameError:
            if not isinstance(filepath, bytes):
                filepath = filepath.encode('utf-8')

        face = Face(filepath)
        face.set_char_size(charsize * 64)
        self.face = face
        self.cache = {}

    def __getitem__(self, value):
        return Basefont.char_patch(self, value)

    def char_patch(self, char):
        # if we have the path in our cache, re-use it
        if char in self.cache:
            return PathPatch(self.cache[char])

        self.face.load_char(char)
        slot = self.face.glyph
        outline = slot.outline

        points_ = np.array(outline.points, dtype=[('x', float), ('y', float)])
        xs, ys = points_['x'], points_['y']

        xmin, ymin = xs.min(), ys.min()
        dx, dy = (xs.max() - xmin), (ys.max() - ymin)

        start, end = 0, 0
        verts, codes = [], []

        for i in range(len(outline.contours)):
            end = outline.contours[i]
            points = outline.points[start:end+1]
            points.append(points[0])
            # normalize points to (0, 1) space
            points = [((x - xmin) / dx, (y - ymin) / dy) for x, y in points]
            tags = outline.tags[start:end+1]
            tags.append(tags[0])

            verts.append(points[0])
            codes.append(Path.MOVETO)

            segments = [[points[0]]]
            for j in range(1, len(points)):
                segments[-1].append(points[j])
                if tags[j] & (1 << 0) and j < (len(points) - 1):
                    segments.append([points[j]])

            for segment in segments:
                if len(segment) == 2:
                    verts.extend(segment[1:])
                    codes.append(Path.LINETO)
                elif len(segment) == 3:
                    verts.extend(segments[1:])
                    codes.extend([Path.CURVE3, Path.CURVE3])
                else:
                    verts.append(segment[1])
                    codes.append(Path.CURVE3)
                    for i in range(1, len(segment)-2):
                        A, B = segment[i], segment[i + 1]
                        C = ((A[0] + B[0]) / 2., (A[1] + B[1]) / 2.)
                        verts.extend([C, B])
                        codes.extend([Path.CURVE3, Path.CURVE3])
                    verts.append(segment[-1])
                    codes.append(Path.CURVE3)

            start = end + 1

        path = Path(verts, codes)
        self.cache[char] = path
        glyph = PathPatch(path, antialiased=True)

        return glyph

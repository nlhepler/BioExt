#!/usr/bin/env python

from __future__ import division, print_function

import os.path
import numpy

from setuptools import Extension, setup

from BioExt import __version__ as _bioext_version
from BioExt.references._factory import _installrefdirs


np_inc = [os.path.join(os.path.dirname(numpy.__file__), 'core', 'include')]

ext_modules = [
    Extension(
        'BioExt.align._align',
        sources=[
            os.path.join('BioExt', 'align', '_align.c'),
            os.path.join('BioExt', 'align', 'alignment.c')
            ],
        include_dirs=np_inc,
        libraries=['m'],
        extra_compile_args=['-O3']
        ),
    Extension(
        'BioExt.graphing._count',
        sources=[os.path.join('BioExt', 'graphing', '_count.c')],
        include_dirs=np_inc,
        extra_compile_args=['-O3']
        ),
    Extension(
        'BioExt.merge._merge',
        sources=[
            os.path.join('BioExt', 'merge', '_merge.c'),
            os.path.join('BioExt', 'merge', 'merge.cpp')
            ],
        extra_compile_args=['-O3']
        ),
    Extension(
        'BioExt.rateclass._rateclass',
        sources=[
            os.path.join('BioExt', 'rateclass', '_rateclass.cpp'),
            os.path.join('BioExt', 'rateclass', 'rateclass.cpp')
            ],
        extra_compile_args=['-O3']
        )
    ]


setup(
    name='bioext',
    version=_bioext_version,
    description='Misc utilities and definitions not included or hidden in BioPython',
    author='N Lance Hepler',
    author_email='nlhepler@gmail.com',
    url='http://github.com/nlhepler/bioext',
    license='GNU GPL version 3',
    packages=[
        'BioExt',
        'BioExt.align',
        'BioExt.args',
        'BioExt.collections',
        'BioExt.errorize',
        'BioExt.freetype',
        'BioExt.freetype.ft_enums',
        'BioExt.graphing',
        'BioExt.io',
        'BioExt.io.BamIO',
        'BioExt.io.LazyAlignIO',
        'BioExt.io.SamIO',
        'BioExt.joblib',
        'BioExt.joblib.test',
        'BioExt.merge',
        'BioExt.misc',
        'BioExt.ndarray',
        'BioExt.optimize',
        'BioExt.orflist',
        'BioExt.phylo',
        'BioExt.quiver',
        'BioExt.rateclass',
        'BioExt.references',
        'BioExt.scorematrices',
        'BioExt.stats',
        'BioExt.uds',
        'BioExt.untranslate'
        ],
    package_dir={
        'BioExt': 'BioExt',
        'BioExt.align': 'BioExt/align',
        'BioExt.args': 'BioExt/args',
        'BioExt.collections': 'BioExt/collections',
        'BioExt.errorize': 'BioExt/errorize',
        'BioExt.freetype': 'BioExt/freetype',
        'BioExt.freetype.ft_enums': 'BioExt/freetype/ft_enums',
        'BioExt.graphing': 'BioExt/graphing',
        'BioExt.io': 'BioExt/io',
        'BioExt.io.BamIO': 'BioExt/io/BamIO',
        'BioExt.io.LazyAlignIO': 'BioExt/io/LazyAlignIO',
        'BioExt.io.SamIO': 'BioExt/io/SamIO',
        'BioExt.joblib': 'BioExt/joblib',
        'BioExt.joblib.test': 'BioExt/joblib/test',
        'BioExt.merge': 'BioExt/merge',
        'BioExt.misc': 'BioExt/misc',
        'BioExt.ndarray': 'BioExt/ndarray',
        'BioExt.optimize': 'BioExt/optimize',
        'BioExt.orflist': 'BioExt/orflist',
        'BioExt.phylo': 'BioExt/phylo',
        'BioExt.quiver': 'BioExt/quiver',
        'BioExt.rateclass': 'BioExt/rateclass',
        'BioExt.references': 'BioExt/references',
        'BioExt.scorematrices': 'BioExt/scorematrices',
        'BioExt.stats': 'BioExt/stats',
        'BioExt.uds': 'BioExt/uds',
        'BioExt.untranslate': 'BioExt/untranslate'
        },
    package_data={
        'BioExt': [
            'data/fonts/ttf/*.ttf',
            'data/scorematrices/*.txt'
            ] + _installrefdirs
        },
    scripts=[
        'scripts/bam2fna',
        'scripts/bam2msa',
        'scripts/bamclip',
        'scripts/bealign',
        'scripts/begraph',
        'scripts/clipedge',
        'scripts/consensus',
        'scripts/msa2bam',
        'scripts/seqmerge',
        'scripts/translate'
        # 'scripts/variants'
        ],
    ext_modules=ext_modules,
    requires=[
        'Bio (>=1.58)',
        'matplotlib (>=1.2.1)',
        'numpy (>=1.6)',
        'pysam (>=0.7.5)',
        ]
    )

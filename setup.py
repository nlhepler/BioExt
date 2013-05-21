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
        'BioExt.aligner._align',
        sources=[
            os.path.join('BioExt', 'aligner', '_align.c'),
            os.path.join('BioExt', 'aligner', 'alignment.c')
            ],
        include_dirs=np_inc,
        libraries=['m'],
        extra_compile_args=['-O3']
        ),
    Extension(
        'BioExt.merger._merge',
        sources=[
            os.path.join('BioExt', 'merger', '_merge.c'),
            os.path.join('BioExt', 'merger', 'merge.cpp')
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
        'BioExt.aligner',
        'BioExt.collections',
        'BioExt.errorize',
        'BioExt.io',
        'BioExt.io.BamIO',
        'BioExt.io.LazyAlignIO',
        'BioExt.io.SamIO',
        'BioExt.joblib',
        'BioExt.joblib.test',
        'BioExt.merger',
        'BioExt.misc',
        'BioExt.ndarray',
        'BioExt.orflist',
        'BioExt.phylo',
        'BioExt.quiver',
        'BioExt.rateclass',
        'BioExt.references',
        'BioExt.scorematrix',
        'BioExt.stats',
        'BioExt.uds',
        'BioExt.untranslate'
        ],
    package_dir={
        'BioExt': 'BioExt',
        'BioExt.aligner': 'BioExt/aligner',
        'BioExt.collections': 'BioExt/collections',
        'BioExt.errorize': 'BioExt/errorize',
        'BioExt.io': 'BioExt/io',
        'BioExt.io.BamIO': 'BioExt/io/BamIO',
        'BioExt.io.LazyAlignIO': 'BioExt/io/LazyAlignIO',
        'BioExt.io.SamIO': 'BioExt/io/SamIO',
        'BioExt.joblib': 'BioExt/joblib',
        'BioExt.joblib.test': 'BioExt/joblib/test',
        'BioExt.merger': 'BioExt/merger',
        'BioExt.misc': 'BioExt/misc',
        'BioExt.ndarray': 'BioExt/ndarray',
        'BioExt.orflist': 'BioExt/orflist',
        'BioExt.phylo': 'BioExt/phylo',
        'BioExt.quiver': 'BioExt/quiver',
        'BioExt.rateclass': 'BioExt/rateclass',
        'BioExt.references': 'BioExt/references',
        'BioExt.scorematrix': 'BioExt/scorematrix',
        'BioExt.stats': 'BioExt/stats',
        'BioExt.uds': 'BioExt/uds',
        'BioExt.untranslate': 'BioExt/untranslate'
        },
    package_data={
        'BioExt': [
            'data/scorematrices/*.txt'
            ] + _installrefdirs
        },
    scripts=[
        'scripts/bam2fna',
        'scripts/bam2msa',
        'scripts/bamclip',
        'scripts/bealign',
        # 'scripts/consensus',
        # 'scripts/variants'
        ],
    ext_modules=ext_modules,
    requires=[
        'Bio (>=1.58)',
        'numpy (>=1.6)',
        ]
    )

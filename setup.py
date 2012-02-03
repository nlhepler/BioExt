#!/usr/bin/env python

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from BioExt import __version__ as _bioext_version

setup(name='bioext',
      version=_bioext_version,
      description='Misc utilities and definitions not included or hidden in biopython',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/nlhepler/bioext',
      license='GNU GPL version 3',
      packages=['BioExt'],
      package_dir={'BioExt': 'lib/BioExt'}
     )

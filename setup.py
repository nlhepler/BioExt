#!/usr/bin/env python

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from BioExt import __version__ as _bioext_version
from BioExt._references._factory import _installrefdirs

setup(name='bioext',
      version=_bioext_version,
      description='Misc utilities and definitions not included or hidden in biopython',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/nlhepler/bioext',
      license='GNU GPL version 3',
      packages=['BioExt', 'BioExt._references'],
      package_dir={
        'BioExt': 'lib/BioExt',
        'BioExt._references': 'lib/BioExt/_references'
      },
      package_data={
        'BioExt': _installrefdirs
      }
     )

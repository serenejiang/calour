#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
from setuptools import find_packages, setup


# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('calour/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'CALOUR: exploratory and interactive microbiome analyses based on heatmap'

with open('README.md') as f:
    long_description = f.read()

keywords = 'microbiome heatmap analysis bioinformatics',

setup(name='calour',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="calour development team",
      maintainer="calour development team",
      url='http://microbio.me/calour',
      test_suite='nose.collector',
      packages=find_packages(),
      package_data={'calour': ['log.cfg', 'calour.config']},
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'matplotlib >= 2.0',
          'scikit-learn',
          'scikit-bio',
          'biom-format',
          'statsmodels',
          'ipywidgets'],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"],
                      'doc': ["Sphinx >= 1.4"]},
      entry_points={
          'console_scripts': [
              'calour=calour.cli:cmd',
          ]})

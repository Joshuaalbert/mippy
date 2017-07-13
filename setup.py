#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

import mippy

setup_requires = ['numpy>=' + mippy.__minimum_numpy_version__]

setup(name='mippy',
      version='0.0.1',
      description='Mixed integer programming interface from python',
      author=['Josh Albert'],
      author_email=['albert@strw.leidenuniv.nl'],
##      url='https://www.python.org/sigs/distutils-sig/',
    setup_requires=setup_requires,  
    tests_require=[
        'pytest>=2.8',
    ],
      packages=['mippy']
     )


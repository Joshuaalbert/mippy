#!/usr/bin/env python

from setuptools import setup
from setuptools import Extension
import numpy as np

from numpy.distutils.misc_util import get_info

#https://docs.python.org/3/distutils/apiref.html#module-distutils.core
lpsolver_interface_module = Extension("_lpsolver_interface", 
sources=['lpsolver_interface.c'],
include_dirs=['../../../local/include',np.get_include()], 
libraries=['m'],
extra_objects = ['../../../local/lib/libglpk.a'],
extra_compile_args=['--std=c99'])

setup(name='lpsolver_interface',
    ext_modules = [lpsolver_interface_module]
     )


#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.extension import Extension

import os
from setuptools import setup
from setuptools.command.install import install
import subprocess

__minimum_numpy_version__ = '1.9.0'

import numpy as np

#https://docs.python.org/3/distutils/apiref.html#module-distutils.core
lpsolver_interface_module = Extension("mippy.interface.C_lpsolver_interface", 
sources=['src/mippy/interface/C_lpsolver_interface.c'],
include_dirs=['local/include',np.get_include()], 
libraries=['m'],
extra_objects = ['local/lib/libglpk.a'],
extra_compile_args=['--std=c99'])

def compile_and_install_software():
    """Use the subprocess module to compile/install the C software."""
    src_path = '.'

    # install the software (into the virtualenv bin dir if present)
    subprocess.check_call('./install-glpk.sh', cwd=src_path, shell=True)

    #src_path = './src/mippy/interface'

    # install the software (into the virtualenv bin dir if present)
    #subprocess.check_call('make', cwd=src_path, shell=True)


class LPSolverInstall(install, object):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_install_software()
        super(LPSolverInstall,self).run()


setup_requires = ['numpy>=' + __minimum_numpy_version__]

setup(name='mippy',
      version='0.0.1',
      description='Mixed integer programming interface from python',
      author=['Josh Albert'],
      author_email=['albert@strw.leidenuniv.nl'],
##      url='https://www.python.org/sigs/distutils-sig/',
    setup_requires=setup_requires,  
    tests_require=[
        'pytest>=2.8', 'astropy'
    ],
    package_dir = {'':'src'},
  #  package_data={'mippy/interface':['lpsolver_interface']},
      packages=find_packages('src'),#['mippy','mippy/interface'],
    cmdclass={'install': LPSolverInstall},
    ext_modules = [lpsolver_interface_module]
     )


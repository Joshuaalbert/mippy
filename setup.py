#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

import mippy

import os
from setuptools import setup
from setuptools.command.install import install
import subprocess


def compile_and_install_software():
    """Use the subprocess module to compile/install the C software."""
    src_path = '.'

    # install the software (into the virtualenv bin dir if present)
    subprocess.check_call('./install-glpk.sh', cwd=src_path, shell=True)

    src_path = './src/mippy/interface'

    # install the software (into the virtualenv bin dir if present)
    subprocess.check_call('make', cwd=src_path, shell=True)


class LPSolverInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_install_software()
        super(LPSolverInstall,self).run()


setup_requires = ['numpy>=' + mippy.__minimum_numpy_version__, 'astropy']

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
    package_dir = {'':'src'},
    package_data={'mippy/interface':['lpsolver_interface']},
      packages=find_packages('src'),#['mippy','mippy/interface'],
    cmdclass={'install': LPSolverInstall}
     )


# mippy
Mixed Integer Programming Solver interface from Python

### Install (Linux)
```
git clone https://github.com/Joshuaalbert/mipp.git
cd mippy
./install-glpk.sh
python setup.py install
export PATH=$(pwd)/src/mippy/interface:${PATH}

#test with
pytest
```
### Status

[![Build Status](https://travis-ci.org/Joshuaalbert/mippy.svg?branch=master)](https://travis-ci.org/Joshuaalbert/mippy)

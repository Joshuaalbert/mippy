#!/bin/bash

# replace with yuor preferred version of glpk
wget -N -nc http://ftp.gnu.org/gnu/glpk/glpk-4.62.tar.gz

#leave thigns under here alone

shopt -s nullglob
array=(glpk*.tar.gz)

newest=${array[-1]}
echo With install ${newest}

dir=${newest%.tar.gz}
tar xzf ${newest}

root=$(pwd)
cd ${dir} && ./configure --prefix=${root}/local
make && make install

#done in setup.py
#echo Compiling interface
#cd ${root}/src/mippy/interface && make


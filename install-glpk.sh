#!/bin/bash
shopt -s nullglob
array=(glpk*.tar.gz)

newest=${array[-1]}
echo With install ${newest}

dir=${newest%.tar.gz}
tar xzf ${newest}

root=$(pwd)
cd ${dir} && ./configure --prefix=${root}/local
make && make install

echo Compiling interface
cd ${root}/mippy/interface && make


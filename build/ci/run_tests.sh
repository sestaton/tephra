#!/bin/bash

dir=`pwd`
git clone https://github.com/bioperl/bioperl-live.git
cd bioperl-live
echo "n" | perl Build.PL
./Build install
cd $dir

perl Makefile.PL && make test



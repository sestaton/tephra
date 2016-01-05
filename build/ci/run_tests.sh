#!/bin/bash

dir=`pwd`
git clone https://github.com/bioperl/bioperl-live.git
cd bioperl-live
echo "n" | perl Build.PL 2>&1 > /dev/null
./Build install 2>&1 > /dev/null
cd $dir

perl Makefile.PL && make && prove -bv t/01-findltrs.t



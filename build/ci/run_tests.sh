#!/bin/bash

set -euo pipefail

##
#dir=`pwd`

#cp build/ci/tephra-deps-ubuntu-precise.tar.bz2 ~/
#cd
#tar xjf tephra-deps-ubuntu-precise.tar.bz2
#cd $dir

perl Makefile.PL --debug
make 
#prove -bv t/04-findtirs.t
#prove -bv t/0[45]*t
make test
#prove -bv t/01-findltrs.t
#ls -l t/test_data
#perl -Mblib blib/bin/tephra findltrs -c t/test_data/tephra_ltr_config.yml -g t/test_data/ref.fas --debug


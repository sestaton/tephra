#!/bin/bash

set -euo pipefail

dir=`pwd`

cp build/ci/tephra-deps-ubuntu-precise.tar.bz2 ~/
cd
tar xjf tephra-deps-ubuntu-precise.tar.bz2
cd $dir

wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar xzf hmmer-3.1b2-linux-intel-x86_64.tar.gz
sudo cp hmmer-3.1b2-linux-intel-x86_64/binaries/* /usr/local/bin/

perl Makefile.PL
make
prove -bv t/04-findtirs.t
#prove -bv t/0[123]*t
#make test
#perl -Mblib blib/bin/tephra findltrs -c t/test_data/tephra_ltr_config.yml -g t/test_data/ref.fas -t t/test_data/trnas.fas -d t/test_data/te.hmm --clean
#prove -bv t/01-findltrs.t



#!/bin/bash

#dir=`pwd`
#git clone https://github.com/bioperl/bioperl-live.git
#cd bioperl-live
#echo "n" | perl Build.PL 2>&1 > /dev/null
#./Build install 2>&1 > /dev/null
#cd $dir

cp build/ci/tephra-deps-ubuntu-precise.tar.bz2 ~/
cd
tar xjf tephra-deps-ubuntu-precise.tar.bz2
cd $dir

#echo "Contents of home: "
#ls -la ~/
#ls -l ~/.tephra

echo "ltrharvest:
  - mintsd: 4
  - maxtsd: 6
  - minlenltr: 100
  - maxlenltr: 6000
  - mindistltr: 1500
  - maxdistltr: 25000
  - seedlength: 30
  - tsdradius: 60
  - xdrop: 5
  - swmat: 2 
  - swmis: -2
  - swins: -3
  - swdel: -3
  - overlaps: best
ltrdigest:
  - pptradius: 30
  - pptlen: 8 30
  - pptagpr: 0.25
  - uboxlen: 3 30
  - uboxutpr: 0.91
  - pbsradius: 30
  - pbslen: 11 30
  - pbsoffset: 0 5
  - pbstrnaoffset: 0 5
  - pbsmaxeditdist: 1
  - pdomevalue: 10E-6
  - pdomcutoff: NONE
  - maxgaplen: 50" > t/test_data/tephra_ltr_config.yml

perl Makefile.PL
make
perl -Mblib blib/bin/tephra findltrs -c t/test_data/tephra_ltr_config.yml -g t/test_data/ref.fas -t t/test_data/trnas.fas -d t/test_data/te.hmm --clean
#prove -bv t/01-findltrs.t



#!/bin/bash

if [ -d data ]; then
    mkdir -p src
elif [ ! -d data ]; then
    mkdir -p data src
fi

vers=`perl -Ilib bin/chloro --version | sed 's/chloro //;s/[()]//g;s/ bin\/chloro//'`
echo -e "========== Getting dependencies for: $vers ==========\n"

cd src

## Velvet
curl -L https://api.github.com/repos/dzerbino/velvet/tarball > velvet.tar.gz
tar xzf velvet.tar.gz
mv dzerbino-velvet* velvet && cd velvet
make MAXKMERLENGTH=99 2> /dev/null
cd ..

## VelvetOptimiser
curl -L https://api.github.com/repos/sestaton/VelvetOptimiser/tarball > VO.tar.gz
tar xzf VO.tar.gz
mv sestaton-VelvetOptimiser* VelvetOptimiser
cd ..

## Pairfq-lite
cd bin
curl -L git.io/pairfq_lite > pairfq_lite
chmod +x pairfq_lite

echo -e "\n========== Done with dependencies, now building package: $vers =========="
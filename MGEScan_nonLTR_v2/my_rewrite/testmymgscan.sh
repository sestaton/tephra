#!/bin/bash

set -euo pipefail

export HMMER2=/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2

dir=`pwd`
genome=genome
pdir=$dir
data=data

perl run_mgescan -d $data -p $pdir -g $genome

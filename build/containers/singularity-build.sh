#!/bin/bash

set -euo pipefail

vers=$(egrep "our.*VERSION" bin/tephra | sed "s/^.* '//;s/'.*$//")

export SINGULARITYENV_LC_ALL="C"
export SINGULARITY_LOCALCACHEDIR=${PWD}/image
export SINGULARITY_TMPDIR=${PWD}/image
export TMPDIR=${PWD}/image

mkdir -p $TMPDIR

# build singularity image
docker run -v \
    /var/run/docker.sock:/var/run/docker.sock \
    -v $TMPDIR:/output \
    --privileged \
    -t \
    --rm \
    singularityware/docker2singularity \
    --name tephra \
    sestaton/tephra

# sign image
singularity sign image/tephra.simg

# push to remote
echo -e "\n=====> Pushing tephra:${vers} and tephra:latest to Sylabs Cloud\n"
singularity push image/tephra.simg library://sestaton/default/tephra:${vers}
singularity push image/tephra.simg library://sestaton/default/tephra:latest

# clean up -- this is important because we do not want to commit any images
rm -rf $TMPDIR

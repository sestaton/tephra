#!/bin/bash

set -euo pipefail

vers=$(egrep "our.*VERSION" bin/tephra | sed "s/^.* '//;s/'.*$//")

export SINGULARITYENV_LC_ALL="C"
export SINGULARITY_LOCALCACHEDIR=${PWD}/image
export SINGULARITY_TMPDIR=${PWD}/image
export TMPDIR=${PWD}/image

mkdir -p $TMPDIR

# build singularity image
echo -e "\n=====> [TEPHRA]: Pulling Docking image and building Singularity image now.\n"
docker run -v \
    /var/run/docker.sock:/var/run/docker.sock \
    -v $TMPDIR:/output \
    --privileged \
    -t \
    --rm \
    quay.io/singularity/docker2singularity \
    --name tephra \
    sestaton/tephra

# -- sign image --
# NB: This is likely to fail because keys expire on Sylabs Cloud after only a short time. Therefore we warn....
echo -e "\n=====> [TEPHRA] Attempting to sign image now. If this fails, log on to: https//cloud.sylabs.io/ to update keys or sign image\n"
singularity sign image/tephra.sif

# push to remote
echo -e "\n=====> [TEPHRA] Pushing tephra:${vers} and tephra:latest to Sylabs Cloud\n"
singularity push image/tephra.sif library://sestaton/default/tephra:${vers}
singularity push image/tephra.sif library://sestaton/default/tephra:latest

# clean up -- this is important because we do not want to commit any images
rm -rf $TMPDIR

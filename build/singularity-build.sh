#!/bin/bash

set -euo pipefail

vers=$(egrep "our.*VERSION" bin/tephra | sed "s/^.* '//;s/'.*$//")

export SINGULARITYENV_LC_ALL=C

mkdir -p image

# build singularity image
for ver in "latest";
do
    sing_vers=$(echo $ver | sed 's/\.[0-9]$//')
    # this is certainly a bug in Singularity to not support semantic versioning

    docker run -v \
	/var/run/docker.sock:/var/run/docker.sock \
	-v ${PWD}/image:/output \
	--privileged \
	-t \
	--rm \
	singularityware/docker2singularity \
	--name tephra:$ver \
	sestaton/tephra

    # sign image
    singularity sign image/tephra:${sing_vers}.simg

    # push to remote
    singularity push image/tephra:${sing_vers}.simg library://sestaton/default/tephra
done

# clean up -- this is important because we do not want to commit any images
rm -rf image/
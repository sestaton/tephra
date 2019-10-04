#!/bin/bash

set -euo pipefail

vers=$(egrep "our.*VERSION" bin/tephra | sed "s/^.* '//;s/'.*$//")
#echo $vers

#--build-arg LC_ALL=C
# build
docker build \
-t sestaton/tephra:$vers .

# tag
docker tag sestaton/tephra:$vers sestaton/tephra:latest

# push
docker push sestaton/tephra:$vers
docker push sestaton/tephra:latest

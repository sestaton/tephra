#!/bin/bash

set -euo pipefail

vers=$(egrep "our.*VERSION" bin/tephra | sed "s/^.* '//;s/'.*$//")
echo "=====> Building Docker image for Tephra v$vers"

#--build-arg LC_ALL=C
# build
docker build \
-t sestaton/tephra:$vers .

echo "=====> Tagging Docker image for Tephra v$vers"
# tag
docker tag sestaton/tephra:$vers sestaton/tephra:latest

echo "=====> Pushing Docker image for Tephra v$vers to Docker Hub"
# push
docker push sestaton/tephra:$vers
docker push sestaton/tephra:latest

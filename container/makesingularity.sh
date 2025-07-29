#!/bin/sh

# Builds a Singularity image from the Docker image.
# Use makedocker.sh first.
# See https://stackoverflow.com/a/60316979

TAG="1.0.3"
REPO="nrlabcruk/invar2:$TAG"
IMAGE="invar2-${TAG}.sif"

sudo rm -rf invar_sandbox $IMAGE

#sudo singularity build --sandbox invar_sandbox docker-daemon://${REPO}
sudo singularity build $IMAGE docker-daemon://${REPO}

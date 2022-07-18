#!/bin/sh

# Builds a Singularity image from the Docker image.
# Use makedocker.sh first.
# See https://stackoverflow.com/a/60316979

TAG="latest"
REPO="nrlabcruk/invar2:$TAG"

sudo rm -rf invar_sandbox invar.sif

#sudo singularity build --sandbox invar_sandbox docker-daemon://${REPO}
sudo singularity build invar.sif docker-daemon://${REPO}


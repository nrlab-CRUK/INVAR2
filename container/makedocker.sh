#!/bin/sh

TAG="latest"
REPO="crukcibioinformatics/invar2:$TAG"

sudo docker build --tag "$REPO" --file Dockerfile .
sudo docker push "$REPO"


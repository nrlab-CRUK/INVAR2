#!/bin/sh

TAG="1.0.3"
REPO="nrlabcruk/invar2:$TAG"

sudo docker build --tag "$REPO" --file Dockerfile .
# sudo docker push "$REPO"

FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG CONDA_VERSION=py39_4.10.3

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -o miniconda3.sh && \
    mkdir -p /opt && \
    sh miniconda3.sh -b -p /opt/conda && \
    rm miniconda3.sh

COPY invar2.yml .
RUN /opt/conda/bin/conda env create -f invar2.yml

RUN /opt/conda/bin/conda clean -a

ENV PATH /opt/conda/envs/invar2/bin:$PATH


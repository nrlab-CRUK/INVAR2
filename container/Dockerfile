FROM    centos:8

LABEL   maintainer="Richard Bowers<richard.bowers@cruk.cam.ac.uk>"

ARG INSTALL_DIR=/usr/local
ARG BUILD_DIR=/var/tmp/invar_software_build

ARG BEDTOOLS=2.27.1
ARG PICARD=2.10.2
ARG HTSLIB=1.9
ARG SAMTOOLS=1.9
ARG BCFTOOLS=1.6

ARG TAROPTS="--no-same-owner --no-same-permissions"

ENV JAVA_HOME=/usr/lib/jvm/jre-11
ENV PATH=$JAVA_HOME/bin:$PATH

RUN dnf install -y dnf-plugins-core epel-release
RUN dnf config-manager --set-enabled powertools
RUN dnf makecache && dnf update -y
RUN dnf install -y \
	bzip2 bzip2-devel gcc-c++ git java-11-openjdk-headless \
	libcurl-devel libxml2-devel make ncurses-devel openssl-devel \
	python2 python3 R-core R-devel unzip wget xz-devel zlib-devel

RUN alternatives --set python /usr/bin/python2
RUN mkdir -p ${INSTALL_DIR} ${BUILD_DIR}

# htslib
# RUN cd ${BUILD_DIR}; \
#     wget https://github.com/samtools/htslib/releases/download/${HTSLIB}/htslib-${HTSLIB}.tar.bz2 && \
#     tar ${TAROPTS} -jxf htslib-${HTSLIB}.tar.bz2
# RUN cd ${BUILD_DIR}/htslib-${HTSLIB}; \
#     ./configure --prefix="${INSTALL_DIR}" && \
#     make install

# Can add "--with-htslib=system" to use a shared htslib in samtools & bcftools, if the same version.

# samtools
RUN cd ${BUILD_DIR}; \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS}/samtools-${SAMTOOLS}.tar.bz2 && \
    tar ${TAROPTS} -jxf samtools-${SAMTOOLS}.tar.bz2
RUN cd ${BUILD_DIR}/samtools-${SAMTOOLS}; \
    ./configure --prefix="${INSTALL_DIR}" && \
    make && \
    make install

# bcftools
RUN cd ${BUILD_DIR}; \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS}/bcftools-${BCFTOOLS}.tar.bz2 && \
    tar ${TAROPTS} -jxf bcftools-${BCFTOOLS}.tar.bz2
RUN cd ${BUILD_DIR}/bcftools-${BCFTOOLS}; \
    ./configure --prefix="${INSTALL_DIR}" && \
    make && \
    make install

# bedtools2
RUN cd ${BUILD_DIR}; \
    wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS}/bedtools-${BEDTOOLS}.tar.gz && \
    tar ${TAROPTS} -xzf ${BUILD_DIR}/bedtools-${BEDTOOLS}.tar.gz
RUN cd ${BUILD_DIR}/bedtools2; \
    make install prefix="" DESTDIR="${INSTALL_DIR}"

# Picard
RUN cd ${INSTALL_DIR}; \
    wget -O picard-${PICARD}.jar https://github.com/broadinstitute/picard/releases/download/${PICARD}/picard.jar

# R packages
RUN /usr/bin/R --vanilla -e \
    "install.packages(pkgs = c('optparse', 'tidyverse'), repos = c('https://cran.ma.imperial.ac.uk/', 'https://www.stats.bris.ac.uk/R/'))"

## Clean up
RUN cd / && rm -rf ${BUILD_DIR}

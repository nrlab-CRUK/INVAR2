#!/bin/sh

rm -r .nextflow* work results analysis
nextflow run testing.nf -profile bigserver

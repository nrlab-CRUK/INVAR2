#!/bin/bash

set -o pipefail

# generate slopped bed file
bedtools slop -i "!{bedFile}" -g "!{genomeFile}" -b !{params.SLOP_BP} > slopped.bed

# Filter lines that are <= 100bp in range.
Rscript --vanilla "!{projectDir}/R/parse/filterSloppedBed.R" slopped.bed "!{filteredBedFile}"

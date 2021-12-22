#!/bin/bash -l

set -euo pipefail

# generate slopped bed file
bedtools slop -i "!{bedFile}" -g "!{genomeFile}" -b !{params.SLOP_BP} > slopped.bed

# Filter lines that are <= 100bp in range.
Rscript --vanilla "!{projectDir}/R/invar1/filterSloppedBed.R" slopped.bed "!{filteredBedFile}"

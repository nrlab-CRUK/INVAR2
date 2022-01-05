#!/bin/bash

set -euo pipefail

# Convert patient CSV to mutation list in BED format.
Rscript --vanilla "!{projectDir}/R/invar1/patientListCsvToBed.R" \
    "!{csvFile}" mutationlist.bed

# generate slopped bed file
bedtools slop -i mutationlist.bed -g "!{genomeFile}" -b !{params.SLOP_BP} > slopped.bed

# Filter lines that are <= 100bp in range.
Rscript --vanilla "!{projectDir}/R/invar1/filterSloppedBed.R" slopped.bed "!{filteredBedFile}"

#!/bin/bash

set -euo pipefail

# Convert patient CSV to mutation list in BED format.
Rscript --vanilla "!{params.projectHome}/R/1_parse/patientListCsvToBed.R" \
    "!{csvFile}" mutationlist.bed

# generate slopped bed file
bedtools slop -i mutationlist.bed -g "!{genomeFile}" -b !{params.SLOP_BP} > slopped.bed

# Filter lines that are <= 100bp in range.
Rscript --vanilla "!{params.projectHome}/R/1_parse/filterSloppedBed.R" slopped.bed "!{filteredBedFile}"

# Clean up
rm -f mutationlist.bed slopped.bed

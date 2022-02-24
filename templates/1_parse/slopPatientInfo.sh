#!/bin/bash

set -euo pipefail

# Convert patient CSV to mutation list in BED format.
Rscript --vanilla "!{params.projectHome}/R/1_parse/patientListCsvToBed.R" \
    "!{csvFile}" mutationlist.bed

# Convert sequence dictionary into simple chromosome<tab>length format.
python3 "!{params.projectHome}/python/1_parse/toSlopIndex.py" \
    "!{dictionaryFile}" slopindex.tsv

# generate slopped bed file
bedtools slop -i mutationlist.bed -g slopindex.tsv -b !{params.SLOP_BASES} > slopped.bed

# Filter lines that are <= 100bp in range.
Rscript --vanilla "!{params.projectHome}/R/1_parse/filterSloppedBed.R" slopped.bed "!{filteredBedFile}"

# Clean up
rm -f mutationlist.bed slopped.bed slopindex.tsv

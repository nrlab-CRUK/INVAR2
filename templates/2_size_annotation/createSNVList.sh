#!/bin/bash

set -euo pipefail

# Convert patient CSV to mutation list in BED format.
# (Can reuse the script from part one.)

Rscript --vanilla "!{params.projectHome}/R/1_parse/patientListCsvToBed.R" \
    "!{csvFile}" mutationlist.bed

# Convert the BED file to the format needed for getFragmentSize.py

cut --output-delimiter ':' --fields '1,3,4,5' mutationlist.bed \
    | uniq > "!{snvListFile}"

# Clean up

##rm -f mutationlist.bed

#!/bin/bash

set -euo pipefail

# If the index symbolic  link doesn't resolve, remove the link.
# samtools will automatically generate the index it needs in that case.

if ! [ -e "!{fastaIndex}" ]
then
    rm -f "!{fastaIndex}"
fi

python3 "!{params.projectHome}/python/1_parse/toRegions.py" "!{mutationFile}" 1 "!{params.CHROMOSOME_ID_PREFIX}" > regions.txt

samtools faidx -r regions.txt "!{fastaReference}" > "!{trinucleotideFile}"

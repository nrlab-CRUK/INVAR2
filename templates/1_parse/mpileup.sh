#!/bin/bash

set -euo pipefail

# If the index symbolic  link doesn't resolve, remove the link.
# samtools will automatically generate the index it needs in that case.

if ! [ -e "!{fastaIndex}" ]
then
    rm -f "!{fastaIndex}"
fi

# Give a larger -d option than the default 8000 as that is too few for some cases.

samtools mpileup -x -d 100000 -q !{params.MAPPING_QUALITY} -Q !{params.BASE_QUALITY} !{dedupFlags} \
    -g -l "!{sloppedBedFile}" \
    -t "DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR" \
    -f "!{fastaReference}" \
    "!{bamFile}" \
| bcftools call -A -m -p 0 -o "!{vcfFile}"

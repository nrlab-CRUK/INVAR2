#!/bin/bash

set -euo pipefail

# If the index symbolic  link doesn't resolve, remove the link.
# samtools will automatically generate the index it needs in that case.

if ! [ -e "!{fastaIndex}" ]
then
    rm -f "!{fastaIndex}"
fi

samtools mpileup -x -d !{params.MPILEUP_MAXIMUM_DEPTH} -q !{params.MAPPING_QUALITY} -Q !{params.BASE_QUALITY} !{dedupFlags} \
    -g -l "!{sloppedBedFile}" \
    -t "DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR" \
    -f "!{fastaReference}" \
    "!{bamFile}" \
| bcftools call -A -m -p 0 -o "!{vcfFile}"

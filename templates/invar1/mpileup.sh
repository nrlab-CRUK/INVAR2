#!/bin/bash

set -euo pipefail

samtools mpileup -x -d !{params.MAX_DP} -q !{params.MAPQ} -Q !{params.BASEQ} !{dedupFlags} \
    -g -l !{sloppedBedFile} \
    -t "DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR" \
    -f !{fastaReference} \
    !{bamFile} \
| bcftools call -A -m -p 0 -o "!{vcfFile}"

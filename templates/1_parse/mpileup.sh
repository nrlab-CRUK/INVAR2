#!/bin/bash

set -euo pipefail

# If the index symbolic  link doesn't resolve, remove the link.
# samtools will automatically generate the index it needs in that case.

if ! [ -e "!{fastaIndex}" ]
then
    rm -f "!{fastaIndex}"
fi

# Give a larger -d option than the default 8000 as that is too few for some cases.
# Seems that the recommended move to "bcftools mpileup" isn't quite as simple
# as just changing the command.
# -t doesn't exist.
# An unexpected contig results in the program stopping. 

samtools mpileup --ignore-overlaps --max-depth 100000 \
    --min-MQ !{params.MAPPING_QUALITY} --min-BQ !{params.BASE_QUALITY} !{dedupFlags} \
    -t "DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR" \
    --positions "!{sloppedBedFile}" \
    --fasta-ref "!{fastaReference}" \
    --BCF "!{bamFile}" \
| bcftools call --keep-alts --multiallelic-caller --pval-threshold 0 \
    --output "!{vcfFile}"

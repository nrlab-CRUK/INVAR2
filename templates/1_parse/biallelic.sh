#!/bin/bash

set -euo pipefail

## split multiallelic sites into biallelic records

bcftools norm -m -both "!{vcfFile}" -o split1.vcf

## get cols of interest

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t[%ADF]\t[%ADR]\t%INFO/MQSB\n' \
    split1.vcf -o split2.tsv

## Add columns, headers and filter

Rscript --vanilla "!{params.projectHome}/R/1_parse/formatAndAnnotateMutations.R" \
    split2.tsv \
    "!{mutationsFile}" \
    !{params.MPILEUP_MINIMUM_DEPTH} \
    "!{sampleId}"

## Clean up

## rm -f split1.vcf split2.tsv

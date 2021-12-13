#!/bin/bash -l

python3 "!{projectDir}/python/invar1/toRegions.py" "!{mutationFile}" 1 chr > regions.txt

samtools faidx -r regions.txt "!{fastaReference}" > "!{trinucleotideFile}"

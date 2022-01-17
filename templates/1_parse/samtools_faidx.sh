#!/bin/bash

set -euo pipefail

python3 "!{params.projectHome}/python/1_parse/toRegions.py" "!{mutationFile}" 1 chr > regions.txt

samtools faidx -r regions.txt "!{fastaReference}" > "!{trinucleotideFile}"

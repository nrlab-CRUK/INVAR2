#!/bin/bash -l

set -euo pipefail

python3 "!{projectDir}/python/invar1/toRegions.py" "!{mutationFile}" > regions.txt

tabix -R regions.txt "!{tabixDatabase}" > "!{tabixFile}"

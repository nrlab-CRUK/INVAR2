#!/bin/bash

set -euo pipefail

cut -f1,2 "!{mutationFile}" | sed '1d;s/^chr//i' > regions.txt

tabix -R regions.txt --separate-regions "!{tabixDatabase}" > "!{tabixFile}"

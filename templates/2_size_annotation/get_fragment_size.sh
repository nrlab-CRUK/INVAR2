#!/bin/bash

set -euo pipefail

if ! [ -e "!{bamIndex}" ]
then
    samtools index "!{bamFile}"
fi

python3 "!{projectDir}/python/2_size_annotation/getFragmentSize.py" \
        "!{bamFile}" \
        "!{snvList}" \
        "!{insertsFile}"

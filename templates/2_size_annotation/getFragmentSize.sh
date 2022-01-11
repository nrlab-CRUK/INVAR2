#!/bin/bash

set -euo pipefail

# The exists test (-e) will be false if a symbolic link is present but
# its target is not.

if ! [ -e "!{bamIndex}" ]
then
    # Remove the symbolic link and create the index as a file.

    rm -f "!{bamIndex}"
    samtools index "!{bamFile}" "!{bamIndex}"
fi

python3 "!{projectDir}/python/2_size_annotation/getFragmentSize.py" \
        "!{bamFile}" \
        "!{snvList}" \
        "!{insertsFile}"

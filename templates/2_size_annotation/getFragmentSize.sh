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

python3 "!{params.projectHome}/python/2_size_annotation/getFragmentSize.py" \
        --min-length=!{params.MINIMUM_FRAGMENT_LENGTH} \
        --max-length=!{params.MAXIMUM_FRAGMENT_LENGTH} \
        --min-mapping-quality=!{params.MAPPING_QUALITY} \
        --min-base-quality=!{params.BASE_QUALITY} \
        "!{bamFile}" \
        "!{snvList}" \
        "!{insertsFile}"

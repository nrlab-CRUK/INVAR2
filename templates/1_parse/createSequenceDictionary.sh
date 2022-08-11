#!/bin/bash

# Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary

picard CreateSequenceDictionary \
    -Xms!{javaMem}m -Xmx!{javaMem}m \
    --REFERENCE "!{fastaReference}" \
    --OUTPUT "!{sequenceDictionary}"

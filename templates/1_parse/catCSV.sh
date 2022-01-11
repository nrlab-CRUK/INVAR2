#!/bin/bash

set -euo pipefail

# Cat files together taking the header from the first file and skipping the
# header in all subsequent files.

head --quiet -n 1 "!{csvFiles[0]}" > "!{combinedFile}"

tail --quiet -n +2 "!{csvFiles.join('" "')}" >> "!{combinedFile}"

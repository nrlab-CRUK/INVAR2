#!/bin/sh

python3 "!{projectDir}/python/invar1/toRegions.py" "!{mutationFile}" | split -l 500 - "regions.txt."

rm -f "!{tabixFile}"

for regions in regions.txt.*
do
    tabix "!{tabixDatabase}" $(cat $regions) >> "!{tabixFile}"
done

rm regions.txt.*

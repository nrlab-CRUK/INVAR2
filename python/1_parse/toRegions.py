#!/usr/bin/env python3

# Build a file of regions for Tabix and samtools.

import csv
import sys

# sys.argv[0] is the script path.
mutationsFile = sys.argv[1]
offset = int(sys.argv[2]) if len(sys.argv) > 2 else 0

with open(mutationsFile, 'r') as fp1:
    csvReader = csv.DictReader(fp1, delimiter = '\t')
    for mutation in csvReader:
        chr = mutation['CHROM']
        pos = int(mutation['POS'])
        print(f"{chr}:{pos - offset}-{pos + offset}")

#!/usr/bin/env python3

# Convert a reference sequencing dictionary to the simple tab separated
# chromosome, length file required by bedtools slop.

import csv
import re
import sys

# sys.argv[0] is the script path.
sequenceDictionary = sys.argv[1]
indexFile = sys.argv[2]

colonSplit = re.compile(r':')

with open(sequenceDictionary, 'r') as fp1:
    csvReader = csv.reader(fp1, delimiter = '\t')
    with open(indexFile, 'w') as fp2:
        for line in csvReader:
            if line[0] == '@SQ':
                map = {}
                for pair in line[1::]:
                    kv = colonSplit.split(pair)
                    map[kv[0]] = kv[1]
                fp2.write("{}\t{}\n".format(map['SN'], map['LN']))

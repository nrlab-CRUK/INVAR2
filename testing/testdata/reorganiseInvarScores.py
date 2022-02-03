#!/usr/bin/env python3

import re
import sys

originalOutput = sys.argv[1]
convertedOutput = sys.argv[2]

with open(originalOutput, "r") as input:
    with open(convertedOutput, "w") as output:
        first = True
        for line in input:
            if 'INVAR_SCORE' in line:
                search = re.search('^(.*):?\\[1\\] "(.*)"$', line.strip())
                if search is not None:

                    parts = search.group(2).split(':')
                    headers = parts[0].strip().split(',')
                    values = parts[1].strip().split(' ')

                    usingSize = 'using_size' in headers[0]

                    headers[0] = 'INVAR_SCORE'
                    headers.append("USING_SIZE")
                    values.append("TRUE" if usingSize else "FALSE")

                    poolSearch = re.search("(SLX-\\d+)_(\\w+)\\.", search.group(1))

                    if poolSearch is not None:
                        headers.append("POOL")
                        headers.append("BARCODE")
                        values.append(poolSearch.group(1))
                        values.append(poolSearch.group(2))

                    if first:
                        output.write(",".join(headers))
                        output.write('\n')
                        first = False

                    output.write(",".join(values))
                    output.write('\n')



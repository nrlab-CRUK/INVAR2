#!/usr/bin/env python3

import re
import subprocess
import sys

from os.path import exists

for source in sys.argv:
    matcher = re.search("SLX-(\\d+)_(\\w+)\\.(\\d+)", source)
    if matcher is not None:
        pool = matcher.group(1)
        barcode = matcher.group(2)
        sample = matcher.group(3)

        outputFile = f"mutation_table.outliersuppressed.SLX{pool}{barcode}.{sample}.rds"

        if not exists(outputFile):
            subprocess.run(['Rscript', '../../reorganiseMutationTable.R', source, outputFile])


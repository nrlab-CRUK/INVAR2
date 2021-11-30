#!/usr/bin/env python3

# Slightly rewritten version of the original ann.py to read from files
# from other processes rather than running commands within this script.

import csv
import json
import sys

mutation           = json.loads(sys.argv[1])
snp_tabix_file     = sys.argv[2]
cosmic_tabix_file  = sys.argv[3]
trinucleotide_file = sys.argv[4]
mutation_file      = sys.argv[5]

# The mutation file needs the following column headings
# CHROM    POS    REF    ALT

chromosome  = mutation['CHROM'].upper()
position    = mutation['POS']
reference   = mutation['REF']
alternative = mutation['ALT']

chromosome_short  = chromosome.replace("CHR", "")

def get_tabix_info(values):
    infoDict = {}
    for pair in values.split(';'):
        try:
            key,value = pair.split('=')
            infoDict[key] = value
        except ValueError:
            ## for boolean field use 1 as the value
            infoDict[pair] = '1'
    return infoDict


k1g_af = 0
with open(snp_tabix_file, 'r') as fp:
    for vcf_line in fp:
        vcf_line = vcf_line.strip()

        vcf_values     = vcf_line.split("\t")
        k1g_chromosome = vcf_values[0]
        k1g_position   = vcf_values[1]
        k1g_reference  = vcf_values[3]
        k1g_alt        = vcf_values[4]
        k1g_filter     = vcf_values[6]

        infoDict = get_tabix_info(vcf_values[7])

        k1g_vt = infoDict['VT']

        print(f"SNP: {chromosome_short}=={k1g_chromosome} {position}=={k1g_position} {reference}=={k1g_reference} {alternative}=={k1g_alt} {k1g_filter}==PASS {k1g_vt}==SNP")

        # sanity check values match
        if (chromosome_short == k1g_chromosome and
            position    == k1g_position and
            reference   == k1g_reference and
            alternative == k1g_alt and
            k1g_filter  == "PASS" and
            k1g_vt      == "SNP"):

            k1g_af = infoDict['AF']
            break


cosmic_sampleN = 0
cosmic_snp = 0
with open(cosmic_tabix_file, 'r') as fp:
    for vcf_line in fp:
        vcf_line = vcf_line.strip()

        cosmic_values    = vcf_line.split("\t")
        cosmic_chromsome = cosmic_values[0]
        cosmic_position  = cosmic_values[1]
        cosmic_reference = cosmic_values[3]
        cosmic_alt       = cosmic_values[4]
        cosmic_filter    = cosmic_values[6]
        
        infoDict = get_tabix_info(cosmic_values[7])

        cosmic_sampleN = infoDict['CNT']

        try:
            cosmic_snp = infoDict['SNP']
        except KeyError:
            cosmic_snp = 0

        print(f"COSMIC: {chromosome_short}=={cosmic_chromsome} {position}=={cosmic_position} {reference}=={cosmic_reference} {alternative}=={cosmic_alt} {cosmic_filter}==PASS")

        # sanity check values match
        if (chromosome_short == cosmic_chromsome and
            position         == cosmic_position and
            reference        == cosmic_reference and
            alternative      == cosmic_alt and
            cosmic_filter    == "PASS"):

            cosmic_sampleN = infoDict['CNT']
            break


trinucleotide = ''
with open(trinucleotide_file, 'r') as fp:
    fp.readline()   # Skip the header line
    trinucleotide = fp.readline().strip()


headers = list(mutation.keys())
headers.extend(['COSMIC_MUTATIONS', 'COSMIC_SNP', '1KG_AF', 'TRINUCLEOTIDE'])
values = list(mutation.values())
values.extend([str(cosmic_sampleN), str(cosmic_snp), str(k1g_af), trinucleotide])

with open(mutation_file, 'w') as fp:
    print('\t'.join(headers), file = fp)
    print('\t'.join(values), file = fp)

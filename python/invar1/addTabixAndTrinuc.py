#!/usr/bin/env python3

import csv
import re
import sys

# See https://stackoverflow.com/a/15063941
csv.field_size_limit(sys.maxsize)

mutation_file      = sys.argv[1]
snp_tabix_file     = sys.argv[2]
cosmic_tabix_file  = sys.argv[3]
trinucleotide_file = sys.argv[4]
annotated_file     = None

if len(sys.argv) > 5 and sys.argv[5] != '-':
    annotated_file = sys.argv[5]

chr_remover = re.compile('^chr', re.IGNORECASE)
fasta_extractor = re.compile(r'^>(chr)?([A-Z0-9]+):([0-9]+)-([0-9]+)$', re.IGNORECASE)

def get_position(mutation):
    chr = chr_remover.sub('', mutation['CHROM']).upper()
    pos = int(mutation['POS'])
    return chr, pos

def get_position_key(mutation):
    chr, pos = get_position(mutation)
    return frozenset([chr, pos])


def get_tabix_info(values):
    infoDict = {}
    for pair in values.split(';'):
        keyValue = pair.split('=')
        # for boolean field (no value) use 1 as the value
        infoDict[keyValue[0]] = keyValue[1] if len(keyValue) > 1 else '1'
    return infoDict

def read_tabix_file(file):
    tabixDict = dict()
    with open(file, 'r') as fp:
        reader = csv.reader(fp, delimiter = '\t')
        for row in reader:
            chr = row[0]
            pos = int(row[1])
            key = frozenset([chr, pos])
            
            info = { 'CHROM':chr, 'POS':pos, 'REF':row[3], 'ALT':row[4], 'FILTER': row[6] }
            info['INFO'] = get_tabix_info(row[7])
            
            tabixRows = tabixDict.get(key)
            if not tabixRows:
                tabixRows = list()
                tabixDict[key] = tabixRows
            tabixRows.append(info)
    return tabixDict

def read_fasta_file(file):
    fastaDict = dict()
    with open(file, 'r') as fp:
        while True:
            posStr = fp.readline()
            trinuc = fp.readline()
            if posStr and trinuc:
                m = fasta_extractor.match(posStr)
                if m:
                    chr = m.group(2)
                    pos = int(m.group(3)) + 1
                    key = frozenset([chr, pos])
                    fastaDict[key] = trinuc.rstrip()
            else:
                break
    return fastaDict


def addSnpInfo(mutation, tabix_info):
    chr, pos = get_position(mutation)
    key = frozenset([chr, pos])
    
    mutation['1KG_AF'] = '0'
    
    snpInfo = tabix_info.get(key)
    if snpInfo is not None:
        for info in snpInfo:
            #print(f"SNP: {chr}=={info['CHROM']} {pos}=={info['POS']} {mutation['REF']}=={info['REF']} {mutation['ALT']}=={info['ALT']} {info['FILTER']}==PASS {info['INFO'].get('VT')}==SNP", file = sys.stderr)

            if (info['CHROM'] == chr and
               info['POS'] == pos and
               info['REF'] == mutation['REF'] and
               info['ALT'] == mutation['ALT'] and
               info['FILTER'] == "PASS" and
               info['INFO']['VT'] == "SNP"):
                mutation['1KG_AF'] = info['INFO'].get('AF', '')
                
                #print(f"Set SNP info for {key} to AF = {mutation['1KG_AF']}", file = sys.stderr)
                
                break

    return mutation

def addCosmicInfo(mutation, tabix_info):
    chr, pos = get_position(mutation)
    key = frozenset([chr, pos])
    
    mutation['COSMIC_MUTATIONS'] = '0'
    mutation['COSMIC_SNP'] = '0'
    
    cosmicInfo = tabix_info.get(key)
    if cosmicInfo is not None:
        for info in cosmicInfo:
            #print(f"COSMIC: {chr}=={info['CHROM']} {pos}=={info['POS']} {mutation['REF']}=={info['REF']} {mutation['ALT']}=={info['ALT']} {info['FILTER']}==PASS", file = sys.stderr)
            
            if (info['CHROM'] == chr and
               info['POS'] == pos and
               info['REF'] == mutation['REF'] and
               info['ALT'] == mutation['ALT']):
                mutation['COSMIC_MUTATIONS'] = info['INFO'].get('CNT', '0')
                mutation['COSMIC_SNP'] = info['INFO'].get('SNP', '0')
                
                #print(f"Set COSMIC info for {key} to count = {mutation['COSMIC_MUTATIONS']} and snp = {mutation['COSMIC_SNP']}", file = sys.stderr)
                
                break

    return mutation


def addTrinucleotideInfo(mutation, trinuc_info):
    key = get_position_key(mutation)
    mutation['TRINUCLEOTIDE'] = trinuc_info.get(key, '')
    return mutation


# Main code.

snpData = read_tabix_file(snp_tabix_file)
print(f"Have {len(snpData)} positions for SNP data.", file = sys.stderr)

cosmicData = read_tabix_file(cosmic_tabix_file)
print(f"Have {len(cosmicData)} positions for COSMIC data.", file = sys.stderr)

trinucData = read_fasta_file(trinucleotide_file)
print(f"Have {len(trinucData)} positions for trinucleotides.", file = sys.stderr)

first = True

with open(mutation_file, 'r') as fp1:
    csv = csv.DictReader(fp1, delimiter = '\t')
    
    fp2 = open(annotated_file, 'w') if annotated_file else sys.stdout
    try:
        for mutation in csv:
            addCosmicInfo(mutation, cosmicData)
            addSnpInfo(mutation, snpData)
            addTrinucleotideInfo(mutation, trinucData)
            
            if first:
                print('\t'.join(mutation.keys()), file = fp2)
                first = False
    
            print('\t'.join(mutation.values()), file = fp2)
    finally:
        if annotated_file:
            fp2.close()

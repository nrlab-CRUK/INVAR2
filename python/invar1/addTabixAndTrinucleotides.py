#!/usr/bin/env python3

import csv
import re
import sys

# See https://stackoverflow.com/a/15063941
csv.field_size_limit(sys.maxsize)

mutationFile      = sys.argv[1]
snpTabixFile      = sys.argv[2]
cosmicTabixFile   = sys.argv[3]
trinucleotideFile = sys.argv[4]
annotatedFile     = None

if len(sys.argv) > 5 and sys.argv[5] != '-':
    annotatedFile = sys.argv[5]

chrRemover = re.compile('^chr', re.IGNORECASE)
fastaExtractor = re.compile(r'^>(chr)?([A-Z0-9]+):([0-9]+)-([0-9]+)$', re.IGNORECASE)


def getPosition(mutation):
    chr = chrRemover.sub('', mutation['CHROM']).upper()
    pos = int(mutation['POS'])
    return chr, pos

def getPositionKey(mutation):
    chr, pos = getPosition(mutation)
    return frozenset([chr, pos])


def getTabixInfo(values):
    infoDict = {}
    for pair in values.split(';'):
        keyValue = pair.split('=')
        # for boolean field (no value) use 'T' for the value, representing true
        infoDict[keyValue[0]] = keyValue[1] if len(keyValue) > 1 else 'T'
    return infoDict

def readTabixFile(file):
    tabixDict = dict()
    with open(file, 'r') as fp:
        reader = csv.reader(fp, delimiter = '\t')
        key = None
        for row in reader:
            if len(row) == 1 and row[0].startswith("#"):
                (chr, region) = row[0][1:].split(":")
                pos = int(region.split("-")[0])
                key = frozenset([chr, pos])
            else:
                chr = row[0]
                pos = int(row[1])
                info = { 'CHROM':chr, 'POS':pos, 'REF':row[3], 'ALT':row[4], 'FILTER': row[6] }
                info['INFO'] = getTabixInfo(row[7])
                tabixRows = tabixDict.get(key)
                if not tabixRows:
                    tabixRows = list()
                    tabixDict[key] = tabixRows
                tabixRows.append(info)
    return tabixDict

def readFastaFile(file):
    fastaDict = dict()
    with open(file, 'r') as fp:
        while True:
            posStr = fp.readline()
            trinuc = fp.readline()
            if posStr and trinuc:
                m = fastaExtractor.match(posStr)
                if m:
                    chr = m.group(2)
                    pos = int(m.group(3)) + 1
                    key = frozenset([chr, pos])
                    fastaDict[key] = trinuc.rstrip().upper()
            else:
                break
    return fastaDict


def addSnpInfo(mutation, tabixInfo):
    chr, pos = getPosition(mutation)
    key = frozenset([chr, pos])
    
    zero = '0'
    mutation['1KG_AF'] = zero
    
    snpInfo = tabixInfo.get(key)
    if snpInfo is not None:
        for info in snpInfo:
            #print(f"SNP: {chr}=={info['CHROM']} {pos}=={info['POS']} {mutation['REF']}=={info['REF']} {mutation['ALT']}=={info['ALT']} {info['FILTER']}==PASS {info['INFO'].get('VT')}==SNP", file = sys.stderr)

            if (info['CHROM'] == chr and
               info['POS'] == pos and
               info['REF'] == mutation['REF'] and
               info['ALT'] == mutation['ALT'] and
               info['FILTER'] == "PASS" and
               info['INFO']['VT'] == "SNP"):
                mutation['1KG_AF'] = info['INFO'].get('AF', zero)
                
                #print(f"Set SNP info for {key} to AF = {mutation['1KG_AF']}", file = sys.stderr)
                
                break

    return mutation

def addCosmicInfo(mutation, tabixInfo):
    chr, pos = getPosition(mutation)
    key = frozenset([chr, pos])
    
    mutation['COSMIC_MUTATIONS'] = '0' # This is an integer, so default to zero.
    mutation['COSMIC_SNP'] = 'F'       # This is a boolean, so use T or F.

    cosmicInfo = tabixInfo.get(key)
    if cosmicInfo is not None:
        for info in cosmicInfo:
            #print(f"COSMIC: {chr}=={info['CHROM']} {pos}=={info['POS']} {mutation['REF']}=={info['REF']} {mutation['ALT']}=={info['ALT']} {info['FILTER']}==PASS", file = sys.stderr)

            #print(f"Set COSMIC info for {key} to count = {mutation['COSMIC_MUTATIONS']} and snp = {mutation['COSMIC_SNP']}", file = sys.stderr)

            # TODO
            # Should we check that the position and alternate allele match?
            # (the original ann.py does not but it looks like a bug)
            # This would need to be handled carefully since COSMIC has mutations
            # such as the following where matching on the position field would
            # not work.
            # tabix ../resources/COSMIC/v82/CosmicCodingMuts.vcf.gz 1:12785297-12785297
            #   1      12785295      COSM5377462    CGG    CAA    .      .       GENE=AADACL3;STRAND=+;CDS=c.176_177GG>AA;AA=p.R59Q;CNT=1
            #   1      12785295      COSM5377461    CGG    CAA    .      .       GENE=AADACL3_ENST00000359318;STRAND=+;CDS=c.386_387GG>AA;AA=p.R129Q;CNT=1

            mutation['COSMIC_MUTATIONS'] = info['INFO'].get('CNT', '0')
            mutation['COSMIC_SNP'] = info['INFO'].get('SNP', 'F')

    return mutation

def addTrinucleotideInfo(mutation, trinucleotideInfo):
    key = getPositionKey(mutation)
    mutation['TRINUCLEOTIDE'] = trinucleotideInfo.get(key, '')
    return mutation


# Main code.

snpData = readTabixFile(snpTabixFile)
print(f"Have {len(snpData)} positions for SNP data.", file = sys.stderr)

cosmicData = readTabixFile(cosmicTabixFile)
print(f"Have {len(cosmicData)} positions for COSMIC data.", file = sys.stderr)

trinucleotideData = readFastaFile(trinucleotideFile)
print(f"Have {len(trinucleotideData)} positions for trinucleotides.", file = sys.stderr)

first = True

with open(mutationFile, 'r') as fp1:
    csv = csv.DictReader(fp1, delimiter = '\t')
    
    fp2 = open(annotatedFile, 'w') if annotatedFile else sys.stdout
    try:
        for mutation in csv:
            addCosmicInfo(mutation, cosmicData)
            addSnpInfo(mutation, snpData)
            addTrinucleotideInfo(mutation, trinucleotideData)
            
            if first:
                print('\t'.join(mutation.keys()), file = fp2)
                first = False
    
            print('\t'.join(mutation.values()), file = fp2)
    finally:
        if annotatedFile:
            fp2.close()

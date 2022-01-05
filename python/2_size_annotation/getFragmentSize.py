#!/usr/bin python3

import argparse
import pysam
import os
import sys

def processPileUps(samfile, output_file, chrom, pos, args):
    for pileup in samfile.pileup(chrom, pos-1, pos, max_depth=args.max_depth, stepper="nofilter"):
        if pileup.reference_pos == pos-1: # filter for position of interest

            if False:
                print("Processing {} reads covering SNV position {}:{} in {} with SNV as {}".format(
                      len(pileup.pileups), chrom, pos, args.input_bam, SNV))

            for read in pileup.pileups:
                if read.query_position:
                    SNV_base = read.alignment.query_sequence[read.query_position]
                    mutant_status = "."

                    if SNV_base == ref:
                        mutant_status = '0'
                    elif SNV_base == alt:
                        mutant_status = '1'

                    output_file.write("{},{},{},{:d}\n".format(SNV_base, SNV.strip(), mutant_status, abs(read.alignment.template_length)))


# required variables
parser = argparse.ArgumentParser(description="getFragmentSize.py - "
    "Extracts reads from a BAM file that cover a specified SNV and "
    "writes the reference and alternate allele containing reads to separate BAM files.")

parser.add_argument("input_bam", help="Input BAM file. Has to be indexed.")

parser.add_argument("SNV_list", help="List of SNV positions (.txt), reference and alternate "
    "allele of interest in chr:position:ref:alt format, e.g. chr21:11106932:A:G. "
    "Coordinates are 1-based.")

parser.add_argument("output", help="Name of your output file")

# optional variables
parser.add_argument("--max-depth", dest="max_depth", type=int, default=1000000,
                    help="Maximum number of reads to process at the specified SNV position")

args = parser.parse_args()

if args.max_depth < 10000:
    print("Specified max_depth is too low - changing to 10000")
    args.max_depth = 10000

nucleotides = ['A', 'C', 'G', 'T']

with pysam.AlignmentFile(args.input_bam, "rb") as samfile:
    with open(args.output, "w") as output_file:

        # start loop to pileup at each of the SNV sites based on SNV_list
        # this is to get the ref/alt status, and the fragment length (TLEN)
        with open(args.SNV_list, "r") as SNVs:
            for SNV in SNVs:

                # do check on the format of the SNV_list
                try:
                    chrom, pos, ref, alt = SNV.strip().split(":")
                except ValueError:
                    print(f"SNV specified '{SNV}' not in chr:position:ref:alt format")
                    sys.exit(1)

                # set case
                SNV_ref = ref.upper()
                SNV_alt = alt.upper()

                # check nucleotide validity
                if ref not in nucleotides:
                    print(f"Reference allele {ref} is not A, C, G or T")
                    continue

                if alt not in nucleotides:
                    print(f"Alternate allele {alt} is not A, C, G or T")
                    continue

                # Check validity of pos - is it an integer?
                try:
                    pos = int(pos)
                except ValueError:
                    print(f"Position {pos} is not valid")
                    continue

                # get the read
                try:
                    processPileUps(samfile, output_file, chrom, pos, args)
                except ValueError as error:
                    # Can happen when there are contigs not in the index.
                    if not 'invalid contig' in error.args[0]:
                        raise error

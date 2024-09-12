#!/usr/bin/env python3

import argparse
import pysam
import csv
import json
import sys
import re




# function adapted from https://github.com/BCCDC-PHL/vcf-diff/blob/main/vcf-diff#L123C1-L169C15
def parse_vcf(vcf_path):

    vcf_dict = []
    with open(vcf_path, mode='r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            else:
                line_split = line.strip().split('\t')
                variant = {}
                variant['chrom'] = line_split[0]
                variant['pos'] = int(line_split[1])
                variant['id'] = line_split[2]
                variant['ref'] = line_split[3]
                variant['alt'] = line_split[4]
                variant['qual'] = float(line_split[5])
                variant['filter'] = line_split[6]

                vcf_dict.append(variant)

    return vcf_dict


def parse_snps(snp_path):

    snp_dict = []
    with open(snp_path, mode='r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            variant = {}
            variant['pos'] = int(line_split[0])
            variant['alt'] = line_split[1]
            snp_dict.append(variant)

    return snp_dict

def detect_snps(vcf, target_snps):
    
    detected_snps = []
    for variant in vcf:
        for target in target_snps:
            # Check if 'position' and 'ref' values match
            if variant['pos'] == target['pos'] and variant['alt'] == target['alt']:
                detected_snps.append(variant)



def main(args):

    vcf = parse_vcf(args.vcf)
    nml_orygis_snps = parse_snps(args.snps)
    result = detect_snps(vcf, nml_orygis_snps)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check tbprofiler whole genome vcf for NML orgyis SNPs')
    parser.add_argument('-v', '--vcf', type=str, help='tbprofiler whole genome vcf')
    parser.add_argument('-s', '--snps', type=str, help='NML orygis snps')
    args = parser.parse_args()
    main(args)
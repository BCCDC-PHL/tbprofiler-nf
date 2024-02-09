#!/usr/bin/env python3


import argparse
import json

def parse_fasta(fasta):
    """
    Parse a fasta file into a header and a sequence.
    """
    with open(fasta, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip()
    sequence = ''.join(lines[1:])
    sequence = sequence.replace('\n', '')
    parsed_fasta = {
        'header': header,
        'sequence': sequence
    }

    return parsed_fasta

def parse_vcf(vcf):
    """
    """
    parsed_vcf = []
    header = []
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.strip().lstrip('#').split('\t')
            else:
                vcf_line_split = line.strip().split('\t')
                vcf_dict = dict(zip(header, vcf_line_split))
                vcf_dict['POS'] = int(vcf_dict['POS'])
                vcf_dict['INFO'] = dict([x.split('=') for x in vcf_dict['INFO'].split(';')])
                if 'ANN' in vcf_dict['INFO']:
                    vcf_dict['INFO']['ANN'] = vcf_dict['INFO']['ANN'].split(',')
                parsed_vcf.append(vcf_dict)
            
    return parsed_vcf


def apply_variants(genome, variants):
    """
    Apply variants to a reference genome.
    """
    for variant in variants:
        if variant['ALT'] == '.':
            continue
        genome['sequence'] = genome['sequence'][:variant['POS']-1] + variant['ALT'] + genome['sequence'][variant['POS']:]

    return genome
    


def main(args):
    genome = parse_fasta(args.genome)
    variants = parse_vcf(args.variants)
    genome = apply_variants(genome, variants)
    with open(args.output, 'w') as f:
        f.write(genome['header'] + '\n')
        f.write(genome['sequence'] + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-g', '--genome', type=str, help='Input reference genome')
    parser.add_argument('-v', '--variants', type=str, help='Variants to apply to the reference genome (vcf format)')
    parser.add_argument('-o', '--output', type=str, help='Output file path, (fasta format)')
    args = parser.parse_args()
    main(args)

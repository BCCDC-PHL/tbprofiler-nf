#!/usr/bin/env python3

import argparse
import csv
import json

from pathlib import Path

def get_genes_from_tbdb_bed(tbdb_bed: Path) -> list:
    """
    Get genes from TBDB BED

    :param tbdb_bed: TBDB BED
    :return: List of genes
    :rtype: list
    """
    genes = set()
    with open(tbdb_bed, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            gene_name = line_split[4]
            if gene_name == 'Rv0678':
                gene_name = 'mmpR5'
            genes.add(gene_name)
    genes = list(genes)

    return genes


def get_genes_from_tbdb_gff(tbdb_gff: Path, resistance_genes: list[str]) -> list:
    """
    Get genes from TBDB GFF

    :param tbdb_gff: TBDB GFF
    :return: List of genes
    :rtype: list
    """
    genes = {}
    with open(tbdb_gff, 'r') as f:
        for line in f:
            if line.startswith('##'):
                pass
            else:
                line_split = line.strip().split('\t')
                if len(line_split) < 9:
                    continue
                start_pos = int(line_split[3]) - 1
                end_pos = int(line_split[4])
                strand = line_split[6]
                description = line_split[8]
                description_split = description.split(';')
                for item in description_split:
                    if item.startswith('Name='):
                        gene_name = item.split('=')[1]
                        if gene_name in resistance_genes:
                            genes[gene_name] = {
                                'start': start_pos,
                                'end': end_pos,
                                'strand': strand,
                            }
                    elif item.startswith('ID='):
                        gene_id = item.split('=')[1].split(':')[1]
                        if gene_id in resistance_genes:
                            genes[gene_id] = {
                                'start': start_pos,
                                'end': end_pos,
                                'strand': strand,
                            }

    return genes
            

def main(args):
    resistance_genes = get_genes_from_tbdb_bed(args.tbdb_bed)
    gene_positions_by_gene_name = get_genes_from_tbdb_gff(args.tbdb_gff, resistance_genes)
    for gene_name, gene_positions in gene_positions_by_gene_name.items():
        print(f"{args.chrom_name}\t{gene_positions['start']}\t{gene_positions['end']}\t{gene_name}\t0\t{gene_positions['strand']}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate BED file of genes from TBDB')
    parser.add_argument('--tbdb-gff', type=Path, help='tbdb genome.gff file', required=True)
    parser.add_argument('--tbdb-bed', type=Path, help='tbdb tbdb.bed file', required=True)
    parser.add_argument('--chrom-name', type=str, default='NC_000962.3', help='Chromosome name')
    args = parser.parse_args()
    main(args)
    

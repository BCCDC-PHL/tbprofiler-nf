#!/usr/bin/env python3

import argparse
import csv
import json

from pathlib import Path


def get_gene_info_from_tbdb(tbdb_bed: Path, tbdb_gff: Path) -> list:
    """
    Get genes and positions from TBDB BED and strand from TBDB GFF

    :param tbdb_gff: TBDB GFF
    :param tbdb_bed: TBDB BED
    :return: List of genes with strand and position info
    :rtype: list
    """
    genes = {}
    
    with open(tbdb_bed, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            res_gene_name = line_split[4]
            start_pos = int(line_split[1]) - 1 #tbdb.bed 1 indexed, subtract 1 for 0 based index
            end_pos = int(line_split[2])
            drug_resistance = line_split[5]
            if res_gene_name == 'Rv0678':
                res_gene_name = 'mmpR5'

            with open(tbdb_gff, 'r') as f:
                for line in f:
                    if line.startswith('##'):
                        pass
                    else:
                        line_split = line.strip().split('\t')
                        if len(line_split) < 9:
                            continue
                        strand = line_split[6]
                        description = line_split[8]
                        description_split = description.split(';')
                        
                        for item in description_split:
                            if item.startswith('Name='):
                                gene_name = item.split('=')[1]
                                if gene_name == res_gene_name:
                                    genes[gene_name] = {
                                        'start': start_pos,
                                        'end': end_pos,  
                                        'strand': strand,
                                        'drugs' : drug_resistance,
                                        }
                            elif item.startswith('ID='):
                                gene_id = item.split('=')[1].split(':')[1]
                                if gene_id == res_gene_name:
                                    genes[gene_id] = {
                                        'start': start_pos,
                                        'end': end_pos,
                                        'strand': strand,
                                        'drugs' : drug_resistance,
                                        }
    return genes
            

def main(args):
    gene_positions_by_gene_name = get_gene_info_from_tbdb(args.tbdb_bed,args.tbdb_gff)
    for gene_name, gene_positions in gene_positions_by_gene_name.items():
        print(f"{args.chrom_name}\t{gene_positions['start']}\t{gene_positions['end']}\t{gene_name}\t0\t{gene_positions['strand']}\t{gene_positions['drugs']}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate BED file of genes from TBDB')
    parser.add_argument('--tbdb-gff', type=Path, help='tbdb genome.gff file', required=True)
    parser.add_argument('--tbdb-bed', type=Path, help='tbdb tbdb.bed file', required=True)
    parser.add_argument('--chrom-name', type=str, default='NC_000962.3', help='Chromosome name')
    args = parser.parse_args()
    main(args)
    

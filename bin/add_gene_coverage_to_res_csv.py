#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path

def process_coverage(coverage_file: Path, threshold: float) -> dict:
    """
    Calculate gene coverage metric info from resistance gene coverage csv

    :param coverage_file: Path to the gene coverage CSV file.
    :type coverage_file: pathlib.Path
    :param threshold: The percentage threshold to consider a gene with complete coverage (default 100%).
    :type threshold: float
    :return: Dictionary wtih number of genes and genes with complete coverage above threshold
    :rtype: dict
    """
    coverage_data = {}

    with coverage_file.open('r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            drug = row['drug']
            percent_coverage = float(row['percent_of_gene_covered_above_depth_threshold'])

            if drug not in coverage_data:
                coverage_data[drug] = {
                    'total_genes': 0, 
                    'fully_covered_genes': 0
                    }

            coverage_data[drug]['total_genes'] += 1
            if percent_coverage >= threshold:
                coverage_data[drug]['fully_covered_genes'] += 1

    return coverage_data

def calculate_metrics(coverage_data: dict) -> dict:
    """
    Calculate coverage metrics for each drug based on gene coverage data.

    :param coverage_data: A dictionary containing total and fully covered genes for each drug.
    :type coverage_data: dict
    :return: A dictionary with calculated metrics including the number of genes assessed, 
             the number of genes with complete coverage, and whether all genes have full coverage.
    :rtype: dict
    """
    metrics = {}
    for drug, data in coverage_data.items():
        num_genes_assessed = data['total_genes']
        num_genes_with_complete_coverage = data['fully_covered_genes']
        all_genes_with_complete_coverage = num_genes_assessed == num_genes_with_complete_coverage

        metrics[drug] = {
            'num_genes_assessed': num_genes_assessed,
            'num_genes_with_complete_coverage': num_genes_with_complete_coverage,
            'all_genes_with_complete_coverage': all_genes_with_complete_coverage
        }

    return metrics

def main(args: argparse.Namespace) -> None:
    # Process the coverage data to get gene coverage metrics
    threshold = args.threshold
    coverage_data = process_coverage(args.coverage, threshold)
    coverage_metrics = calculate_metrics(coverage_data)

    # Add metrics to resistance file
    output_file = args.output
    with args.resistance.open('r') as res_file, output_file.open('w', newline='') as out_file:
        res_reader = csv.DictReader(res_file)
        fieldnames = res_reader.fieldnames + ['all_genes_with_complete_coverage', 'num_genes_assessed', 'num_genes_with_complete_coverage']
        writer = csv.DictWriter(out_file, fieldnames=fieldnames)
        writer.writeheader()

        for row in res_reader:
            drug = row['drug'].lower()  
            if drug in coverage_metrics:
                metrics = coverage_metrics[drug]
                row['all_genes_with_complete_coverage'] = metrics['all_genes_with_complete_coverage']
                row['num_genes_assessed'] = metrics['num_genes_assessed']
                row['num_genes_with_complete_coverage'] = metrics['num_genes_with_complete_coverage']
            else:
                row['all_genes_with_complete_coverage'] = 'NA'
                row['num_genes_assessed'] = 'NA'
                row['num_genes_with_complete_coverage'] = 'NA'
            
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate resistance gene coverage metrics.')
    parser.add_argument('--resistance', required=True, type=Path, help='Path for resistance CSV file.')
    parser.add_argument('--coverage', required=True, type=Path, help='Path for gene coverage CSV file.')
    parser.add_argument('--threshold', type=float, default=100.0, help='Threshold for complete gene coverage (default: 100).')
    parser.add_argument('--output', required=True, type=Path, help='Output CSV file name.')
    args = parser.parse_args()
    main(args)
#!/usr/bin/env python3

import argparse
import csv
import json
import sys
import re


def parse_snpit(snpit_path):
    """
    Parse the snpit output. We assume that the first line is
    the header, the second line is the results, and there
    is only one line of results (for a single sample).

    :param snpit_path: Path to the snpit output
    :type snpit_path: str
    :return: header, snpit
    :rtype: list, dict
    """
    snpit = {}
    header = []
    with open(snpit_path, 'r') as f:
        header = next(f).strip().split('\t')
        
    with open(snpit_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            snpit = row
            
    return header, snpit


def main(args):
    snpit_header, snpit_results = parse_snpit(args.snpit)

    if snpit_results['Species'] == 'N/A':
        tbprofiler_report = json.load(open(args.tbprofiler_report))
        tbprofiler_lineages = tbprofiler_report['lineage']
        for lineage in tbprofiler_lineages:
            if re.search('H37Rv', lineage['family']) and lineage['frac'] > args.min_frac:
                snpit_results['Species'] = 'M. tuberculosis'
                snpit_results['Lineage'] = tbprofiler_report['main_lin'].replace('lineage', 'Lineage ')
                snpit_results['Percentage'] = 'N/A'
                break

    output_fieldnames = snpit_header
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerow(snpit_results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a lineage file for H37Rv')
    parser.add_argument('-s', '--snpit', type=str, help='Output from snpit')
    parser.add_argument('-t', '--tbprofiler-report', type=str, help='TBProfiler full report (json)')
    parser.add_argument('-f', '--min-frac', type=float, default=0.85, help='Minimum fraction for TBProfiler lineage')
    args = parser.parse_args()
    main(args)

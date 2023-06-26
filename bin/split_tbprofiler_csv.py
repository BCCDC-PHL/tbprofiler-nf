#!/usr/bin/env python3

import argparse
import csv
import datetime
import json
import re
import os
import sys

from pathlib import Path


def parse_summary(tbprofiler_csv_path: Path):
    """
    Parse the tbprofiler summary csv report.

    :param tbprofiler_csv_path: Path to tbprofiler_summary.csv file.
    :type tbprofiler_csv_path: pathlib.Path
    :return: 2-element tuple. First element is list of headers, second element is a list containing a single dict
             keys in dict should match headers. Values: ["id", "date", "strain", "drug_resistance", "median_depth"] 
    :rtype: tuple[list[str], list[dict]]
    """
    summary_lines = []
    with open(tbprofiler_csv_path, 'r') as f:
        while True:
            line = f.readline()
            if line.strip() == "Summary":
                break

        next(f) # Skip '-----' line below section header
        while True:
            line = f.readline().strip().split(',')
            summary_lines.append(line)
            
            if line[0] == "":
                break

    summary_lines = summary_lines[:-1]

    summary_record = {}
    headers = []
    for summary_line in summary_lines:
        header = summary_line[0].lower().replace('-', '_')
        value = summary_line[1]
        if header == 'strain' and value.startswith('lineage'):
            summary_record[header] = value.replace('lineage', '')
        else:
            summary_record[header] = value
        headers.append(header)

    return (headers, [summary_record])


def parse_resistance_report(tbprofiler_csv_path: Path):
    """
    Parse the tbprofiler resistance csv report.

    :param tbprofiler_csv_path:
    :type tbprofiler_csv_path: pathlib.Path
    """
    resistance = []
    with open(tbprofiler_csv_path, 'r') as f:
        while True:
            line = f.readline()
            if line.strip() == "Resistance report":
                break

        next(f) # Skip '-----' line below section header
        while True:
            line = f.readline().strip()
            resistance.append(line)
            
            if line == "":
                break

    resistance = resistance[:-1]
    
    resistance_fields = [x.lower().replace(' ', '_') for x in resistance[0].split(',')]

    parsed_resistance_records = []
    for drug in resistance[1:]:
        resistance_record = {}
        for idx, field in enumerate(resistance_fields):
            resistance_record[field] = re.split(',', drug)[idx].replace('"', '')
        mutations = []
        for mutation in resistance_record['mutations'].split(','):
            if mutation:
                mutation = mutation.strip()
                mutation = mutation.split(' ')
                gene = mutation[0]
                mutation_detail = ' '.join(mutation[1:-1])
                estimated_fraction = float(mutation[-1].strip('(').strip(')'))
                mutation = {
                    'gene': gene,
                    'estimated_fraction': estimated_fraction,
                    'mutation': mutation_detail,
                }
                mutations.append(mutation)
                
        resistance_record['mutations'] = mutations
        
        parsed_resistance_records.append(resistance_record)

    return parsed_resistance_records


def parse_section(path: Path, header: str):
    """
    Parse a specific section of the tbprofiler report.

    :param path: Path to the tbprofiler report.
    :type path: pathlib.Path
    :param header: Header for the section of the tbprofiler report to be parsed.
                   eg. 'Lineage report', 'Resistance report'
    :type header: str
    :return: 2-element tuple. First element is the field names for the records
             of that section of the report. Second element is a list of the records
             that section of the report.
    :rtype: tuple[list[str], list[dict]]
    """
    section_lines = []
    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if line.strip() == header:
                break

        next(f) # Skip '-----' line below section header
        while True:
            line = f.readline().strip()
            section_lines.append(line)
            
            if line == "":
                break

    section_lines = section_lines[:-1]

    section_fields = [x.lower().replace(' ', '_') for x in section_lines[0].split(',')]

    parsed_section_records = []

    if header == "Resistance report":
        for line in section_lines[1:]:
            record = {}

            line_list = re.split(',', line)
            record['drug'] =line_list[0]
            record['genotypic_resistance'] = line_list[1]
            record['mutations'] = line_list[2:]
            
            parsed_section_records.append(record)    

    if header == "Lineage report":
        for line in section_lines[1:]:
            record = {}
            for idx, field in enumerate(section_fields):
                record[field] = re.split(',', line)[idx].replace('"', '')
            parsed_section_records.append(record)
    
        for record in parsed_section_records:
            if record['lineage'].startswith('lineage'):
                record['lineage'] = record['lineage'].replace('lineage', '')

    return (section_fields, parsed_section_records)


def parse_lineage_report(path: Path, sample_id: str):
    """
    Parse the 'Lineage report' section of the tbprofiler report.

    :param path: Path to tbprofiler report.
    :type path: pathlib.Path
    :param sample_id: Sample ID
    :type sample_id: str
    :return: 2-element tuple. First element is list of fieldnames for each
             record of the report. Second element is list of records.
    :rtype: tuple[list[str], list[dict]]
    """
    lineage_fields, parsed_lineage = parse_section(path, "Lineage report")
    for l in parsed_lineage:
        l['sample_id'] = sample_id
    lineage_fields = ['sample_id'] + lineage_fields

    return (lineage_fields, parsed_lineage)


def parse_resistance_report(path, sample_id):
    """
    Parse the 'Resistance report' section of the tbprofiler report.

    :param path: Path to tbprofiler report.
    :type path: pathlib.Path
    :param sample_id: Sample ID
    :type sample_id: str
    :return: 2-element tuple. First element is list of fieldnames for each
             record of the report. Second element is list of records.
    :rtype: tuple[list[str], list[dict]]

    """
    resistance_fields, parsed_resistance = parse_section(path, "Resistance report")
    for r in parsed_resistance:
        r['sample_id'] = sample_id
    resistance_fields = ['sample_id'] + resistance_fields

    return resistance_fields, parsed_resistance    


def create_two_resistance_tables(resistance_report):
    """
    Take the resistance report and create two tables with the desired headers - one to report resistance and one to report mutations.

    :param resistance_report: tbprofiler full report resistance rows
    :type resistance_report: list[dict]
    :return: Two lists that report drug resistance and predicted drug mutations.
             Keys for resistance table: ['sample_id', 'drug', 'genotypic_resistance']
             Keys for mutations table: ['sample_id', 'drug', 'gene', 'mutation', 'estimated_fraction']
    :rtype: tuple[list[dict], list[dict]]
    """
  
    resistance_table = []
    mutation_table = []

    for row in resistance_report:
   
        resistance_table_row = {}
        resistance_table_row['drug'] = row['drug']
        resistance_table_row['genotypic_resistance'] = row['genotypic_resistance']
        resistance_table_row['sample_id'] = row['sample_id']

        resistance_table.append(resistance_table_row)

        for mutation in row['mutations']:
            mutation_table_row = {}

            mutation_table_row['drug'] = row['drug']
            mutation_table_row['sample_id'] = row['sample_id']

            mutation = mutation.strip()

            mutation = mutation.split(' ')

            mutation_table_row['gene'] = mutation[0]
            mutation_detail = ' '.join(mutation[1:-1])
            mutation_table_row['mutation'] = mutation_detail

            estimated_fraction = mutation[-1].strip('(').strip(')')
            mutation_table_row['estimated_fraction'] = estimated_fraction

            mutation_table.append(mutation_table_row)

    return resistance_table, mutation_table


def main(args):

    (summary_fields, summary) = parse_summary(args.tbprofiler_csv)
    
    (lineage_fields, lineage_report) = parse_lineage_report(args.tbprofiler_csv, args.sample_id)

    (resistance_fields, resistance_report) = parse_resistance_report(args.tbprofiler_csv, args.sample_id)
    
    (resistance_table, mutation_table) = create_two_resistance_tables(resistance_report)
    

    with open(args.prefix + '_tbprofiler_summary.csv', 'w', newline=os.linesep) as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        for row in summary:
            writer.writerow(row)

    with open(args.prefix + '_tbprofiler_lineage.csv', 'w', newline=os.linesep) as f:
        writer = csv.DictWriter(f, fieldnames=lineage_fields)
        writer.writeheader()
        for row in lineage_report:
            writer.writerow(row)

                        
    resistance_fields = ['sample_id', 'drug', 'genotypic_resistance']
    mutation_fields = ['sample_id','drug' , 'gene','mutation', 'estimated_fraction']       
    

    with open(args.prefix + '_tbprofiler_resistance.csv', 'w', newline=os.linesep) as f:
        writer = csv.DictWriter(f, fieldnames=resistance_fields)
        writer.writeheader()
        for row in resistance_table:
            writer.writerow(row)

    with open(args.prefix + '_tbprofiler_resistance_mutations.csv', 'w', newline=os.linesep) as f:
        writer = csv.DictWriter(f, fieldnames=mutation_fields)
        writer.writeheader()
        for row in mutation_table:
            writer.writerow(row)
  

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tbprofiler_csv', type=Path)
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-p', '--prefix')
    args = parser.parse_args()
    main(args)

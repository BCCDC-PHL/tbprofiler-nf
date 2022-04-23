#!/usr/bin/env python3

import argparse
import csv
import datetime
import json
import re
import os
import sys

def parse_summary(tbprofiler_csv_path):
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


def parse_resistance_report(tbprofiler_csv_path):
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





def parse_section(path, header):
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
    for line in section_lines[1:]:
        record = {}
        for idx, field in enumerate(section_fields):
            record[field] = re.split(',', line)[idx].replace('"', '')
        parsed_section_records.append(record)

    if header == "Lineage report":
        for record in parsed_section_records:
            if record['lineage'].startswith('lineage'):
                record['lineage'] = record['lineage'].replace('lineage', '')

    return (section_fields, parsed_section_records)


def parse_lineage_report(path, sample_id):
    lineage_fields, parsed_lineage = parse_section(path, "Lineage report")
    for l in parsed_lineage:
        l['sample_id'] = sample_id
    lineage_fields = ['sample_id'] + lineage_fields
    return (lineage_fields, parsed_lineage)

def parse_resistance_report(path, sample_id):
    resistance_fields, parsed_resistance = parse_section(path, "Resistance report")
    for r in parsed_resistance:
        r['sample_id'] = sample_id
    resistance_fields = ['sample_id'] + resistance_fields
    return resistance_fields, parsed_resistance


def main(args):

    (summary_fields, summary) = parse_summary(args.tbprofiler_csv)
    (lineage_fields, lineage_report) = parse_lineage_report(args.tbprofiler_csv, args.sample_id)
    (resistance_fields, resistance_report) = parse_resistance_report(args.tbprofiler_csv, args.sample_id)
    

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

    with open(args.prefix + '_tbprofiler_resistance.csv', 'w', newline=os.linesep) as f:
        writer = csv.DictWriter(f, fieldnames=resistance_fields)
        writer.writeheader()
        for row in resistance_report:
            writer.writerow(row)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tbprofiler_csv')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-p', '--prefix')
    args = parser.parse_args()
    main(args)

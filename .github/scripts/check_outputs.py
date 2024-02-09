#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import urllib.request

from jsonschema import validate
import yaml


def check_provenance_format_valid(provenance_files, schema):
    """
    Check that the provenance files are valid according to the schema.
    """
    for provenance_file in provenance_files:
        with open(provenance_file) as f:
            try:
                provenance = yaml.load(f, Loader=yaml.BaseLoader)
                validate(provenance, schema)
            except Exception as e:
                return False

    return True


def check_expected_mutations(resistance_mutations_files, expected_mutations_by_sample_id):
    """
    Check that the resistance mutations files contain the expected mutations.
    """
    found_mutations_by_sample = {}
    for resistance_mutations_file in resistance_mutations_files:
        with open(resistance_mutations_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row['sample_id']
                gene = row['gene']
                mutation = row['mutation']
                if sample_id not in found_mutations_by_sample:
                    found_mutations_by_sample[sample_id] = set([])
                if mutation != '':
                    found_mutations_by_sample[sample_id].add(':'.join([gene, mutation]))

    for sample_id, expected_mutations in expected_mutations_by_sample_id.items():
        if sample_id not in found_mutations_by_sample:
            return False
        if expected_mutations != found_mutations_by_sample[sample_id]:
            return False

    return True


def main(args):
    provenance_schema_url = "https://raw.githubusercontent.com/BCCDC-PHL/pipeline-provenance-schema/main/schema/pipeline-provenance.json"
    provenance_schema_path = ".github/data/pipeline-provenance.json"
    urllib.request.urlretrieve(provenance_schema_url, provenance_schema_path)

    provenance_schema = None
    with open(provenance_schema_path) as f:
        provenance_schema = json.load(f)

    provenace_files_glob = f"{args.pipeline_outdir}/**/*_provenance.yml"
    provenance_files = glob.glob(provenace_files_glob, recursive=True)

    resistance_mutations_files_glob = f"{args.pipeline_outdir}/**/*tbprofiler_resistance_mutations.csv"
    resistance_mutations_files = glob.glob(resistance_mutations_files_glob, recursive=True)

    expected_mutations_by_sample_id = {
        'NC000962.3': set([]),
        'ERR1664619': set([
            'inhA:p.Ile194Thr',
            'embA:c.-16C>T',
            'embB:p.Met306Val',
            'embB:p.Met423Thr',
            'gyrA:p.Asp94Ala',
            'rrs:n.1401A>G',
        ]),                  
    }       

    tests = [
        {
            "test_name": "provenance_format_valid",
            "test_passed": check_provenance_format_valid(provenance_files, provenance_schema),
        },
        {
            "test_name": "expected_mutations",
            "test_passed": check_expected_mutations(resistance_mutations_files, expected_mutations_by_sample_id),
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = args.output
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction='ignore')
        writer.writeheader()
        for test in tests:
            if test["test_passed"]:
                test["test_result"] = "PASS"
            else:
                test["test_result"] = "FAIL"
            writer.writerow(test)

    for test in tests:
        if not test['test_passed']:
            exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--pipeline-outdir', type=str, help='Path to the pipeline output directory')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(args)

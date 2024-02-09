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
    

def main(args):
    provenance_schema_url = "https://raw.githubusercontent.com/BCCDC-PHL/pipeline-provenance-schema/main/schema/pipeline-provenance.json"
    provenance_schema_path = ".github/data/pipeline-provenance.json"
    urllib.request.urlretrieve(provenance_schema_url, provenance_schema_path)

    provenance_schema = None
    with open(provenance_schema_path) as f:
        provenance_schema = json.load(f)

    provenace_files_glob = f"{args.pipeline_outdir}/**/*_provenance.yml"
    provenance_files = glob.glob(provenace_files_glob, recursive=True)

    tests = [
        {
            "test_name": "provenance_format_valid",
            "test_result": check_provenance_format_valid(provenance_files, provenance_schema),
        }
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    output_path = args.output
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields)
        writer.writeheader()
        for test in tests:
            if test["test_result"]:
                test["test_result"] = "PASS"
            else:
                test["test_result"] = "FAIL"
            writer.writerow(test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--pipeline-outdir', type=str, help='Path to the pipeline output directory')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(args)

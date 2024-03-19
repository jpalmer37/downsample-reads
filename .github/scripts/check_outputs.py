#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import math
import sys
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


def check_expected_coverage(fastp_file, expected_coverage_by_sample_id, tolerance=0.1):
    """
    Check that the estimated coverage matches the expected coverage (within some tolerance).
    """
    estimated_coverage_by_sample = {}

    with open(fastp_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row['sample_id']
            target_coverage = row['target_coverage']
            if target_coverage == 'original':
                continue
            target_coverage = float(target_coverage)
            estimated_coverage = float(row['estimated_coverage'])
            estimated_coverage_by_sample[sample_id] = estimated_coverage

    for sample_id, expected_coverage in expected_coverage_by_sample_id.items():
        if sample_id not in estimated_coverage_by_sample:
            return False
        if not math.isclose(estimated_coverage_by_sample[sample_id], expected_coverage, rel_tol=tolerance):
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

    fastp_file = f"{args.pipeline_outdir}/fastp.csv"

    expected_coverage_by_sample_id = {
        'NC000913': 5.0,
        'NC000962': 10.0,
    }       

    tests = [
        {
            "test_name": "provenance_format_valid",
            "test_passed": check_provenance_format_valid(provenance_files, provenance_schema),
        },
        {
            "test_name": "expected_coverage",
            "test_passed": check_expected_coverage(fastp_file, expected_coverage_by_sample_id),
        },
    ]

    output_fields = [
        "test_name",
        "test_result"
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, extrasaction='ignore')
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
    args = parser.parse_args()
    main(args)

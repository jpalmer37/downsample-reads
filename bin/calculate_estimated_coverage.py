#!/usr/bin/env python3

import argparse
import csv
import sys

def parse_genome_size(genome_size_string):
    """
    Genome size string is expected to be in a format like: '30k', '5.0m' or '3.0g'
    """
    genome_size_string = genome_size_string.lower()
    genome_size_bases  = None
    if genome_size_string.endswith("k"):
        genome_size_bases = float(genome_size_string.strip("k")) * 1_000
    elif genome_size_string.endswith("m"):
        genome_size_bases = float(genome_size_string.strip("m")) * 1_000_000
    elif genome_size_string.endswith("g"):
        genome_size_bases = float(genome_size_string.strip("g")) * 1_000_000_000
    else:
        genome_size_bases = float(genome_size_string)

    return genome_size_bases

def main(args):

    output = []
    reader = csv.DictReader(sys.stdin)
    for row in reader:
        genome_size = parse_genome_size(row['genome_size'])
        total_bases = int(row['total_bases'])
        estimated_coverage = total_bases / genome_size
        row['estimated_coverage'] = round(estimated_coverage, 3)
        output.append(row)

    output_fieldnames = reader.fieldnames + ['estimated_coverage']
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    for row in output:
        writer.writerow(row)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate estimated coverage")
    args = parser.parse_args()
    main(args)

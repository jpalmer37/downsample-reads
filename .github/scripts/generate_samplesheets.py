#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import random
import sys




def main(args):
    output = []
    downsampling_params_by_sample_id = {
        'NC000962': {
            'GENOME_SIZE': '4.4m',
            'COVERAGE': 10,
        },
        'NC000913': {
            'GENOME_SIZE': '4.6m',
            'COVERAGE': 5,
        },
    }
    r1_files = glob.glob(f"{args.fastq_dir}/*_R1*.fastq.gz")
    for r1_file in r1_files:
        sample_id = os.path.basename(r1_file).split("_")[0]
        r2_file = r1_file.replace("_R1", "_R2")
        output_record = {
            "ID": sample_id,
            "R1": os.path.abspath(r1_file),
            "R2": os.path.abspath(r2_file),
        }
        output_record.update(downsampling_params_by_sample_id[sample_id])
        output.append(output_record)

    output_fieldnames = [
        "ID",
        "R1",
        "R2",
        "GENOME_SIZE",
        "COVERAGE",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    for row in output:
        writer.writerow(row)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('fastq_dir', help='Directory containing fastq files.')
    args = parser.parse_args()
    main(args)

#!/bin/bash

set -eo pipefail

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --samplesheet_input .github/data/samplesheet.csv \
	 --outdir .github/data/test_output \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv

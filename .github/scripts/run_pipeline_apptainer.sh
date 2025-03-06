#!/bin/bash

set -eo pipefail

nextflow run main.nf \
	 -profile apptainer \
	 --cache ${HOME}/.apptainer/cache/images \
	 --samplesheet_input .github/data/samplesheet.csv \
	 --outdir .github/data/test_output \
	 --collect_outputs \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv

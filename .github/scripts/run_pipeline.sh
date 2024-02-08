#!/bin/bash

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --outdir .github/data/test_output \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv

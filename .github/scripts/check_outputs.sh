#!/usr/bin/env bash

source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate check-outputs


.github/scripts/check_outputs.py --pipeline-outdir .github/data/test_output -o artifacts/check_outputs_results.csv

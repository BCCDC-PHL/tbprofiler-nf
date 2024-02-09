#!/usr/bin/env bash

set -eo pipefail

source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate check-outputs


.github/scripts/check_outputs.py --pipeline-outdir .github/data/test_output -o artifacts/check_outputs_results.csv

grep -v 'FAIL' .github/artifacts/check_outputs_results.csv

#!/bin/bash

set -eo pipefail

echo "Prepare artifacts .." >> artifacts/test.log

mkdir -p artifacts/fastq

mv .github/data/fastq/*.fastq.gz artifacts/fastq


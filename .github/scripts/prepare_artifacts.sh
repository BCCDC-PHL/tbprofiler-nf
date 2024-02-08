#!/bin/bash

set -eo pipefail

artifacts_dir="artifacts-nextflow-${NXF_VER}"

echo "Prepare artifacts .." >> ${artifacts_dir}/test.log

mkdir -p artifacts/fastq

mv .github/data/fastq/*.fastq.gz ${artifacts_dir}/fastq


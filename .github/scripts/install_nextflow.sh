#!/bin/bash

set -eo pipefail

artifacts_dir="artifacts-nextflow-${NXF_VER}"

echo Install Nextflow .. >> ${artifacts_dir}/test.log

wget -qO- https://get.nextflow.io | bash

sudo mv nextflow /usr/local/bin/

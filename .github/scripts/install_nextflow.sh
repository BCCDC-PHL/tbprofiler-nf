#!/bin/bash

set -eo pipefail

echo Install Nextflow .. >> artifacts/test_artifact.log

wget -qO- https://get.nextflow.io | bash

sudo mv nextflow /usr/local/bin/

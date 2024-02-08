#!/bin/bash

set -eo pipefail

echo Install Nextflow .. >> artifacts/test.log

wget -qO- https://get.nextflow.io | bash

sudo mv nextflow /usr/local/bin/

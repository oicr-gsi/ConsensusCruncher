#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -regex '.*\.tsv$' -exec wc -l {} \;
ls *.tsv | sed 's/.*\.//' | sort | uniq -c

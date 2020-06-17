#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null
find . -regex '.*\.bam$' -exec samtools flagstat {} \;
ls | sed 's/.*\.//' | sort | uniq -c

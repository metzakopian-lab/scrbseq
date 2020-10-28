#!/bin/bash

BAMFILE=$1
OUTFILE=$2

OUTPUT="$OUTFILE/$(dirname $BAMFILE).gz"
echo $BAMFILE
samtools view $BAMFILE |  awk -F'\t' '{print $1 "," $18}' | awk -F'_' '{print $3}' | gzip -c > $OUTPUT
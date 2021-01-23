#!/bin/bash

BAMFILE=$1
OUTFILE=$2

samtools view $BAMFILE |  awk -F'\t' '{print $1 "," $18}' | awk -F'_' '{print $3}' | gzip -c > $OUTFILE


#!/bin/bash

while read accno; do
	fastq-dump --split-3 --gzip $accno 2>>fastq-dump.err
done < $1
gzip -d *.fastq.gz

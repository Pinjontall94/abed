#!/bin/bash

#q2aReformat(){
	# DESCRIPTION: Reformats all fastq files in current directory to fasta

	# Deps: bbtools 'reformat' tool

	echo "Reformatting merged fastqs files to fasta format..."
	for i in ${in_file[@]}; do
		reformat.sh in=$i out=${i%.fastq}.fasta
	done
#}

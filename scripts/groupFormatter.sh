#!/bin/bash

#groupFormatter
# Description: Scans 'fastas' directory for .fasta files, and
#              formats them for use with Mothur's 'make.group()' command


# Set Study Author Name (for snakemake compatibility)
AUTHOR=$1

# Set Output file
OUTPUT=$2

joinBy(){
	# Takes first argument as new delimiter, then echoes all following
	# 	arguments with that delimiter
	# Usage: joinBy <separator> <arg1> [...] <argN>

	OIFS=$IFS
	IFS=$1

	shift
	echo "$*"

	IFS=$OIF
}


# Move fastas into working directory and add to 'raw' array
mv fastas/*.fasta .

for i in *.fasta; do
	raw_fastas+=("$i")
done


# Create hyphen-delimited lists for fastas and groups
formatter_fastas=$(joinBy - ${raw_fastas[*]})
for i in ${raw_fastas[@]}; do
	group_array+=("$AUTHOR"_"${i%%.*}")
done
formatter_groups=$(joinBy - ${group_array[*]})


# Set mothur command parameters and run mothur
group_params=("fasta=$formatter_fastas," "groups=$formatter_groups,"\
	"output=$OUTPUT")
mothur <(./scripts/mothurBatch.sh make group "${group_params[*]}")

# Move fastas back to the 'fastas' folder
mv *.fasta fastas

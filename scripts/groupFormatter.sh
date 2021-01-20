#!/bin/bash

#groupFormatter
# Description: Scans 'fastas' directory for .fasta files, and
#              formats them for use with Mothur's 'make.group()' command
# Usage:
# 	groupFormatter.sh <author> <output dir/output file> <fasta 1> [...] <fasta N>


# Set Study Author Name (for snakemake compatibility)
AUTHOR=$1

# Set output dir and output filename from 2nd argument
# 	(mothur needs them specified separately)
OUTPUT_DIR=${2%/*} # Take everything left of the first '/'
OUTPUT_FILE=${2##*/} # Take everything after the last '/'

# Shift to allow use of $@ argument for fasta inputs
shift 2

# Populate the raw_fastas array with the remaining arguments
for i in $@; do
	raw_fastas+=("$i")
done

joinBy(){
	# Description: Takes first argument as new delimiter, then echoes all
	# 	following arguments with that delimiter (by temporarily manipulating
	# 	the Input Field Separator)
	# Usage: joinBy <separator> <arg 1> [...] <arg N>

	OIFS=$IFS
	IFS=$1

	shift
	echo "$*"

	IFS=$OIF
}


# Create hyphen-delimited lists for fastas and groups
formatter_fastas=$(joinBy - ${raw_fastas[*]})
for i in ${raw_fastas[@]}; do
	group_array+=("$AUTHOR"_"${i%%.*}")
done
formatter_groups=$(joinBy - ${group_array[*]})


# Set mothur command parameters and run mothur
group_params=("fasta=$formatter_fastas," "groups=$formatter_groups,"
	"outputdir=$OUTPUT_DIR,"  "output=$OUTPUT_FILE")


# Run Mothur make.group on selected fastas
mothur <(./scripts/mothurBatch.sh make group "${group_params[*]}")

sed -i "s/$AUTHOR\_fastas\///" $OUTPUT_DIR/$OUTPUT_FILE

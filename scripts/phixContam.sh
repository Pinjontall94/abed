#!/bin/bash
# phixContam.sh
# Usage: Use in a snakemake shell block to generate list of accession numbers
#   with PhiX contamination. NOTE: This creates an empty output file to satisfy
#   Snakemake if the bowtie2 aligned outputs are zero-length
#   (i.e. no contamination was found).
#
# Ex:
#   input: "PhiX_out/{sample}.merged.PhiX"
#   output: "PhiX_out/PhiX.accnos"
#   shell: "./scripts/phixContam.sh {output} {input}"

# Eat the first argument and save as $OUTPUT (lets us use $@ later)
OUTPUT=$1
shift

for i in $@; do
	if [[ -s $i ]] && [[ ! -f $OUTPUT ]]; then
	    # Search for accnos in input, and list them in output file
	    grep -Eo "[A-Z]{3,6}[0-9]+\.[0-9]+" $i | \
	        awk -F: '{ print $OUTPUT }' >> $OUTPUT
	else
	    # Create/Update empty output file
	    touch $OUTPUT
	fi
done

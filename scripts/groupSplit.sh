#!/bin/bash

#groupSplit
# Splits the screened fasta to prepare for fasta header relabeling,
# 	necessary for vsearch later

GOOD_FASTA=$1
GOOD_GROUPS=$2
OUTPUT_DIR=$3


echo "GOOD_FASTA=$GOOD_FASTA, GOOD_GROUPS=$GOOD_GROUPS"

split_params=("fasta=$GOOD_FASTA,"\
	"group=$GOOD_GROUPS," "outputdir=$OUTPUT_DIR")
mothur <(./scripts/mothurBatch.sh split groups "${split_params[*]}")

# Rename output fastas to reasonable names
for i in $OUTPUT_DIR; do
	if [[ $i =~ ".$AUTHOR" ]]; then
		# Remove everything up to "concat."
		# (i.e. result: author_year_accno.fasta)
		echo "Renaming output files to"\
			"format: author_year_accno.fasta"
		mv -v $i ${i//.good./_}
	fi
done

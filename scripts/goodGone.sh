#!/bin/bash

# goodGone.sh
# Usage: Renames the output from mothur's group.split to satisfy snakemake
# 			and later, vsearch's relabeling requirements
#
# Ex:
# goodGone.sh <output dir for group.split>
 
for i in $@; do
	if [[ $i =~ "good" ]]; then
		# Remove everything up to "concat."
		# (i.e. result: author_year_accno.fasta)
		echo "Renaming output files to"\
			"format: author_year_accno.fasta"
		mv -v $i ${i//.good./_}
	fi
done

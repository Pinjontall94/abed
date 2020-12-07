#!/bin/bash

#groupSplit(){
	# Splits the screened fasta to prepare for fasta header relabeling,
	# 	necessary for vsearch later

	for i in *.fasta; do
		if [[ $i =~ good ]]; then
			GOOD_FASTA=$i
		fi
	done; 
#	GOOD_FASTA=$AUTHOR.concat.fasta
	for i in *.groups; do
		if [[ $i =~ good ]]; then
			GOOD_GROUPS=$i
		fi
	done
#	GOOD_GROUPS=$AUTHOR.groups
	
	echo "GOOD_FASTA=$GOOD_FASTA, GOOD_GROUPS=$GOOD_GROUPS"
	
	split_params=("fasta=$GOOD_FASTA,"\
		"group=$GOOD_GROUPS")
	mothur <(mothurBatch split groups "${split_params[*]}")

	# Rename output fastas to reasonable names
	for i in *.fasta; do
		if [[ $i =~ ".$AUTHOR" ]]; then
			# Remove everything up to "concat."
			# (i.e. result: author_year_accno.fasta)
			echo "Renaming output files to"\
				"format: author_year_accno.fasta"
			mv -v $i ${i#*good.}
		fi
	done
#}

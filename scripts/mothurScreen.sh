#!/bin/bash

#mothurScreen
# Runs mothur against screening_batch.txt to screen concat'd fasta

# INPUT: mothur(), PhiX.accnos, $AUTHOR.concat.fasta,
#		$AUTHOR.groups
# OUTPUT: *.good.* *.bad.*
#
# usage: mothurscreen <concat fasta> <groupfile> <output filename>

# Move screen fasta and groupfile into working directory
mv $1 $2 -t .

concat_fasta=${1##*/}
screen_groups=${2##*/}

# Remove PhiX contamination
if [[ -s PhiX_out/PhiX.accnos ]]; then
	remove_params=("fasta=$concat_fasta,"\
		"group=$screen_groups,"\	# Not needed by mothur?
		"accnos=PhiX/PhiX.accnos")
	mothur <(./scripts/mothurBatch remove seqs "${remove_params[*]}")
	rm -v phix_batch.txt
else
	echo "PhiX.accnos is empty, skipping..."
fi

if [[ ! $concat_fasta ]] || [[ ! $screen_groups ]]; then
	echo "$0: ERROR: Zero-length input variables detected!"
	exit 1
else
	screen_params=("fasta=$concat_fasta," "group=$screen_groups,"\
		"minlength=200," "maxlength=300," "maxambig=0,"\
		"maxhomop=8")
	summary_params=("fasta=current," "processors=$THREADS")
	./scripts/mothurBatch screen seqs "${screen_params[*]}" > screening_batch.txt
	./scripts/mothurBatch summary seqs "${summary_params[*]}" >> screening_batch.txt
	./scripts/mothurBatch count groups "group=current" >> screening_batch.txt 

	echo "cat'ing screening_batch.txt"
	cat screening_batch.txt
	mothur screening_batch.txt 
	rm screening_batch.txt
fi

# Rename output file with 3rd argument
#mv 

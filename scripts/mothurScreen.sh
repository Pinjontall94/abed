#!/bin/bash

#mothurScreen
# Runs mothur against screening_batch.txt to screen concat'd fasta

# INPUT: mothur(), PhiX.accnos, $AUTHOR.concat.fasta,
#		$AUTHOR.groups
# OUTPUT: *.good.* *.bad.*
#
# Usage:
# ./scripts/mothurScreen.sh <concat fasta> <groupfile> \
# 		<PhiX accnos path> <output filename>

# Move screen fasta, groupfile, and PhiX.accnos to working directory
mv $1 $2 $3 -t .

concat_fasta=${1##*/}
screen_groups=${2##*/}

# Remove PhiX contamination
if [[ -s $3 ]]; then
	remove_params=("fasta=$concat_fasta,"\
		"group=$screen_groups,"\	# Not needed by mothur?
		"accnos=PhiX.accnos")
	mothur <(./scripts/mothurBatch.sh remove seqs "${remove_params[*]}")
	#rm -v phix_batch.txt
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
	./scripts/mothurBatch.sh screen seqs "${screen_params[*]}" > screening_batch.txt
	./scripts/mothurBatch.sh summary seqs "${summary_params[*]}" >> screening_batch.txt
	./scripts/mothurBatch.sh count groups "group=current" >> screening_batch.txt 

	echo "cat'ing screening_batch.txt"
	cat screening_batch.txt
	mothur screening_batch.txt 
	rm screening_batch.txt
fi

# Move input fastas and group file back to original directories
mv $1 -t screened
mv $2 -t mothur_in
mv $3 -t PhiX_out

# Rename output file with 4th argument
mv *.good.fasta $4

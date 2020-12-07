#!/bin/bash

#mothurScreen(){
	# Runs mothur against screening_batch.txt to screen concat'd fasta

	# INPUT: mothur(), PhiX.accnos, $AUTHOR.concat.fasta,
	#		$AUTHOR.groups
	# OUTPUT: *.good.* *.bad.*

	screen_fasta=$AUTHOR.concat.fasta
	screen_groups=$AUTHOR.groups

	# Remove PhiX contamination
	if [[ -s PhiX_outs/PhiX.accnos ]]; then
		remove_params=("fasta=$screen_fasta,"\
			"group=$screen_groups,"\	# Not needed by mothur?
			"accnos=PhiX/PhiX.accnos")
		mothur <(mothurBatch remove seqs "${remove_params[*]}")
		rm -v phix_batch.txt
	else
		echo "PhiX.accnos is empty, skipping..."
	fi

	if [[ ! $screen_fasta ]] || [[ ! $screen_groups ]]; then
		echo "$0: ERROR: Zero-length input variables detected!"
		exit 1
	else
		screen_params=("fasta=$screen_fasta," "group=$screen_groups,"\
			"minlength=200," "maxlength=300," "maxambig=0,"\
			"maxhomop=8")
		summary_params=("fasta=current," "processors=$THREADS")
		mothurBatch screen seqs "${screen_params[*]}" > screening_batch.txt
		mothurBatch summary seqs "${summary_params[*]}" >> screening_batch.txt
		mothurBatch count groups "group=current" >> screening_batch.txt 
	
		echo "cat'ing screening_batch.txt"
		cat screening_batch.txt
		mothur screening_batch.txt 
		rm screening_batch.txt
	fi
#}

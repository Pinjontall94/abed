#!/bin/bash

#phixScreen(){
	# DESCRIPTION: Takes in a fasta file, phix databases in a directory named
	# 	$PWD/phix_db/ and screens for PhiX contamination

	# INPUT: bowtie2(), (trimmed_fastas/|offset_fastas/), $PHIX_DB
	# OUTPUT: *.screened.fasta, *.merged.PhiX, *.merged.bowtie
	# Deps: bowtie2
	
	# Set default PhiX db directory (if none provided)
	if [[ ! $PHIX_DB ]]; then
		PHIX_DB=phix_db/PhiX_bowtie_db
	fi

	# Loop bowtie2 over processed fastas
	echo "Screening for PhiX contamination"
	for i in ${in_file[@]}; do
		bowtie2 -f $i -x $PHIX_DB\
			-S ${i%%.*}.merged.bowtie\
			--un ${i%%.*}.screened.fasta\
			--al ${i%%.*}.merged.PhiX\
			--local -p $THREADS
	done

	# Once complete, generate list of reads to be screened out
	for i in *.merged.PhiX; do
		if [[ -s $i ]]; then
			echo "generating PhiX accession number file..."
			# Old grep command, probably won't work on multiple
			# 	files without generating "<filename>:<accno>"
			# grep ">" *merged.PhiX >> PhiX.accnos

			# [Do accnos in *.merged.PhiX start with ">"?]
			# 	also, this works with *.merged.PhiX,
			#	not just a single file "merged.PhiX"
			grep -Eo "[A-Z]{3,6}[0-9]+\.[0-9]+" *.merged.bowtie \
			| awk -F: '{ print $2 }' >> PhiX.accnos
		else
			rm $i
		fi
	done

	# Concatenate PhiX screened fastas to single file
	# (identified by $AUTHOR)
	cat *.screened.fasta > $AUTHOR.concat.fasta

	# Create output folder for screened fastas (separate from other PhiX outputs)
	mkdir screened_fastas && mv *.screened.fasta screened_fastas

	# Remove unneeded ".bowtie" outs
	rm *.merged.bowtie
#}

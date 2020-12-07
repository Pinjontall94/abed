#!/bin/bash

#mergeSeqs(){
	# merge all fastqs in the current directory with bbmerge.sh

	# INPUT: sufRegex(), bbmerge.sh(), fastqs/
	# OUTPUT: merged_fastqs/

	# Deps: bbmerge.sh

	echo "Checking forward and reverse fastq read counts"\
		"prior to merging"

	# Merge all _1 & _2 fastqs
	for i in ${in_file[@]}; do
		if [[ $i =~ _1.fastq ]]; then
			# Select forward reads in the for loop,
			#	then designate the reverse from $i
			REV=${i//_1.fastq/_2.fastq}
			READ_NAME=${i%_1.fastq}
			OUT_FILE=${i%_1.fastq}.merged.fastq
			# Make sure the read numbers are identical between
			#	forward and reverse
			echo "Checking: $i $REV"
			F_READ_NUM=$(grep -c "@$READ_NAME" $i)
			R_READ_NUM=$(grep -c "@$READ_NAME" $REV)
	
			if [[ $F_READ_NUM == $R_READ_NUM ]]; then
				echo "Read counts match; "\
					"Now merging: $i and $REV"
				echo "Merge out file: $OUT_FILE"
				bbmerge.sh in1="$i" in2="$REV" out="$OUT_FILE"
			else
				echo "Mismatch found, delete reads for $READ_NAME and "\
					"re-run fastq-dump before proceeding"
				exit 0
			fi
		elif [[ ! $i =~ _(1|2).fastq ]]; then # Matches reads w/o "_1" or "_2"
			echo "Renaming pre-merged read $i to ${i//fastq/merged.fastq}"
			mv -v $i ${i%.fastq}.merged.fastq
		else
			continue
		fi
	done

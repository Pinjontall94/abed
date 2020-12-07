#!/bin/bash

offsetTrim(){
	# Apply offset (if applicable)

	for i in ${in_file[@]}; do
		PRE_OFFSET=${i//trimmed/preoffset}
		OFFSET_TEMP=${i//trimmed/offset.temp}
		OUTPUT=${i//trimmed/offset}
		if [[ $FOR_OFFSET ]] && [[ $REV_OFFSET ]]; then
			echo "Trimming $FOR_OFFSET nucleotides"\
				"from 5' end of $i"
			cutadapt -u $FOR_OFFSET -o $OFFSET_TEMP $i

			echo "Trimming $REV_OFFSET nucleotides"\
				"from 3' of $OFFSET_TEMP"
			cutadapt -u -$REV_OFFSET -o $OUTPUT $OFFSET_TEMP
		elif [[ $FOR_OFFSET ]]; then
			echo "Trimming $FOR_OFFSET nucleotides"\
				"from 5' end of $i"
			cutadapt -u $FOR_OFFSET -o $OUTPUT $i
		else
			echo "Trimming $REV_OFFSET nucleotides"\
				"from 3' of $i"
			cutadapt -u -$REV_OFFSET -o $OUTPUT $i
		fi
	done

	# Remove temp files
	rm -v *.temp.fasta 
#}


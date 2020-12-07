#!/bin/bash

#trimLoop(){
	# Takes the primer patterns from command line,
	# inputs to cutadapt, and loops through all fasta # files in $(pwd)

	# INPUT: cutadapt(), $FORPRIME, $REVPRIME, *.merged.fasta
	# OUTPUT: *.trimmed.fasta
	
	echo "Trimming primers..."
	echo "FORPRIME = $FORPRIME \t REVPRIME = $REVPRIME"

	for i in ${in_file[@]}; do
		m=${i%.fasta}
		if [[ $FORPRIME ]]; then
			echo "FORPRIME = $FORPRIME"
			cutadapt -g $FORPRIME\
				-o $m.temp.fasta $i\
				--discard-untrimmed
		fi

		if [[ $REVPRIME ]]; then
			echo "REVPRIME = $REVPRIME"
			cutadapt -a $REVPRIME\
				-o $m.trimmed.fasta $m.temp.fasta\
				--discard-untrimmed
		fi
	done
	rm -v *.temp.fasta
#}

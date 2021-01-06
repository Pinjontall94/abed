#!/bin/bash

#fastaHeaderrelabel
	# Takes a fasta file labeled <author>_<year>_<accno>.fasta
	# 	and relabels headers for use with vsearch
	#
	# e.g.:
	# ">SRR10007909.1201 1201 length=251"
	# 		|
	# 		V
	# ">Li_2019_SRR10007909_1;barcodelabel=Li_2019_SRR10007909;"
		

	# INPUT:  $AUTHOR_*.fasta
	# OUTPUT: *.barcoded.fasta

OUTPUT=$1
shift

for i in $@; do
	awk -v NAME_STRIPPED="${i%.fasta}" '{
		# Count every header line
		if (/^>/) COUNT+=1 

		# Define new header format
		VSEARCH_HEADER=">"NAME_STRIPPED"_"COUNT"\;barcodelabel="NAME_STRIPPED"\;"

		# Relabel everything after ">" with new header
		gsub(/^>.*/, VSEARCH_HEADER)
		print;
	}' $i >> $OUTPUT		#${i%.fasta}.bar.fasta 2>/dev/null
done

# Concatenate output files to single barcoded fasta
# (also identified by $AUTHOR)
#cat *.bar.fasta > $AUTHOR.barcoded.fasta

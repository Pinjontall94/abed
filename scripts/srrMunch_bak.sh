#!/bin/bash

#srrMunch

# Short script to loop through SRRs in an SRA list, running fastq-dump, and
# 	nicely formatting them with bbtools' reformat.sh (for now, with manual
#	read number relabeling).

# Deps: fastq-dump from NCBI SRA Toolkit

# Usage: srrMunch <SraRunTable.txt> <SRR_Acc_List.txt> <Region Filter String>

# Check Run Table first for a list of 16S seqs
# TODO: check that the following snippet accurately grabs the accession
#    16S accession numbers from the correct column

# Determine Run Table delimiter (commas or tabs)
comma_count=$(grep -o "," $1 | wc -c)
tab_count=$(grep -o "\t" $1 | wc -c)
if ((comma_count > tab_count)); then
	# commas > tabs
	#Setting delimiter to commas
	DELIMITER=","
elif ((comma_count < tab_count)); then
	# commas < tabs
	# Setting delimiter to tabs
	DELIMITER="\t"
else
	echo "$0: Warning, Run Table not comma- or tab-separated"
	echo "Defaulting to whitespace..."
	DELIMITER=" "
fi

# Extract all accession numbers for a specified region (default: "16S")
if [[ $3 ]]; then
	Region="$3"
else
	Region="16S"
fi

awk -F "$DELIMITER" -v region="$Region" '
	# Loop over column 2, set column with valid accession number as target
	NR==2 {
		for(i=1;i<=NF;i++){
			if($i ~ "[A-Z]RR[0-9]{4,6}"){ target=i }
		}
	}

	# Print all entries in the target column if they contain the filter string
	# in their row (skipping the header row, of course)
	NR>1 && /$region/ { print $target }
	' $1 > tr_list.txt.temp

# Set which list to loop fastq-dump over
if [[ -s tr_list.txt.temp ]]; then
	echo "$0: awk filter successful, setting 16S_list.txt.temp to \$READ_LIST"
	READ_LIST=tr_list.txt.temp	
else
	echo "$0: awk filter failed to find 16S seqs,"\
		"leaving filtering duties to cutadapt..."
	READ_LIST=$2
fi

echo "Looping through $READ_LIST with sratoolkit and" \
	"unpacking fastqs..."

# Loop fastq-dump over Accession Number list
while read accno; do
	fastq-dump --split-3 --gzip $accno 2>>fastq-dump.err
done < $READ_LIST
gzip -d *.fastq.gz
[[ ! -f 16S_list.txt.temp ]] || rm -v 16S_list.txt.temp

# Parse fastq-dump.err files for fastq-dump download errors,
# store matching accession numbers in array
fastq_err=()
while IFS= read -r -d $'\0' line; do
	fastq_err+=("$line")
done < <(awk '/fastq-dump.*err:/ { print $NF }' fastq-dump.err)

# If fastq-dump error accnos are found, delete the messed up fastqs
# 	re-download, and re-decompress them
if [[ ${#fastq_err[@]} > 0 ]]; then
	for i in ${fastq_err[*]}; do
		rm $i*.fastq
		fastq-dump --split-3 --gzip $i
		gzip -d $accno*.fastq.gz
	done
else
	echo "No download errors found for fastq-dump"
	rm -v fastq-dump.err
fi

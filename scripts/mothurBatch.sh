#!/bin/bash

#mothurBatch
# Wrapper script for generating simple mothur commands, prints
# 	command to stdout to create batch files

# Note: If you have many parameters, put them in an array!
# Ex: EX_PARAMS+=("fasta=blah.fasta, group=blah.groups")

# Just make sure to refer to them with * expansion, not @ expansion,
#	and enclose in double quotes!
# e.g. mothurBatch doggo play "${EX_PARAMS[*]}"
# 	|==> doggo.play(fasta=blah.fasta, group=blah.groups)


COMMAND=$1.$2
PARAMS="$3"
echo "$COMMAND($PARAMS)"

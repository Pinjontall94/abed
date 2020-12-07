#!/bin/bash

joinBy(){
	# Takes first argument as new delimiter, then echoes all following
	# 	arguments with that delimiter
	# Usage: joinBy <separator> <arg1> [...] <argN>

	OIFS=$IFS
	IFS=$1

	shift
	echo "$*"

	IFS=$OIF
}

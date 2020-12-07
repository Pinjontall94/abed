#!/bin/bash

#garbageDay
	# Delete unneeded outputs and directories to clean up final output

	rm $AUTHOR.*.bad.*

	# Remove PhiX_outs entirely, if it's empty
	# (equivalent to "rmdir PhiX_outs" without unnecessary errors)
	find PhiX_outs -maxdepth 0 -empty -exec rmdir {} \;

	# Move all remaining author files to outputs dir
	mkdir asdMetagen_final && mv $AUTHOR.* asdMetagen_final

	# Move logs to log dir
	mkdir asdMetagen_logs && mv *.log* asdMetagen_logs

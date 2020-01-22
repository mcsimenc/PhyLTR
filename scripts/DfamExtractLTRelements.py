#!/usr/bin/env python3

import sys


# this script takes the Dfam multi-nhmm as input and outputs a multi-nhmm file
# containing only hmm profiles which are annotated as being in the LTR class
hmmlines = ''
gethmm = False
for line in sys.stdin:
    # recognize line with annotation information and set the variable hmmlines
    # to true, which signals to output lines of the current profile after they
    # have all been collected
	if line.startswith('CT'):
		if 'Class;' in line:
			if 'LTR;' in line:
				hmmlines += line
				gethmm = True
			else:
				hmmlines = ''
				gethmm = False
	elif line.startswith('//'):
        # at the start of the next profile, print the collected lines if the
        # profile has a LTR annotation, otherwise reset the variable hmmlines
        # to prepare to collect the next profile's lines
		if gethmm:
			hmmlines += line
			print(hmmlines, end='')
			gethmm = False
		hmmlines = ''
    # start collecting lines from the start of a profile until the CT line
	else:
		hmmlines += line

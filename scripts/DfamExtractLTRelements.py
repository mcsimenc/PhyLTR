#!/usr/bin/env python3

import sys

hmmlines = ''
GETHMM = False
for line in sys.stdin:
	if line.startswith('CT'):
		if 'Class;' in line:
			if 'LTR;' in line:
				hmmlines += line
				GETHMM = True
			else:
				hmmlines = ''
				GETHMM = False
	elif line.startswith('//'):
		if GETHMM:
			hmmlines += line
			print(hmmlines, end='')
			GETHMM = False
		hmmlines = ''
	else:
		hmmlines += line

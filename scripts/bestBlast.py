#!/usr/bin/env python3

# This script takes a blast output format 7 table on stdin and prints the first (presumably highest scoring/lowest evalue) hit for each query on stdout

import sys

if '-h' in sys.argv:
	print('''
 usage:
	python3 bestBlast.py < blastTabularOutput

 description:
	this script outputs the identity of the first match for each
   	query in the output (which is the highest scoring match in
	blast outfmt 6 or 7) and prints it.
''', file=sys.stderr)
	sys.exit()

last_query = '!@#--sjsjsj-BI29(*YFen230*012kr3jnr)FU)R(u$%^&*()_QWERTYUIOP{'
for line in sys.stdin:
	if not line.startswith('#') and not line.startswith('{0}\t'.format(last_query)):
		fields = line.strip().split('\t')
		print('\t'.join(fields[:2]))
		last_query = line.strip().split()[0]

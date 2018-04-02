#!/usr/bin/env python3


# This script takes a blast output format 7 table on stdin and prints the first (presumably highest scoring/lowest evalue) hit for each query on stdout


import sys



last_query = '!@#--sjsjsj-BI29(*YFen230*012kr3jnr)FU)R(u$%^&*()_QWERTYUIOP{'
for line in sys.stdin:
	if not line.startswith('#') and not line.startswith('{0}\t'.format(last_query)):
		fields = line.strip().split('\t')
		if fields[12] == 'N/A':
			fields[12] = '?'
		print('\t'.join(fields[:2]+[fields[12]]))
		last_query = line.strip().split()[0]


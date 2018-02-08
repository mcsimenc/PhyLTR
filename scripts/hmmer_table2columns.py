#!/usr/bin/env python3

import sys

if '-h' in sys.argv:
	print('''
		 usage:
			HMMER.tbl.parser.py < input.tbl > output.tbl

		keeps best hit, within e-value threshold. Default 1e-5
		''', file=sys.stderr)
	sys.exit()

def update_dct(line, dct={}, key_field=1, value_field=3, evalue_field=5, score_field=6, evalue_cutoff=1e-5):

	'''
	Adds info from line to dict. Keeps key with highest score, within e-value thresold. Defaults for hmmsearch --tblout
	'''

	fields = line.strip().split()

	if float(fields[int(evalue_field)-1]) < float(evalue_cutoff):
		if fields[int(key_field) -1] in dct:
			if float(fields[int(score_field)-1]) > dct[fields[int(key_field) -1]][1]:
				dct[fields[int(key_field) -1]] = (fields[int(value_field)-1], float(fields[int(score_field)-1]))

		else:
			dct[fields[int(key_field) -1]] = (fields[int(value_field)-1], float(fields[int(score_field)-1]))

	return dct


dct = {}

for line in sys.stdin:
	if not line.startswith('#'):

		dct = update_dct(line, dct)

for key in sorted(list(dct.keys())):

	print('{0}\t{1}'.format(key, dct[key][0]))
##                                                                     --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
## target name             accession  query name           accession     E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
##     ------------------- ---------- --------------------  ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
#RM_24157_rnd-6_family-183  -          5S                   DF0000012.4   2.6e-41  141.2   0.0   2.6e-12   47.4   0.0   5.1   6   0   0   6   6   6   3 -
#RM_24151_rnd-6_family-627  -          5S                   DF0000012.4   4.2e-31  108.1   0.0   5.5e-11   43.1   0.0   5.6   7   0   0   7   7   7   3 -
#RM_24152_rnd-6_family-141  -          5S                   DF0000012.4   6.2e-24   85.0   0.0     1e-07   32.6   0.0   4.5   4   0   0   4   4   4   4 -
#RM_24160_rnd-5_family-1629 -          5S                   DF0000012.4   9.1e-19   68.3   0.0    0.0048   17.5   0.0   4.5   4   0   0   4   4   4   4 -
#RM_24145_rnd-5_family-1467 -          5S                   DF0000012.4   6.9e-18   65.4   0.0    0.0076   16.8   0.0   4.4   4   0   0   4   4   4   4 -
#RM_24144_rnd-5_family-1209 -          5S                   DF0000012.4   9.5e-18   65.0   0.0     0.015   15.9   0.0   4.5   4   0   0   4   4   4   4 -
#RM_24138_rnd-5_family-1432 -          5S                   DF0000012.4   1.7e-17   64.1   0.0     0.029   14.9   0.0   4.5   4   0   0   4   4   4   4 -

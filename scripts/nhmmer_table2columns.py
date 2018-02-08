#!/usr/bin/env python3

import sys

if '-h' in sys.argv:
	print('''
		 usage:
			nhmmer_table2columns.py < input.tbl > output.tbl

		keeps best hit, within e-value threshold. Default 1e-5
		''', file=sys.stderr)
	sys.exit()

def update_dct(line, dct={}, key_field=1, value_field=3, evalue_field=13, score_field=14, evalue_cutoff=1e-5, strand_field=12):

	'''
	Adds info from line to dict. Keeps key with highest score, within e-value thresold. Defaults for hmmsearch --tblout
	'''

	fields = line.strip().split()

	if float(fields[int(evalue_field)-1]) < float(evalue_cutoff):
		if fields[int(key_field) -1] in dct:
			if float(fields[int(score_field)-1]) > dct[fields[int(key_field) -1]][1]:
				dct[fields[int(key_field) -1]] = (fields[int(value_field)-1], float(fields[int(score_field)-1]), fields[int(strand_field)-1])

		else:
			dct[fields[int(key_field) -1]] = (fields[int(value_field)-1], float(fields[int(score_field)-1]), fields[int(strand_field)-1])

	return dct


dct = {}

for line in sys.stdin:
	if not line.startswith('#'):

		dct = update_dct(line, dct)

for key in sorted(list(dct.keys())):

	print('{0}\t{1}\t{2}'.format(key, dct[key][0], dct[key][2]))

## nhmmer output
##
## target name            accession  query name           accession   hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
##    ------------------- ---------- --------------------  ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
#LTR_retrotransposon3327  -          ACCORD2_I            DF0001530.0    4416    5277    4041    4905    4021    4926    7028    +       9e-24   82.6  34.4  -
#LTR_retrotransposon10954 -          ACCORD2_I            DF0001530.0    4416    5277    2992    2128    3012    2108    6655    -     7.4e-23   79.6  32.7  -
#LTR_retrotransposon3256  -          ACCORD2_I            DF0001530.0    4246    5277    3965    4984    3945    5004    7108    +     1.7e-22   78.4  35.8  -
#LTR_retrotransposon12749 -          ACCORD2_I            DF0001530.0    4416    5277    2988    2124    3008    2104    7109    -     2.1e-22   78.0  34.2  -
#LTR_retrotransposon4233  -          ACCORD2_I            DF0001530.0    4416    5276    4121    4984    4101    5005    7113    +     2.5e-22   77.8  35.6  -
#LTR_retrotransposon5179  -          ACCORD2_I            DF0001530.0    4416    5277    3190    2326    3210    2306   14241    -     2.7e-22   77.7  34.1  -
#LTR_retrotransposon11608 -          ACCORD2_I            DF0001530.0    4416    5277    2979    2115    2999    2095    6935    -     2.8e-22   77.6  35.5  -
#LTR_retrotransposon13259 -          ACCORD2_I            DF0001530.0    4416    5277    4094    4958    4074    4978    7082    +       3e-22   77.5  32.8  -

#!/usr/bin/env python3

import sys

best_orfs = {}
query = None
for line in sys.stdin:
	if line.startswith('LTR'):
		query, subject, pid, alnlen, mis, gapo, qstart, qend, sstart, send, evalue, bitscore = line.strip().split('\t')
		alnlen = int(alnlen)
		if query in best_orfs:
			best_orfs[query].append([alnlen])
		else:
			best_orfs[query] = [[alnlen]]
	elif line.startswith('>'):
		desc = ' '.join(line.strip().split(' ')[1:])
		best_orfs[query][-1].append(desc)

for orf in best_orfs:
	best_orfs[orf].sort(reverse=True, key=lambda x:x[0])
	i = 0
	while i < 6:
		print('{0}\t{1}\t{2}'.format(orf, best_orfs[orf][i][0], best_orfs[orf][i][1]))
		i += 1
	print('###')
	print('###')


## Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#LTR_retrotransposon2709.ORF.02	XP_013727080.1	61.97	71	27	0	27	97	20	90	6e-23	102
#>XP_013727080.1 PREDICTED: uncharacterized protein LOC106430829 [Brassica napus]
#
## Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#LTR_retrotransposon2709.ORF.02	EMS66121.1	75.68	37	9	0	33	69	76	112	5e-11	68.2
#>EMS66121.1 hypothetical protein TRIUR3_20723 [Triticum urartu]

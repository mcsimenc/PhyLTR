#!/usr/bin/env python3

import sys


def help():
	print('''
	usage:
		blast2nrheaders.py -nr <path> -blast <path> -min_pid <int|float>

		-nr		Path to nr
		-blast		Path to blast results 
		-min_pid	Minimum percent id to include in output. default 60.0
	''', file=sys.stderr)


args = sys.argv

if '-nr' not in args or '-blast' not in args or len(args) < 5:
	help()
	sys.exit()

nr = args[args.index('-nr')+1]
blast = args[args.index('-blast')+1]
if '-min_pid' in args:
	min_pid = float(args[args.index('-min_pid')+1])
else:
	min_pid = 60.0

orfs = {}
with open(blast, 'r') as blastFl:
	for line in blastFl:
		if line.startswith('#'):
			continue
		query, subject, pid, alnlen, mis, gapo, qstart, qend, sstart, send, evalue, bitscore = line.strip().split('\t')
		pid = float(pid)
		if pid >= min_pid:
			if query in orfs:
				orfs[query][subject] = {'blastresult':line.strip()}
			else:
				orfs[query] = {subject:{'blastresult':line.strip()}}
		
with open(nr, 'r') as nrFl:
	for line in nrFl:
		if line.startswith('>'):
			acc = line.split(' ')[0][1:]
			for orf in orfs:
				if acc in orfs[orf]:
					orfs[orf][acc]['nrdesc'] = line.strip()
outFlPth = '{0}.with_nr_descriptions.tab'
fields = '# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n'
with open(outFlPth, 'w') as outFl:
	for orf in sorted(list(orfs.keys())):
		for hit in orfs[orf]:
			outFl.write(fields)
			outFl.write(orfs[orf][hit]['blastresult']+'\n')
			outFl.write(orfs[orf][hit]['nrdesc']+'\n\n')

print('Done.', file=sys.stderr)

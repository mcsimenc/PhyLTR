#!/usr/bin/env python3

# Takes a two-column list of scaffold names and lengths and outputs a ideogram file for circos

import sys

if '-h' in sys.argv:
	print('''
			usage: ideogramFromLengths.py < lengths.txt > karyotype.circos

			
			description: takes a two-col input file where first col is a seq name
				     and second col is the seq length and outputs a karyotype
				     file for circos.
		''')
	sys.exit()

counter = 1
for line in sys.stdin:
	line=line.strip().split('\t')
	if counter%2==0:
		color = 'greys-6-seq-2'
	else:
		color = 'greys-6-seq-4'
	#chrNum = int(line[0].split('s')[1])
	chrNum = counter
	print('chr - {0}\t{1}\t0\t{2}\t{3}'.format(line[0], chrNum, line[1], color))
	counter += 1

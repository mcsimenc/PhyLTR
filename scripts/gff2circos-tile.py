#!/usr/bin/env python3

# If you have questions about this script or something doesn't work right you can email Matt at mcsimenc@gmail.com

import sys


def help():
	print('''

		Usage:
		------------
		gff2circos-tile.py -gff <file> [-scafList <file> -valueDef <file> ]

		Description:
		------------
		Takes a GFF3 file as input and outputs a Circos data file for a tile
		track. The fourth column in the output (value) tells circos how to color
		the feature just like it does with heatmap tracks. The default value is 0
		but a two-column file can be provided using -valueDef to customize values.


		Required parameters:
		------------
		-gff	  	<path>	Path to input maker gff3 file (mandatory).

		#-scafList 	<path>	Limit the number of scaffolds printed out. NOT IMPLEMENTED

		-valueDef 	<path>	Tab-delimited file where first col is a string
					to search for in each GFF3 line and the second
					col is the value to assign if the string is found.
					The second col needs to contain numbers only.

		Output:
		------------
		scaf	start	stop	value

''', file=sys.stderr)

def gff2circosTileTrack(gffFl, valueDef=None):
	with open(gffFl) as fl:
		for line in fl:
			if not line.startswith('#'):
				fields = line.strip().split('\t')
				seq = fields[0]
				start = fields[3]
				end = fields[4]
				value = 0
				if not valueDef == None:
					FOUND_DEF = False
					MatchValue = 0
					for definition in valueDef:
						if definition[0] in line:
							value = definition[1]
							if FOUND_DEF == True and value != MatchValue:
								print('Found more than one definition with different values for GFF3 line:\n{0}'.format(line.strip()), file=sys.stderr)
							FOUND_DEF = True
							MatchValue = value

				print('{0}\t{1}\t{2}\t{3}'.format(seq, start, end, value))

				

args = sys.argv

if '-h' in args or len(args) < 3 or '-gff' not in args:
	help()
	sys.exit()

gff_filepath = args[args.index('-gff') +1]

#if '-scafList' in args:
#	scafListFl = args[args.index('-scafList') +1]
#	scafList = open(scafListFl).read().strip().split('\n')
#else:
#	scafList = None


if '-valueDef' in args:
	valueDefFl = args[args.index('-valueDef') +1]
	valueDef = set([ (definition.split('\t')[0], float(definition.split('\t')[1])) for definition in open(valueDefFl).read().strip().split('\n') ])
else:
	valueDef = None


gff2circosTileTrack(gff_filepath, valueDef)

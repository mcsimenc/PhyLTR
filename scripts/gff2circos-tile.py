#!/usr/bin/env python3

# If you have questions about this script or something doesn't work right you can email Matt at mcsimenc@gmail.com

import sys

class GFF3_line:
	'''
	Attributes:
			field0, ... , field8 - string
			attributes - dictionary

	Methods:    str()		prints gff line
		    refreshAttrStr()	updates gff line attributes (needed if changes to attributes were made)

	kwargs is a dictionary
	'''
	def __init__(self, line=None, **kwargs):

		if line == None:
			(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes_str) = [None]*9
			self.attributes_order = []
			self.attributes = {}
			
		else:
			(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes_str) = line.strip().split('\t')
			self.start = int(self.start)
			self.end = int(self.end)
			attributes_list = self.attributes_str.split(';')
			self.attributes_order = [ attr.split('=')[0] for attr in attributes_list ]
			self.attributes = { attr.split('=')[0]:attr.split('=')[1] for attr in attributes_list }

		self.line_number = None

		if 'line_number' in kwargs:	
			self.line_number = kwargs['line_number']

		# rename the name attribute so it conforms to GFF3 specifications, where Name is a reserved attribute key. The version of LTRDigest makes attribute key name
		if 'name' in self.attributes:
			self.attributes['Name'] = self.attributes.pop('name')
			self.attributes_order[self.attributes_order.index('name')] = 'Name'

	def __repr__(self): # for when str() or repr() are called

		return '\t'.join( [ str(self.seqid), str(self.source), str(self.type), str(self.start), str(self.end), str(self.score), str(self.strand), str(self.phase), str(self.attributes_str) ] )

	def refreshAttrStr(self, attrOrder=None):
		'''
		If the attributes have been changed this needs to be called
		before the changes are reflected on a str() call on this object.

		If an attribute has been added it should have also been added to self.attributes_order
		'''
		self.attributes_str = ';'.join([ '='.join([ attr, self.attributes[attr] ]) for attr in self.attributes_order ])



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

					If valueDef == 'LTR' then the ID attribute is expected to be
						LTR_retrotransposonN where N is some integer. N is used
						as the value.

		Output:
		------------
		scaf	start	stop	value

''', file=sys.stderr)
def gff2circosTileTrack(gffFl, valueDef=None):


	with open(gffFl) as fl:
		for line in fl:
			if not line.startswith('#'):
				if valueDef=='LTR':
					gffLine = GFF3_line(line)
					start = gffLine.start
					end = gffLine.end
					value = gffLine.attributes['ID'][19:]

				else:
					fields = line.strip().split('\t')
					seq = fields[0]
					start = fields[3]
					end = fields[4]
					value = 0

				if not valueDef == None and not valueDef == 'LTR':
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
	if valueDefFl == 'LTR':
		valueDef = 'LTR'
	else:
		valueDef = set([ (definition.split('\t')[0], float(definition.split('\t')[1])) for definition in open(valueDefFl).read().strip().split('\n') ])
else:
	valueDef = None


gff2circosTileTrack(gff_filepath, valueDef)

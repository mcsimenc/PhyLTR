#!/usr/bin/env python3

import sys

class GFF3_line:
	'''
	Attributes:
			seqid, source, type, start, end, score, strand, phase, attributes

			attributes_order: A list declaring the order of attributes desired.
					  If new attributes are added to object, they must also
					  be added to attributs_order and then run refreshAttrStr()

			kwargs: line_number

	Methods:
			str()			Prints gff line
			
			refreshAttrStr()	Updates gff line attributes (needed if changes to attributes were made)
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

for line in sys.stdin:
	if not line.startswith('#'):
		l = GFF3_line(line)
		if l.type == 'ORF':
			print('>{}'.format(l.attributes['ID']))
			print(l.attributes['translated_seq'])

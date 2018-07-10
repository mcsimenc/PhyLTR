#!/usr/bin/env python3

import sys


def help():
	print('''
 Usage:
 ------------
 addGOandIPS.py <domains> <go_ips> <mapped_descs> <membership>
 
 Description:
 ------------
 Puts GO and IPS terms and descriptions of the Blast2GO mapped protein name into
 a single TSV file.

 Input:
 -----------
 <domains>	2-col tsv	domain	geneName
 <go_ips>	GFF3		with attr. (col 9) of: GO_terms and InterProScan_terms
 <mapped_descs>	2-col tsv	geneName	desc
 <membership>	3-col tsv	LTRRTnum	classification	cluster
 
 ''')


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

args=sys.argv

if len(args) < 5:
	help()
	sys.exit()


IPS = {}
GO = {}
with open(args[2],'r') as inFl:
	for line in inFl:
		if not line.startswith('#'):
			gl = GFF3_line(line)
			if gl.type == 'gene':
				g = gl.attributes['ID']
				go = gl.attributes['GO_terms']
				ips = gl.attributes['InterProScan_terms']
				
				if g in IPS:
					sys.exit('two lines for one gene: {0}'.format(g))
				if g in GO:
					sys.exit('two lines for one gene: {0}'.format(g))
				
				GO[g] = go
				IPS[g] = ips

			
DESC = {}
with open(args[3], 'r') as inFl:
	for line in inFl:
		g, d = line.strip().split('\t')
		DESC[g] = d

#Azfi_s4729	maker	gene	5003	7129	.	-	.	ID=Azfi_s4729.g121372;Name=Azfi_s4729.g121372;GO_terms=None;InterProScan_terms=noIPSmatch

INFO = {}
with open(args[4], 'r') as inFl:
	for line in inFl:
		e, classif, clust = line.strip().split()
		e = 'LTR_retrotransposon'+e
		INFO[e] = [classif, clust]

print('gene\tblast2goMapped\tGOterms\tIPSterms\tdomainHit\telement\tclassification\tcluster')
with open(args[1], 'r') as inFl:
	for line in inFl:
		contents = line.strip().split()
		orf, g = contents[:2]
		if g not in DESC:
			DESC[g] = 'None'
		el = orf.split('.')[0]
		try:
			INFO[el]
			print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(g, DESC[g], GO[g], IPS[g], orf, el, INFO[el][0], INFO[el][1]))
		except KeyError:
			print('{0} not found in cluster membership file.'.format(el))

#element	classification	cluster
#11113	Ngaro	0
#12956	Ngaro	0
#7768	Ngaro	0
#13094	Ngaro	1
#14761	Ngaro	1
#5903	Ngaro	2
#16057	Ngaro	2
#10498	Ngaro	3
#12561	Ngaro	4
#Azfi_s0001.g000026	RVT_1
#Azfi_s0001.g000029	ORF
#Azfi_s0001.g000030	ORF
#Azfi_s0001.g000076	ORF
#Azfi_s0001.g000116	ORF
#Azfi_s0001.g000141	ORF
#Azfi_s0001.g000141	zf-RVT
#Azfi_s0001.g000146	ORF
#Azfi_s0001.g000205	ORF
#Azfi_s0001.g000210	ORF


#!/usr/bin/env python3

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


args = sys.argv

def geneconv2circoslinks(geneconvfile, ltrharvestgff, outfile):
	'''
	Converts GI tract pairs from geneconvClusters() output and writes a links file for Circos
	The GFF3 is needed to get the scaffold name.
	seqlengths needs to be a dictionary with the lengths of the sequences whose names correspond
	to the sequence names in the gff for the features with gene conversion tracts.
	Assumes LTR_retrotransposon features were used.
	'''
	global paths

	#lengths = { l[0]:int(l[1]) for  l in open(seqlengths).read().strip().split('\n')}

	seqs = {}
	starts = {}
	with open(ltrharvestgff, 'r') as inFl:
		for line in inFl:
			if line.startswith('#'):
				continue
			if '\tLTR_retrotransposon\t' in line:
				gffLine = GFF3_line(line)
				element = gffLine.attributes['ID']
				scaf = gffLine.seqid
				seqs[element] = scaf
				starts[element] = int(gffLine.start)
	
	with open(outfile, 'w') as outFl:
		with open(geneconvfile, 'r') as inFl:
			for line in inFl:
				if line.startswith('GI'):
					rec = line.strip().split('\t')
					el1, el2 = [ 'LTR_retrotransposon{0}'.format(e[1:]) for e in rec[1].split(';') ]
					el1start  = int(rec[7]) + starts[el1] - 1
					el1end  = int(rec[8]) + starts[el1] - 1
					el2start  = int(rec[10]) + starts[el2] - 1
					el2end  = int(rec[11]) + starts[el2] - 1
					el1seq = seqs[el1]
					el2seq = seqs[el2]
					outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end))


geneconv2circoslinks(args[1], args[2])

#hs1 100 200 hs2 250 300
#hs1 400 550 hs3 500 750
#hs1 600 800 hs4 150 350

## Frag	Sequence	Sim	BC KA		Aligned Offsets			In Seq1			In Seq2		Num	Num	Total	Mismatch	
##   type	Names	Pvalue	Pvalue	Begin	End	Length	Begin	End	Length	Begin	End	Length	Poly	Dif	Diffs	Penalty	
#GI	S554;S877	0.0016	> 1.0	be	en	le	be	en	len	be	en	le	pol	nd	ndd	pen	Other	0
#GI	S1356;S999	0.0016	> 1.0	43	62	20	31	50	20	41	60	20	20	0	10	None	Other	0
#GI	S600;S860	0.0018	> 1.0	45	60	16	18	33	16	43	58	16	16	0	12	None	Other	0
#GI	S1106;S999	0.0018	> 1.0	41	62	22	39	60	22	39	60	22	22	0	9	None	Other	0
#GI	S1380;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
#GI	S1433;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
#GI	S1515;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
#GI	S521;S999	0.0018	> 1.0	41	62	22	27	48	22	39	60	22	22	0	9	None	Other	0

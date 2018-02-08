#!/usr/bin/env python3

import sys


def help():
	print('''
			usage:
				gffAddAttr.py -gff <filepath> -attr <string> -map <filepath> -mapKey <string> > output.gff


			description: Adds a new attribute to field 9 of input gff. 

			-attr	<string>	New attribute name. Throws an error if the name is taken.

			-gff	<filepath>	Input GFF3 file

			-map	<filepath>	Two column tab delimited file where first field corresponds
						to the value of the attribute specified by -mapKey and the second
						to the value for the new attribute.

			-mapKey	<string>	The attribute for which the first field of the -map file is derived.

			-restrictType <string>	Only add new attribute to features of this type

			-replace		Replace existing attribute if present

			-replaceIfNone		Replace existing attribute if present and is 'None'

			-v			Verbose reporting on stderr

		''', file=sys.stderr)
	sys.exit()

class GFF3_line:
	'''
	Attributes:
			field0, ... , field8 - string
			attributes - dictionary

	Methods:    str()		prints gff line
		    refreshAttrStr()	updates gff line attributes (needed if changes to attributes were made)

	kwargs is a dictionary
	'''
	def __init__(self, line, **kwargs):

		(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes_str) = line.strip().split('\t')
		self.line_number = None

		if 'line_number' in kwargs:	
			self.line_number = kwargs['line_number']

		attributes_list = self.attributes_str.split(';')
		self.attributes_order = [ attr.split('=')[0] for attr in attributes_list ]
		self.attributes = { attr.split('=')[0]:attr.split('=')[1] for attr in attributes_list }

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

args =sys.argv

if len(sys.argv) < 9 or '-h' in args:
	help()

attrMap = {}
with open(sys.argv[sys.argv.index('-map') + 1]) as mapFl:
	for line in mapFl:
		try:
			contents = line.strip().split('\t')
			gene = contents[0]
			annot = contents[1]
			attrMap[gene] = annot
		except:
			continue

#for k in attrMap:
#	print('key:\t{0}\nvalue:\t{1}'.format(k, attrMap[k]))
#sys.exit()

mapKey = sys.argv[sys.argv.index('-mapKey')+1]
newAttrKey = sys.argv[sys.argv.index('-attr')+1]
gffFilepath = sys.argv[sys.argv.index('-gff')+1]
if '-restrictType' in args:
	restrictType = sys.argv[sys.argv.index('-restrictType')+1]


with open(gffFilepath) as gffFl:
	for line in gffFl:

		if line.startswith('#'):
			print(line, end='')
		else:
			
			gffLine = GFF3_line(line)

			if '-restrictType' in args:
				if gffLine.type != restrictType:
					print(line, end='')
					continue

			if newAttrKey in gffLine.attributes:
				if not '-replace' in args and ('-replaceIfNone' in args and gffLine.attributes[newAttrKey] != 'None'):
					if '-v' in args:
						print('{0}\nAbove line in GFF3 input already has attribute {1}. Continuing'.format(str(gffLine), newAttrKey), file=sys.stderr)
					print(line, end='')
					continue

			try:
				featureIdentifier = gffLine.attributes[mapKey]

			except KeyError:
				if '-v' in args:
					print('{0}\n-mapKey {1} not in above GFF3 line.'.format(str(gffLine), mapKey), file=sys.stderr)
					print('--------', file=sys.stderr)

				newAttrValue = "None"
				gffLine.attributes[newAttrKey] = newAttrValue
				if newAttrKey not in gffLine.attributes_order:
					gffLine.attributes_order.append(newAttrKey)
				try:
					gffLine.refreshAttrStr()
				except:
					print(gffLine.attributes, file=sys.stderr)
					print(gffLine.attributes_order, file=sys.stderr)
					sys.exit()
				print(str(gffLine))
				continue

			if featureIdentifier in attrMap:
				newAttrValue = attrMap[featureIdentifier]
				gffLine.attributes[newAttrKey] = newAttrValue
				if newAttrKey not in gffLine.attributes_order:
					gffLine.attributes_order.append(newAttrKey)
				gffLine.refreshAttrStr()
				print(str(gffLine))
			else:
				if '-v' in args:
					print('{0}\nNo item in -map matching {1} from the above GFF3 line'.format(str(gffLine), featureIdentifier), file=sys.stderr)
					print('--------', file=sys.stderr)

				newAttrValue = "None"
				gffLine.attributes[newAttrKey] = newAttrValue
				if newAttrKey not in gffLine.attributes_order:
					gffLine.attributes_order.append(newAttrKey)
				gffLine.refreshAttrStr()
				print(str(gffLine))



		



#Acas.4051_0001	maker	gene	12849	16043	.	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12
#Acas.4051_0001	maker	mRNA	12849	16043	3195	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1;_AED=0.35;_eAED=0.35;_QI=0|-1|0|1|-1|0|1|0|1064
#Acas.4051_0001	maker	exon	12849	16043	.	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1:exon:0;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1
#Acas.4051_0001	maker	CDS	12849	16043	.	+	0	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1:cds;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.12-mRNA-1
#Acas.4051_0001	maker	gene	20802	23126	.	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8
#Acas.4051_0001	maker	mRNA	20802	23126	2325	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1;_AED=0.39;_eAED=0.40;_QI=0|-1|0|1|-1|0|1|0|774
#Acas.4051_0001	maker	exon	20802	23126	.	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1:exon:1;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1
#Acas.4051_0001	maker	CDS	20802	23126	.	+	0	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1:cds;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.8-mRNA-1
#Acas.4051_0001	maker	gene	16070	18857	.	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.6;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.6
#Acas.4051_0001	maker	mRNA	16070	18857	2040	+	.	ID=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.6-mRNA-1;Parent=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.6;Name=maker-Acas.4051_0001-exonerate_protein2genome-gene-0.6-mRNA-1;_AED=0.41;_eAED=0.41;_QI=0|0|0|1|0|0|2|0|679

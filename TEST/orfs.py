#!/usr/bin/env python3

import sys
import subprocess
from Bio import SeqIO

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


def Overlaps(j, k):
	'''
	Inputs, j and k, are tuples/lists with a start and a end coordinate as in: (start, end)
	If they overlap True is returned, if they do not overlap, False is returned.
	'''
	j = sorted([int(i) for i in j])
	k = sorted([int(i) for i in k])
	jk = sorted([j,k], key=lambda x:x[0])
	if jk[0][1] >= jk[1][0]:
		return True
	else:
		return False



def bestORFs(fasta, outdir, gff, minLen=30):
	'''
	Finds all ORFs in fasta using EMBOSS getorf and writes the best ones to a GFF3 and protein FASTA.
	The best ORFs are the set of non-overlapping ORFs containing the longest ORF out of sets on that
	strand,
	given by records in gff, or, if strand is unknown, the best set out of all.

	The coordinates output by getorf are 1-based

	Only ORFs with nucleotide sequences longer than minLen are kept. (default 30 bp)
	'''
	# Get strandedness for each element into a dictionary
	strands = {}
	with open(gff, 'r') as gffFl:
		for line in gffFl:
			if '\tLTR_retrotransposon\t' in line:
				gffLine = GFF3_line(line)
				el = gffLine.attributes['ID']
				strand = gffLine.strand
				strands[el] = strand

	# Run EMBOSS getorf
	outseq = '{0}/{1}.orfs'.format(outdir, fasta.split('/')[-1])
	#getorf_call = [ executables['getorf'], '-sequence', fasta, '-outseq', outseq ]
	#makecall(getorf_call)
	subprocess.call('/home/derstudent/software/EMBOSS-6.6.0/bin/getorf -sequence {0} -outseq {1}'.format(fasta, outseq), shell=True)

	# Read getorf output and put info into dicts
	orfs = list(SeqIO.parse(outseq, 'fasta'))
	orfs_ordered_lengths = {}
	orfs_ordered_coords = {}
	orfs_seqs_dct = {}
	orfs_coords = {}
	coords2lenkey = {}
	for orf in orfs:
		element = '_'.join(orf.id.split('_')[:-1])
		strand_gff = strands[element]
		orfnum = orf.id.split('_')[-1]
		desc = orf.description
		seq = str(orf.seq)
		if 'REVERSE' in desc:
			if strand_gff == '+': # don't consider reverse strand ORFs if the element is already assigned to the forward strand.
				continue
			strand = '-'
			start, end = sorted([int(desc.split()[1][1:]), int(desc.split()[3][:-1])])
			length = end - start + 1
			if length < minLen:
				continue
		else:
			if strand_gff == '-': # don't consider forward strand ORFs if the element is already assigned to the reverse strand.
				continue
			strand = '+'
			start, end = sorted([int(desc.split()[1][1:]), int(desc.split()[3][:-1])])
			length = end - start + 1
			if length < minLen:
				continue

		#print('strand:\t{0}'.format(strand))
		#print('desc:\t{0}'.format(desc))
		#print('orfnum:\t{0}'.format(orfnum))
		#print('element:\t{0}'.format(element))

		if element in orfs_coords:
			if strand in orfs_coords[element]:
				orfs_coords[element][strand][orfnum] = (start, end) # name, length
			else:
				orfs_coords[element][strand] = {orfnum:(start,end)} # name, length
		else:
			orfs_coords[element] = {strand:{orfnum:(start,end)}} # name, length

		if element in coords2lenkey:
			if strand in coords2lenkey[element]:
				coords2lenkey[element][strand][(start, end)] = (orfnum, length) # name, length
			else:
				coords2lenkey[element][strand] = {(start,end):(orfnum, length)} # name, length
		else:
			coords2lenkey[element] = {strand:{(start,end):(orfnum, length)}} # name, length

		if element in orfs_ordered_lengths:
			if strand in orfs_ordered_lengths[element]:
				orfs_ordered_lengths[element][strand].append((orfnum, length)) # name, length
			else:
				orfs_ordered_lengths[element][strand] = [(orfnum, length)] # name, length
		else:
			orfs_ordered_lengths[element] = {strand:[(orfnum, length)]} # name, length

		if element in orfs_ordered_coords:
			if strand in orfs_ordered_coords[element]:
				orfs_ordered_coords[element][strand].append((start, end)) # name, length
			else:
				orfs_ordered_coords[element][strand] = [(start, end)] # name, length
		else:
			orfs_ordered_coords[element] = {strand:[(start, end)]} # name, length

		if element in orfs_seqs_dct:
			if strand in orfs_seqs_dct[element]:
				orfs_seqs_dct[element][strand][orfnum] = seq # name, length
			else:
				orfs_seqs_dct[element][strand] = {orfnum:seq}
		else:
			orfs_seqs_dct[element] = {strand:{orfnum:seq}}

	# For each element, find ORFs for its strand or both strands if strand = ?

	# Order ORFs names in list from element/strand longest to shortest.

	# Order another ORF name list with coordinate position of starts, smallest to largest

	# Make a dict with the sequences to use after selecting which orfs to keep

	sorted_elements = sorted(list(orfs_ordered_lengths.keys()))
	for element in sorted_elements:
		strand = strands[element]
		if strand == '?' or strand == '.':
			strand = ['+', '-']
		else:
			strand = [strand]
		best_orf_sets = {'+':None, '-':None}
		for s in strand:
			i = 0
			orfs_ordered_lengths[element][s].sort(reverse=True, key=lambda x:x[1])
			orfs_ordered_coords[element][s].sort(key=lambda x:x[0])
			while len(orfs_ordered_lengths[element][s]) > i+1: # orfs_ordered_length is a list that is modified. i gets incremented
				orf = orfs_ordered_lengths[element][s][i][0]
				coord = orfs_coords[element][s][orf]
				j = orfs_ordered_coords[element][s].index(coord) # current largest orf
				k = j+1 # check for overlaps with next in proximity
				# Compare j with successively further away orfs until a non-overlap is reached
				if not k > len(orfs_ordered_coords[element][s])-1:
					while Overlaps( orfs_ordered_coords[element][s][j], orfs_ordered_coords[element][s][k] ):
						coord_removed = orfs_ordered_coords[element][s].pop(k)
						lenkey = coords2lenkey[element][s][coord_removed]
						orfs_ordered_lengths[element][s].remove(lenkey)
						j = orfs_ordered_coords[element][s].index(coord) # current largest orf
						k = j+1 # check for overlaps with next in proximity
						if k > len(orfs_ordered_coords[element][s])-1:
							break
				j = orfs_ordered_coords[element][s].index(coord) # current largest orf
				m = j-1 # check for overlaps with previous in proximity
				# Compare j with successively further away orfs until a non-overlap is reached
				if m > 0:
					while Overlaps( orfs_ordered_coords[element][s][j], orfs_ordered_coords[element][s][m] ):
						coord_removed = orfs_ordered_coords[element][s].pop(m)
						lenkey = coords2lenkey[element][s][coord_removed]
						orfs_ordered_lengths[element][s].remove(lenkey)
						j = orfs_ordered_coords[element][s].index(coord) # current largest orf
						m = j-1 # check for overlaps with previous in proximity
						if m < 0:
							break

				i += 1
			
			orf_ids = [ p[0] for p in orfs_ordered_lengths[element][s] ]
			best_orfs = [ [element, s, orf_id, orfs_seqs_dct[element][s][orf_id], orfs_coords[element][s][orf_id]]  for orf_id in orf_ids ]
			best_orf_sets[s] = best_orfs

		if best_orf_sets['+'] == None and best_orf_sets['-'] == None:
			sys.exit('bestORFs() did not populate best_orf_sets')
		if best_orf_sets['+'] == None:
			best_orfs = best_orf_sets['-']
			strand_used = '-'
		elif best_orf_sets['-'] == None:
			best_orfs = best_orf_sets['+']
			strand_used = '+'
		else:
			lengths_plus = sum([ i[4][1]-i[4][0]+1 for i in best_orf_sets['+']])
			lengths_minus = sum([ i[4][1]-i[4][0]+1 for i in best_orf_sets['-']])
			if lengths_plus > lengths_minus:
				best_orfs = best_orf_sets['+']
				strand_used = '+'
			else:
				best_orfs = best_orf_sets['-']
				strand_used = '-'
			
		
		best_orfs.sort(key=lambda x:x[4][0])
		for orf in best_orfs:
			element, s, orf_id, seq, coords = orf
			start, end = coords
			outgff = '{0}/{1}.orfs.gff'.format(outdir, fasta.split('/')[-1])
			with open(outgff, 'a') as outFl:
				outFl.write('{0}\tgetorf\tORF\t{1}\t{2}\t.\t{3}\t.\tParent={4};translated_seq={5}\n'.format(element, start, end, strand_used, element, seq))

	# ORFs in dict with el:coords



#Chr3	LTRdigest	protein_match	13625394	13625634	9.6e-16	-	.	Parent=LTR_retrotransposon997;reading_frame=1;name=rve
#Chr3	LTRdigest	protein_match	13626881	13627349	6.5e-15	-	.	Parent=LTR_retrotransposon997;reading_frame=2;name=RVT_1
#Chr3	LTRdigest	protein_match	13628385	13628553	1.1e-07	-	.	Parent=LTR_retrotransposon997;reading_frame=1;name=Retrotrans_gag
#Chr3	LTRdigest	protein_match	13645536	13645815	1.4e-25	+	.	Parent=LTR_retrotransposon998;reading_frame=0;name=Retrotrans_gag
#Chr3	LTRdigest	protein_match	13649572	13649680	8.3e-09	+	.	Parent=LTR_retrotransposon998;reading_frame=1;name=ATHILA
#Chr3	LTRdigest	protein_match	13649675	13650884	0	+	.	Parent=LTR_retrotransposon998;reading_frame=2;name=ATHILA


args = sys.argv
outdir = '.'
strandsDct = {}
fasta = args[1]
gff = args[2]
bestORFs(fasta, outdir, gff, minLen=100)

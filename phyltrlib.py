#!/usr/bin/env python3

import re
import sys
import os
import subprocess
from Bio import SeqIO, AlignIO


# Things to implement later: decide how to store multiple non-overlapping domains with the same name in the same LTR-RT

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



class LTR_RT:

	'''

	Holds information about LTR Retrotransposons and links
	features to their corresponding line in GFF3 file.

	Attributes:

	name # LTR_retrotransposon1
	number # 1
	domain_count # {domain1:1, domain2:2, ... }
	domain_architecture # { domain1_1:gff_line_number, domain2_1:gff_line_number, ... }
	ltrs # { left:gff_line_number, right:gff_line_number }
	internal_distance # 3000

	'''

	def __init__(self, GFF3_line_obj): # meant to be initialized when a GFF3 line with type long_terminal_repeat

		assert GFF3_line_obj.type == 'long_terminal_repeat', 'GFF3_line obj being used to initiate LTR_RT obj not for feature type long_terminal_repeat'

		assert not GFF3_line_obj.line_number == None, 'GFF3 line number needs to be passed to the initialization of GFF3_line obj for {0} using the kwarg line_number= '.format(GFF3_line_obj.name)

		self.name = GFF3_line_obj.attributes['Parent']
		self.number = int(self.name.lstrip('LTR_retrotransposon')) # in LTRDigest output, the parent LTR RT feature has IDs like LTR_retrotransposon1, LTR_retrotransposon2, ...
		self.domain_count = {}
		self.domain_architecture = []
		self.domain_architecture_coords = []
		self.domains = {}
		self.ltrs = {}
		self.ltr_coords = {}
		self.internal_distance = None

	def add_ltr(self, GFF3_line_obj):

		'''

		Requires GFF3_line object has been given the line number to the kwargs line_number.

		'''

		assert not len(self.ltrs) >= 2, 'More than two long terminal repeats found for this LTR-RT, see {0}, line {1} in input GFF3'.format(self.name, GFF3_line_obj.line_number)

		ltr = 'R' if 'L' in self.ltrs else 'L' 

		self.ltrs[ltr] = GFF3_line_obj.line_number
		self.ltr_coords[ltr] = (GFF3_line_obj.start, GFF3_line_obj.end)

		if ltr == 'R':
			self.internal_distance = int(self.ltr_coords['R'][0]) + 1 - int(self.ltr_coords['L'][1])

	def add_domain(self, GFF3_line_obj):


		domain = GFF3_line_obj.attributes['Name']
		domain_coords = (GFF3_line_obj.start, GFF3_line_obj.end)
		self.domain_count[domain] = self.domain_count[domain] + 1 if domain in self.domain_count else 1

		if domain in self.domains:

			if GFF3_line_obj.score < self.domains[domain]['score']:
				self.domains[domain] = { 'name':domain, 'domain_number':self.domain_count[domain], 'line_number':GFF3_line_obj.line_number, 'score':GFF3_line_obj.score }

		else:

			self.domains[domain] = { 'name':domain, 'domain_number':self.domain_count[domain], 'line_number':GFF3_line_obj.line_number, 'score':GFF3_line_obj.score }


		if self.domain_architecture == [] and self.domain_architecture_coords == []:
			self.domain_architecture.append(domain)
			self.domain_architecture_coords.append(domain_coords)

		elif domain == self.domain_architecture[-1] and int(domain_coords[0]) < int(self.domain_architecture_coords[-1][1]): # don't add domain to this element's domain arch if it overlaps the last one
			pass

		else:

			self.domain_architecture.append(domain)
			self.domain_architecture_coords.append(domain_coords)


def rename_fasta_seq_headers(in_flpath, in_header_pattern, header_to_name_map, out_flpath):

	'''
	in_fasta_filepath - the fasta file with headers to be renamed.
	in_header_pattern - a regex that will match the part of the header that corresponds to keys in header_to_name_map
	header_to_name_map - a dict with keys as part of in_fasta_filepath seq headers that match in_header_pattern and values as the abbreviated LTR-RT ID (e.g. 24_1, which corresponds to LTR_retrotransposon24)

	This function depends on the re and Bio.SeqIO modules
	'''

	
	fasta_file = list(SeqIO.parse(in_flpath, 'fasta'))

	in_header_pattern = re.compile(in_header_pattern)

	for seq in fasta_file:

		seq_id = re.search(in_header_pattern, seq.id).group(1)
		seq.id = header_to_name_map[seq_id]
		seq.description = ''
	
	SeqIO.write(fasta_file, out_flpath, 'fasta')





def write_ltrs_gff3(data):

	'''
	Writes items from a list data to a file (intended to be a GFF3 file in ltran.py)
	For use with -ltrs option in ltran.py
	'''

	if len(data) < 3:
		print('Less than 2 LTRs...?  ...GFF3 not written. DATA:\n{0}'.format('\n'.join(data)), file=sys.stderr)

	elif len(data) > 3:
		print('More than 2 LTRs...?  ...GFF3 not written. DATA:\n{0}'.format('\n'.join(data)), file=sys.stderr)

	else:
		with open(data[-1], 'w') as out_fl: # data[-1] should have the output filepath while data[0:2] are the GFF3 lines
			out_fl.write(''.join(data[0:2]))




def call_process(call_str):

	'''
	Just runs subprocess.call. For simplifying parallel execution of ltran.py
	using the multiprocessing module.
	'''

	subprocess.call(call_str, shell=True)





def count_end_gaps(aln):

	'''
	- Input is a two-element list of seqs of identical length (a pairwise alignment).

	- Returns the number of bases to subtract from the alignment length. The reason
		being that since this is meant for alignments between two LTRs, LTR boundaries
		may be off and part of LTR may have been left in the genome. If LTR boundaries
		were accurate it's possible there would be no end gaps. The length is used in
		ltran.py for estimating sequence divergence.

	- (The number of gaps from the left end in the seq with the most gaps on the 
		left end plus the number from the right end)

	Algorithm: walk inward from each seq end one char at a time and count how many gaps.
		Stop a walk when the char is not the gap char. When all walks are non-gap chars,
		return the sum of the highest left and highest right count.
	'''

	if len(aln[0]) != len(aln[1]): # make sure input seqs are the same length

		raise ValueError('Sequences in alignment have different length') 


	# Keys are all same in these dicts. The number from the key can be stripped and used as the index with which to access the sequence.
	end_gaps = {'L0':0, 'L1':0, 'R0':0, 'R1':0}			# keeps track of how many gaps on each end of each sequence
	chars = {'L0':'-', 'L1':'-', 'R0':'-', 'R1':'-'}		# keeps track of the current character in each walk
	pos = {'L0':0, 'L1':0, 'R0':-1, 'R1':-1}			# keeps track of the current position of each walk
	#non_gap = {'L0':False, 'L1':False, 'R0':False, 'R1':False}	# keeps track of if no gap. DON"T THINK I NEED THIS

	while chars['L0'] == '-' or chars['L1'] == '-' or chars['R0'] == '-' or chars['R1'] == '-':

		for k in end_gaps:

			#if non_gap[k] == False:
			if chars[k] == '-':

				seq = int(k[1]) # Strip num from key k to get index of sequence

				chars[k] = aln[seq][pos[k]] # get next position in current walk

				#if chars[k] != '-':
				#	non_gap[k] = True
				#else:
				#	end_gaps[k] += 1

				if chars[k] == '-': # NEW IMPROVED, YET TO TEST # add one to count if this char is a gap char

					end_gaps[k] += 1 # NEW IMPROVED, YET TO TEST

			if 'L' in k: # get next index to continue walk from left side

				pos[k] += 1

			if 'R' in k: # get next indext to continue walk from right side

				pos[k] -= 1

	return ( max([end_gaps['L0'],end_gaps['L1']]) + max([end_gaps['R0'],end_gaps['R1']]) ) # return gap count







def fastas2supermatrix(**kwargs):

	'''
	needs biopython installed

	kwargs expected:
	input_dir='path to dir containing fasta files to be concatenated'
	input_dir is expected to have only fasta files to be concatenated
	1. Check that lengths of all seqs in each input fasta alignment are the same, report seq lengths and quit if they're different.
	2. For each fasta file input, add each sequence identifier as a key in a dictionary with the value of the key a dictionary with the name of the fasta file the key and the value the sequence.
	3. Make a set of unique sequence ids.
	4. Build a new dicitonary with keys as sequence identifiers and values are concatenated sequences, with Ns or gaps if no sequence is available for that identifier for a given fasta file input.
	5. Write new fasta of concatenated sequences.
	'''


	input_fastas = [ fl for fl in os.listdir(kwargs['input_dir']) if not fl.startswith('.') ] # get names of input files. dir with input files passed to this function as kwarg input_dir.

	seqs_dct = {} # will become a nested dict holding sequences from each fasta file for each sequence header

	seq_lengths = {} # will hold the length of each alignment


	# process each fasta seq
	for fasta in input_fastas:

		fasta_contents = list(SeqIO.parse('{0}/{1}'.format(kwargs['input_dir'], fasta), 'fasta')) # read in alignment

		# make sure seqs in the alignment are the same length. If they aren't then print seq ids and lengths to stderr and quit
		all_seqs_lengths = { seq.id:len(seq) for seq in fasta_contents }
		if len(set(all_seqs_lengths.values())) != 1:
			print('Sequences in file: {0} not all of equal length'.format(fasta), file=sys.stderr)
			print('\n'.join( ['sequence\tlength'] + [ '{0}\t{1}'.format(seq, all_seqs_lengths[seq]) for seq in all_seqs_lengths ] ), file=sys.stderr)
			sys.exit()

		seq_lengths[fasta] = len(fasta_contents[0])


		
		# Add sequences to seqs_dct
		for seq in fasta_contents:

			if seq.id in seqs_dct:

				# check for duplicate seq id and stop program if any duplicates are found
				if fasta in seqs_dct[seq.id]:
					print('Duplicate sequence identifier in {0}: {1}'.format(fasta, seq.id), file=sys.stderr)
					sys.exit()

				seqs_dct[seq.id][fasta] = str(seq.seq) # if no duplicate add sequence to nested dict with key as fasta file name

			else:
				seqs_dct[seq.id] = { fasta:str(seq.seq) } # add sequence to nested dict with key as fasta file name 

	# print concatenated fasta supermatrix

	if 'output_fl' in kwargs:

		out_fl = open(kwargs['output_fl'], 'a')

	for seq in seqs_dct:

		output_seq = ''

		for fasta in input_fastas:

			if fasta not in seqs_dct[seq]:

				output_seq += 'N' * seq_lengths[fasta]

			else:

				output_seq += seqs_dct[seq][fasta]
		
		if 'output_fl' in kwargs:

			out_fl.write('>{0}\n'.format(seq))

			out_fl.write('{0}\n'.format(output_seq))

		else:

			print('>{0}'.format(seq))

			print(output_seq)

	if 'output_fl' in kwargs:

		out_fl.close()



def bedtoolsid2attr(gff_flpath, attr='ID', strand=False, lstrip=None):

	'''
	bedtoolsid2attr.py -gff <gff> [-attr <str>] [-strand]

	'''

	STRAND = strand
	attr_pat = re.compile('{0}=(.+?)(;|$)'.format(attr))
	map_dct = {}

	with open(gff_flpath) as in_fl:
		for line in in_fl:
			if not line.startswith('#'):
				fields = line.strip().split('\t')
				seqid = fields[0]
				start = str(int(fields[3]) - 1)
				end = fields[4]
				strand = fields[6]
				attr_value = re.search(attr_pat, fields[8]).group(1)
				
				if not lstrip == None:
					attr_value = attr_value.lstrip(lstrip)

				if STRAND:
					map_dct['{0}:{1}-{2}({3})'.format(seqid, start, end, strand)] = attr_value

				else:
					map_dct['{0}:{1}-{2}'.format(seqid, start, end)] = attr_value

	return map_dct

def rename_fasta_seq_headers(in_flpath, in_header_pattern, header_to_name_map, out_flpath):

	'''
	in_fasta_filepath - the fasta file with headers to be renamed.
	in_header_pattern - a regex that will match the part of the header that corresponds to keys in header_to_name_map
	header_to_name_map - a dict with keys as part of in_fasta_filepath seq headers that match in_header_pattern and values as the abbreviated LTR-RT ID (e.g. 24_1, which corresponds to LTR_retrotransposon24)

	This function depends on the re and Bio.SeqIO modules
	'''

	
	fasta_file = list(SeqIO.parse(in_flpath, 'fasta'))

	in_header_pattern = re.compile(in_header_pattern)

	for seq in fasta_file:

		seq_id = re.search(in_header_pattern, seq.id).group(1)
		seq.id = header_to_name_map[seq_id]
		seq.description = ''
	
	SeqIO.write(fasta_file, out_flpath, 'fasta')

def ChangeFastaHeaders(inputFastaPath, inputGFFpath, attribute='ID'):
	'''
	Creates map of bedtools getfasta-style features, reads in inputFasta, writes newFasta, deletes inpuFasta, renames newFasta as inputFasta
	'''

	bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
	newFasta = '{0}.new.tmp'.format(inputFastaPath)
	header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
	rename_fasta_seq_headers(inputFastaPath, header_pattern, bedtoolsIDmap, newFasta)
	os.replace(newFasta, inputFastaPath)

def ChangeFastaHeadersMultiprocessing(bundle):
	'''
	Creates map of bedtools getfasta-style features, reads in inputFasta, writes newFasta, deletes inpuFasta, renames newFasta as inputFasta
	'''
	inputFastaPath = bundle[0]
	inputGFFpath = bundle[1]
	attribute=bundle[2]

	bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
	newFasta = '{0}.new.tmp'.format(inputFastaPath)
	header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
	rename_fasta_seq_headers(inputFastaPath, header_pattern, bedtoolsIDmap, newFasta)
	os.replace(newFasta, inputFastaPath)

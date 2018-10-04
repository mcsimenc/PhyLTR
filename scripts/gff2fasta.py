#!/usr/bin/env python3
import re
import os
import sys
import subprocess
from Bio import SeqIO

def help():
	print('''
 Usage:
 ------------
 gff2fasta.py -gff <file> -fasta <file> [-attr <attr>]
 
 Description:
 ------------
 Requires BEDtools. Extracts seqs from -fasta for features in -gff
 using the attribute from -gff given by -attr (default: ID)
 
 Options:
 ------------
 -attr <attr>	Attribute from which to derive fasta seq name
		''')
	sys.exit(0)

def bedtoolsid2attr(gff_flpath, attr='ID', strand=False, lstrip=None):
	'''
	This takes a GFF and creates a map of the current attribute specified by
	attr, to the expected names that will be output by bedtools getfasta (a tool
	for extracting sequences from a FASTA which correspond to features in a GFF3.

	Also as a standalone script:

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
	I/O for changing fasta sequence headers. Uses regex matching and dictionary associating current with new name.
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

def makecall(call, stdout=None, stderr=None, stdin=None):
	'''
	Handles running subprocess.call. Used when making calls without multiprocessing
	'''
	if stdout == None and stderr == None and stdin == None:
		subprocess.call(call)
	elif stdout != None:
		with open(stdout, 'w') as outfl:
			if stderr != None:
				with open(stderr, 'w') as errfl:
					if stdin == None:
						subprocess.call(call, stdout=outfl, stderr=errfl)
					else:
						with open(stdin, 'r') as inFl:
							subprocess.call(call, stdin=inFl, stdout=outfl, stderr=errfl)
			elif stderr == None:
				if stdin == None:
					subprocess.call(call, stdout=outfl)
				else:
					with open(stdin, 'r') as inFl:
						subprocess.call(call, stdin=inFl, stdout=outfl)

	elif stderr != None and stdout == None:
		with open(stderr, 'w') as errfl:
			if stdin == None:
				subprocess.call(call, stderr=errfl)
			else:
				with open(stdin, 'r') as inFl:
					subprocess.call(call, stdin=inFl,  stderr=errfl)
		
	elif stdin != None and stderr == None and stdout == None:
		with open(stdin, 'r') as inFl:
			subprocess.call(call, stdin=inFl)

def ChangeFastaHeaders(inputFastaPath, inputGFFpath, attribute='ID'):
	'''
	Creates map of bedtools getfasta-style features, reads in inputFasta, writes newFasta, deletes inpuFasta, renames newFasta as inputFasta
	'''
	bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
	newFasta = '{0}.new.tmp'.format(inputFastaPath)
	header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
	rename_fasta_seq_headers(inputFastaPath, header_pattern, bedtoolsIDmap, newFasta)
	os.replace(newFasta, inputFastaPath)

def removeRedundant(fastaPth):
	seqnames = set()
	tmp = '{0}.nonredundant'.format(fastaPth)
	if os.path.isfile(tmp):
		os.remove(tmp)
	with open(tmp, 'w') as outFl:
		with open(fastaPth, 'r') as inFl:
			SKIP = False
			for line in inFl:
				if line.startswith('>'):
					if line in seqnames:
						SKIP = True
						continue
					else:
						SKIP = False
						seqnames.add(line)
						outFl.write(line)
				else:
					if SKIP:
						continue
					else:
						outFl.write(line)
	os.rename(tmp, fastaPth)

def getfasta(inGFFpth, fastaRefPth, outFastaPth, headerKey='ID'):
	'''
	generic function to run bedtools getfasta
	'''
	call = [  'bedtools', 'getfasta', '-fi', fastaRefPth, '-s', '-bed', inGFFpth ]
	makecall(call, stdout=outFastaPth)
	ChangeFastaHeaders(outFastaPth, inGFFpth, attribute=headerKey)
	removeRedundant(outFastaPth)

if __name__ == '__main__':
	# Parse args
	args = sys.argv
	if "-help" in args or "-h" in args or len(args) < 5:
		help()
	if '-attr' in args:
		attr = args[args.index('-attr')+1]
	else:
		attr = 'ID'
	fastaPth = args[args.index('-fasta')+1]
	gffPth = args[args.index('-gff')+1]
	# Get sequences
	getfasta(gffPth, fastaPth, '{0}.fasta'.format(gffPth), headerKey=attr)

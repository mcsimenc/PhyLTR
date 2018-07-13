#!/usr/bin/env python3

import re
import os
import sys
from Bio import SeqIO
import subprocess

def help():
	print('''
 Usage:
 ------------
 domainSearch.py -gff <gff3> -ref <fasta> -prot <fasta> -procs <int>
 
 Description:
 ------------
 This script was written for comparing putative protein-coding domains
 and ORFs in a GFF3 as output from PhyLTR (LTRharvest+LTRdigest+AnnotateORFs).
 The required programs are bedtools getfasta, and makeblastdb and blastp from
 BLAST+. Their paths should be in a text file named CONFIG file located in the
 same directory as this script using the same format used for PhyLTR.

 Process:
 -----------
 1. Putative domain sequences encoded in -gff are extracted from -ref and
    translated
 2. Translated sequences are compared to protein sequences in -prot, which could
    be, for example, the protein sequences for a
 
 Mandatory flags:
 -----------
 -gff	 <path>	GFF3 file with LTRdigest-format LTR retrotransposon features
 
 -ref    <path>	Nucleotide FASTA file that is the reference associated with -gff
 
 -prot	 <path>	Protein FASTA file to convert to a blast database and

 -procs	 <int>	Number of processors to use for blastp

 -out	 <path> Output file path - you may need to give full path

 Optional flags:
 -----------
 -evalue <int>	max E-value for blastp. Default 1e-2
 
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


def removeRedundant(fastaPth):
	seqnames = set()
	tmp = '{0}.nonredundant'.format(fastaPth)
	if os.path.isfile(tmp):
		os.remove(tmp)
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
					with open(tmp, 'a') as outFl:
						outFl.write(line)
			else:
				if SKIP:
					continue
				else:
					with open(tmp, 'a') as outFl:
						outFl.write(line)
	os.rename(tmp, fastaPth)


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
	header_pattern='(.+?:\d+?-\d+?)\D'
	rename_fasta_seq_headers(inputFastaPath, header_pattern, bedtoolsIDmap, newFasta)
	os.replace(newFasta, inputFastaPath)


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


def getfasta(inGFFpth, fastaRefPth, outFastaPth, headerKey='ID'):
	'''
	generic function to run bedtools getfasta
	'''
	call = [  executables['bedtools'], 'getfasta', '-fi', fastaRefPth, '-s', '-bed', inGFFpth ]
	call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2}'.format(fastaRefPth, inGFFpth, outFastaPth)
	makecall(call, stdout=outFastaPth)
	ChangeFastaHeaders(outFastaPth, inGFFpth, attribute=headerKey)
	removeRedundant(outFastaPth)


def runblast(query, subject, out, evalue, outfmt, percid=None, blast='blastn', procs=1):
	'''
	runs blast programs
	if subject does not have a blastdb one will be created
	expects procs to be a global variable
	'''
	dbtypes = { 'blastn':'nucl', 'blastp':'prot', 'blastx':'prot', 'tblastx':'nucl', 'tblastn':'nucl' }

	assert blast in dbtypes,'runblast() blast param must be one of blastn, blastp, blastx, tblastn, tblastx'

	makeblastdb_call_string = '{0}/makeblastdb -in {1} -dbtype {2} -parse_seqids'.format(executables['blast'], subject, dbtypes[blast] )
	makeblastdb_call = [ '{0}/makeblastdb'.format(executables['blast']), '-in', subject, '-dbtype', dbtypes[blast], '-parse_seqids' ]
	makecall(makeblastdb_call)
	if percid == None or blast != 'blastn': # -perc_identity is only a blastn option
		blast_call = [ '{0}/{1}'.format(executables['blast'], blast), '-db', subject, '-query', query, '-evalue', str(evalue), '-outfmt', str(outfmt), '-num_threads', str(procs)]
		blast_call_string = '{0}/{1} -db {2} -query {3} -evalue {4} -outfmt {5} -num_threads {6} 1>{7} 2>{8}'.format(executables['blast'], blast, subject, query, evalue, outfmt, procs, out, '{0}/runblast.{1}.err'.format('/'.join(out.split('/')[:-1]), blast))
	else:
		blast_call = [ '{0}/{1}'.format(executables['blast'], blast), '-db', subject, '-query', query, '-evalue', str(evalue), '-outfmt', str(outfmt), '-num_threads', str(procs), '-perc_identity', str(percid) ]
		blast_call_string = '{0}/{1} -db {2} -query {3} -evalue {4} -outfmt {5} -num_threads {6} -perc_identity {7} 1>{8} 2>{9}'.format(executables['blast'], blast, subject, query, evalue, outfmt, procs, percid, out, '{0}/runblast.{1}.err'.format('/'.join(out.split('/')[:-1]), blast))
	makecall(blast_call, stdout=out, stderr='{0}/runblast.{1}.err'.format('/'.join(out.split('/')[:-1]), blast))


def addID2LTRdigestDomains(inGFF):
	outPth = inGFF + '.new'
	with open(inGFF) as inFl:
		with open(outPth, 'w') as outFl:
			domainTracker = {}
			for line in inFl:
				if '\tprotein_match\t' in line:
					gffLine = GFF3_line(line)
					el = gffLine.attributes['Parent']
					domain = gffLine.attributes['Name']
					if el in domainTracker:
						if domain in domainTracker[el]:
							ID = '{0}.{1}.{2:02d}'.format(el, domain, domainTracker[el][domain])
							domainTracker[el][domain] += 1
						else:
							domainTracker[el][domain] = 1
							ID = '{0}.{1}.{2:02d}'.format(el, domain, domainTracker[el][domain])
					else:
						domainTracker[el] = {domain:1}
						ID = '{0}.{1}.{2:02d}'.format(el, domain, domainTracker[el][domain])
					gffLine.attributes['ID'] = ID
					gffLine.attributes_order.insert(0,'ID')
					gffLine.refreshAttrStr()
					outFl.write(str(gffLine)+'\n')
				else:
					outFl.write(line)
	os.rename(inGFF, inGFF+'.OLD')
	os.rename(outPth, inGFF)

def writeDomainsGFF(inGFF):
	'''
	Outfile is called <inFl>.domains and includes features from inGFF of type
	ORF and protein_match
	'''
	with open(inGFF) as inFl:
		with open(inGFF+'.domains', 'w') as outFl:
			for line in inFl:
				if '\tprotein_match\t' in line or '\tORF\t' in line:
					outFl.write(line)

def translate(inFASTA):
	getorf_call = [ executables['getorf'], '-sequence', inFASTA, '-outseq', 'tmp.getorf.out' ]
	makecall(getorf_call)
	lastDomain = None
	longestLength = 0
	longestSeq = ''
	thisDomain = None
	thisLength = 0
	WASLASTLINESEQ = False
	fasta_file = list(SeqIO.parse('tmp.getorf.out', 'fasta'))
	outSeqs = {}
	with open(inFASTA+'.prot', 'w') as outFl:
		for S in fasta_file:
			domain = '_'.join(S.id.split('_')[:-1])
			if domain in outSeqs:
				if len(S.seq) > len(outSeqs[domain]):
					outSeqs[domain] = S.seq
			else:
				outSeqs[domain] = S.seq
		for domain in sorted(list(outSeqs.keys())):
			outFl.write('>{0}\n{1}\n'.format(domain, outSeqs[domain]))
						
							
		

args = sys.argv

if not '-gff' in args or not '-ref' in args or not '-out' in args or not '-prot' in args or not '-procs' in args:
	help()
	print('\n MISSING ONE OR MORE OPTIONS\n', file=sys.stderr)
	sys.exit()


# get executable paths from CONFIG file, which should be in the same directory as this script
executables = {}
commentPattern = re.compile('#.*$')
with open('{0}/CONFIG'.format(os.path.dirname(os.path.realpath(__file__)))) as config_file:
	paths = [ re.sub(commentPattern, '', line) for line in config_file.read().strip().split('\n') ]
	for path in paths:
		if not path == '':
			p = path.split('=')
			executables[p[0]] = p[1]

paths = {}
paths['ref'] = args[args.index('-ref') + 1]
paths['gff3'] = args[args.index('-gff') + 1]
paths['prot'] = args[args.index('-prot') + 1]
paths['outfile'] = args[args.index('-out') + 1]
procs = args[args.index('-procs') + 1]
if '-evalue' in args:
	evalue = args[args.index('-evalue') + 1]
else:
	evalue = '1e-2'


addID2LTRdigestDomains(paths['gff3']) # Add unique IDs for protein domains (LTRdigest doesn't do this)
writeDomainsGFF(paths['gff3'])
getfasta(paths['gff3']+'.domains', paths['ref'], paths['gff3']+'.domains.fasta')
translate(paths['gff3']+'.domains.fasta')
runblast(paths['gff3']+'.domains.fasta.prot', paths['prot'], paths['outfile'], evalue, '7', blast='blastp', procs=procs)

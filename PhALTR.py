#!/usr/bin/env python3

import sys
import re
import os
import time
import os.path
import subprocess
import random
from shutil import copyfile, copytree
from datetime import datetime
from math import ceil
from Bio import SeqIO, AlignIO
from multiprocessing import Pool, Manager
from phaltrlib import *
from copy import copy

def mergeCoords(A,B):
	'''
	takes two tuples and outputs two tuples, which will be identical if the original overlap otherwise will be the originals

	let A = (a1, a2), B = (b1, b2) | a1<=b1, a1<=a2, b1<=b2

	case 1: a2<=b1 ---> output originals

	case 2: b1<a2 && b2>a2 ---> output (a1, b2)

	case 3: b2<=a2 ---> output A

	used in ltr_divergence()
	'''
	assert min(A) <= min(B), "tuples given to mergeCoords in wrong order: A={0}, B={1}".format(A,B)

	if min(B) >= max(A):
		return ((A,B), 0)
	elif min(B) < max(A) and max(B) > max(A):
		output = (min(A),max(B))
		return ((output, output), 1)
	elif max(B) <= max(A):
		return ((A,A), 2)
	else:
		raise Exception("Unexpected result from mergeCoords(A,B) using A={0}, B={1}".format(A,B))

def append2logfile(directory, logfilename, content):
	'''
	Appends string 'content' to directory/logfilename
	'''
	with open('{0}/{1}'.format(directory, logfilename), 'a') as logfile:
		logtime = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
		logfile.write('{0}\n'.format(logtime))
		logfile.write('{0}\n\n'.format(content))


def write2summary(text):
	'''
	Write text to a hardcoded summary file at output_dir/summary
	'''
	append2logfile(paths['output_top_dir'], mainlogfile, 'Writing to summary file at {0}'.format('{0}/summary'.format(paths['output_top_dir'])))
	with open('{0}/summary'.format(paths['output_top_dir']), 'a') as summaryFl:
		summaryFl.write('{0}\n'.format(text))


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


def makecallMultiprocessing(callBundle):
	'''
	Handles running subprocess.call when using multiprocessing.
	Used to iterate over calls with Pool().map
	'''
	call = callBundle[0]
	stdout = callBundle[1]
	stderr = callBundle[2]
	stdin = callBundle[3]
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


def MakeDir(pathsname, path):
	'''
	Makes a directory path and stores it in the global dict paths under the key pathsname and writes to logfile.
	'''
	global paths

	paths[pathsname] = path
	if not os.path.exists(path): 
		os.makedirs(path)
		append2logfile(paths['output_top_dir'], mainlogfile, 'Created dir for {0}:\n{1}'.format(pathsname, path))


def addStrandToGFF(strandDct, GFFpth):
	'''
	Updates strand field for element with ? as strand based on Dfam and Repbase results.
	Provide a dictionary with LTR RT # (e.g. 4 for LTR_retrotransposon4) as keys and
	strand as values. A new GFF will be written, and the old one removed.
	'''
	newgff = '{0}.updatingstrandinprocess'.format(GFFpth)

	with open(GFFpth) as inFl:
		for line in inFl:
			if not line.startswith('#'):
				gffLine = GFF3_line(line)
				if gffLine.strand == '?':
					if 'ID' in gffLine.attributes:
						elNum = re.search('\d*$', gffLine.attributes['ID']).group(0)
					elif 'Parent' in gffLine.attributes:
						elNum = re.search('\d*$', gffLine.attributes['Parent']).group(0)
				
					else:
						print('WARNING:\taddStrandToGFF() found the following GFF file to contain the following line where strand is unknown and the attributes lack ID or Parent keys\n{0}\n{1}'.format(GFFpth, line.strip()), file=sys.stderr)

					if elNum in strandDct:
						gffLine.strand = strandDct[elNum]
						with open(newgff) as outGFFfl:
							outGFFfl.write(str(gffLine))
					else:
						with open(newgff) as outGFFfl:
							outGFFfl.write(line)
			else:
				with open(newgff) as outGFFfl:
					outGFFfl.write(line)
	os.rename(newgff, GFFpth)


def RemoveNonLTRretrotransposons(LTRdigestGFFfl, annotAttr2DbDict, outputFlName, REPORTCONFLICTS=True, KEEPCONFLICTS=False, KEEPNOCLASSIFICATION=False, logFilePth='conflictingAnnotations.log'):
	'''
	Removes non-LTR retrotransposons from LTRdigest GFF3 input based on 
	classifications from e.g. Repbase or Dfam

	Parses LTRdigestGFFfl and checks attributes specified in annotAttr2DbDict
	for annotations and looks for them in the Db specified in annotAttr2DbDict.
	If found, that element is output. If conflicting annotations are found,
	REPORTCONFLICTS=True writes conflicts to logfile at logFilePth; if KEEPCONFLICTS=True,
	the element is retained in the output GFF.

	Structure of annotAttr2DbDict should be:

		{ attr1:DbFlName1, ... }

	where attr are keys in the attributes field of the input GFF3 that contain
	classification information derived from the database of which DbFlName should
	be the subset of classifications in the database that are LTR retrotransposons.

	If an attribute's value is "None", it is treated as "no information". If all
	attributes specified in annotAttr2DbDict are "None" and KEEPNOCLASSIFICATION=True,
	that element is output.
	'''
	LTR_retrotransposon_GFF_lines = {}
	FOUNDNONLTR = False
	Db_LTR_retrotransposon_features = {}
	Db_files = list(annotAttr2DbDict.values())

	for dbFlPth in Db_files:
		with open(dbFlPth) as dbFl:
			for line in dbFl:
				flName = dbFlPth.split('/')[-1]
				if flName in Db_LTR_retrotransposon_features:
					Db_LTR_retrotransposon_features[flName].add(line.strip())
				else:
					Db_LTR_retrotransposon_features[flName] = set([line.strip()])

	logfile = open(logFilePth, 'a')
	with open(LTRdigestGFFfl) as gffFl:
		for line in gffFl:
			if line.startswith('#'):
				continue
			gffLine = GFF3_line(line)
			if gffLine.type == 'repeat_region': # First line in LTR retrotransposon block in LTRdigest GFF3
				FOUNDNONLTR = False
				if gffLine.type in LTR_retrotransposon_GFF_lines:
					sys.exit('Line\n{0}\nout of order'.format(line.strip()))
				LTR_retrotransposon_GFF_lines[gffLine.attributes['ID']] = [ line.strip() ]
				continue
			elif gffLine.type == 'LTR_retrotransposon':
				## Check if LTR
				NoClassification = []
				LTRmatching = {}
				for attr in annotAttr2DbDict:
					db = annotAttr2DbDict[attr].split('/')[-1]
					annot = gffLine.attributes[attr]
					if not annot == 'None':
						if annot in Db_LTR_retrotransposon_features[db]:
							LTRmatching[attr] = True
						else:
							LTRmatching[attr] = False
					else:
						NoClassification.append('None')

				LTRmatching_values = list(LTRmatching.values())
				if True in LTRmatching_values and False in LTRmatching_values:

					# Found conflicting annotations: LTR and Non-LTR
					if REPORTCONFLICTS:
						logfile.write('{0}\thas conflicting annotations LTR and Non-LTR\t{1}\n'.format(gffLine.attributes['Parent'], '\t'.join([ '{0}={1}'.format(attr, gffLine.attributes[attr]) for attr in LTRmatching ])))
					if KEEPCONFLICTS:
						LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']].append(line.strip())
					else:
						FOUNDNONLTR = True
						del(LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']])
						continue

				elif 'None' in NoClassification and LTRmatching_values == []:
						if KEEPNOCLASSIFICATION:
							LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']].append(line.strip())
							continue

						FOUNDNONLTR = True
						del(LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']])
						continue

				elif NoClassification == [] and LTRmatching_values == []:
					logfile.write('{0}\thas no classification\n'.format(gffLine.attributes['Parent']))
					FOUNDNONLTR = True
					del(LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']])
					continue
				elif False in LTRmatching_values:
					FOUNDNONLTR = True
					del(LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']])
					continue
				else:
					LTR_retrotransposon_GFF_lines[gffLine.attributes['Parent']].append(line.strip())
					continue
			else:
				if FOUNDNONLTR:
					continue

				parent = gffLine.attributes['Parent']

				if 'LTR_retrotransposon' in parent:
					parent = 'repeat_region{0}'.format(parent.lstrip('LTR_retrotransposon'))
				if not parent in LTR_retrotransposon_GFF_lines:
					sys.exit('Line\n{0}\nParent attribute not in LTR retrotransposon dictionary'.format(line.strip()))

				LTR_retrotransposon_GFF_lines[parent].append(line.strip())

	with open(outputFlName, 'w') as outFl:
		for element in sorted(list(LTR_retrotransposon_GFF_lines.keys())):
			outFl.write('\n'.join(LTR_retrotransposon_GFF_lines[element]))
			outFl.write('\n###\n')

	logfile.close()


def writeLTRretrotransposonInternalRegions(inputGFFpth, outputGFFpth, elementSet=None, truncateParent=False):
	''' 
	Requires Class GFF3_line
	Writes GFF3 for region between two LTRs from a LTRharvest-type file
	Only for elements in elementSet if provided, if elementSet == None, all elements are written.
	if truncateParent=True, Parent attribute has 'LTR_retrotranspson' trimmed from it
	'''
	with open(inputGFFpth, 'r') as inGFF:
		currentNewElement = GFF3_line()
		for line in inGFF:
			if '\tlong_terminal_repeat\t' in line:
				gffLine = GFF3_line(line)
				if elementSet == None or (elementSet != None and gffLine.attributes['Parent'] in elementSet):
						if not currentNewElement.start == None:
							if truncateParent == True:
								gffLine.attributes['Parent'] = gffLine.attributes['Parent'].lstrip('LTR_retrotransposon')
							currentNewElement.end = gffLine.start - 1

							if currentNewElement.end - currentNewElement.start + 1 <= 0:
								currentNewElement = GFF3_line()
								continue

							assert currentNewElement.attributes['Parent'] == gffLine.attributes['Parent'], 'GFF long_terminal_repeats may be out of order. Check near {0} or {1} in {2}'.format(currentNewElement.attributes['Parent'], gffLine.attributes['Parent'], inputGFFpth)

							currentNewElement.seqid = gffLine.seqid
							currentNewElement.source = 'PhALTR'
							currentNewElement.type = 'LTR_retrotransposon_InternalRegion'
							currentNewElement.score = '.'
							currentNewElement.strand = gffLine.strand
							currentNewElement.phase = '.'

							with open(outputGFFpth, 'a') as outGFFfl:
								outGFFfl.write('{0}\n'.format(str(currentNewElement)))
								currentNewElement = GFF3_line()
						else:
							currentNewElement.start = gffLine.end + 1
							if truncateParent == True:
								currentNewElement.attributes['Parent'] = gffLine.attributes['Parent'].lstrip('LTR_retrotransposon')
							else:
								currentNewElement.attributes['Parent'] = gffLine.attributes['Parent']
							currentNewElement.attributes_order.append('Parent')
							currentNewElement.refreshAttrStr()
								

def writeLTRretrotransposonGFF(inputGFFpth, outputGFFpth, elementSet=None):
	''' 
	Requires Class GFF3_line
	Writes GFF3 for LTR_retrotransposon type features from a LTRharvest-type file.
	Only for elements in elementSet if provided, if elementSet == None, all elements are written.
	If truncateParent=True, Parent attribute has 'LTR_retrotranspson' trimmed from it.
	'''
	global paths

	if os.path.isfile(outputGFFpth):
		os.remove(outputGFFpth)

	append2logfile(paths['output_top_dir'], mainlogfile, 'Writing LTR_retrotransposon features:\n{0}\nfrom:\n{1}\nto:\n{2}'.format(','.join(sorted(list(elementSet))),inputGFFpth, outputGFFpth))
	with open(inputGFFpth, 'r') as inGFF:
		currentNewElement = GFF3_line()
		for line in inGFF:
			if '\tLTR_retrotransposon\t' in line:
				gffLine = GFF3_line(line)
				if elementSet == None or (elementSet != None and gffLine.attributes['ID'] in elementSet):
					with open(outputGFFpth, 'a') as outFl:
						gffLine.attributes['ID'] = gffLine.attributes['ID'][19:]
						gffLine.refreshAttrStr()
						outFl.write(str(gffLine)+'\n')


def writeLTRsGFF(inputGFFpth, outputGFFpth, elementSet=None):
	''' 
	Requires Class GFF3_line
	Writes one GFF3 for each pair of LTRs from a LTRharvest-type file.
	Only for elements in elementSet if provided, if elementSet == None, all elements are written.
	if truncateParent=True, Parent attribute has 'LTR_retrotranspson' trimmed from it
	'''
	global paths

	if os.path.isfile(outputGFFpth):
		os.remove(outputGFFpth)

	append2logfile(paths['output_top_dir'], mainlogfile, 'Writing long_terminal_repeat features:\n{0}\nfrom:\n{1}\nto:\n{2}'.format(','.join(sorted(list(elementSet))),inputGFFpth, outputGFFpth))
	with open(inputGFFpth, 'r') as inGFF:
		currentNewElement = GFF3_line()
		LTR_counts = {}
		for line in inGFF:
			if '\tlong_terminal_repeat\t' in line:
				gffLine = GFF3_line(line)
				if elementSet == None or (elementSet != None and gffLine.attributes['Parent'] in elementSet):
					with open(outputGFFpth, 'a') as outFl:
						if gffLine.attributes['Parent'] in LTR_counts:
							if LTR_counts[gffLine.attributes['Parent']] > 2:
								sys.exit('writeLTRsGFF found element {0} to contain more that 2 LTRs in\n{1}'.format(gffLine.attributes['Parent'], inputGFFpth))
							gffLine.attributes['ID'] = gffLine.attributes['Parent'] + '_R'
							LTR_counts[gffLine.attributes['Parent']] += 1
						else:
							gffLine.attributes['ID'] = gffLine.attributes['Parent'] + '_L'
							LTR_counts[gffLine.attributes['Parent']] = 1
						gffLine.attributes_order = ['ID', 'Parent']
						gffLine.refreshAttrStr()
						outFl.write(str(gffLine)+'\n')
							
def ltrharvest():
	'''
	Runs LTRharvest. LTRharvest options can be specified on the command line.
	See phaltr -h for defaults and phaltr -help for explanation.
	'''
	global paths
	global filenames

	if LTRHARVEST:
		if not 'inputFastaSuffixArray' in paths: # If this is in paths this step has been completed. Skip
			MakeDir('suffixerator_dir', '{0}/suffixerator'.format(paths['output_top_dir']))
			paths['inputFastaSuffixArray'] = '{0}/{1}.index'.format(paths['suffixerator_dir'], paths['inputFasta'])
			gt_suffixerator_call_string = 'gt suffixerator -db {1} -indexname {0} -dna -tis -suf -lcp -des -ssp 1>suffixerator.stdout 2>suffixerator.stderr'.format(paths['inputFastaSuffixArray'], paths['inputFasta'])
			gt_suffixerator_call = [ executables['genometools'], 'suffixerator',  '-db',  paths['inputFasta'], '-indexname', paths['inputFastaSuffixArray'], '-dna', '-tis', '-suf', '-lcp', '-des', '-ssp' ]
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began creating suffix array for  {0} using the call:\n{1}'.format(paths['inputFasta'], gt_suffixerator_call_string))

			# Run suffixerator
			makecall(gt_suffixerator_call, '{0}/suffixerator.stdout'.format(paths['suffixerator_dir']), '{0}/suffixerator.stderr'.format(paths['suffixerator_dir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished gt suffixerator' )
			paths['suffixeratorInputFastaCopy'] = '{0}/{1}'.format(paths['suffixerator_dir'], filenames['inputFasta'])
			copyfile(paths['inputFasta'], paths['suffixeratorInputFastaCopy'])
			# Add suffix array path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('inputFastaSuffixArray\t{0}\n'.format(paths['inputFastaSuffixArray']))

		if not 'LTRharvestGFF' in paths: # If this is in paths this step has been completed. Skip
			# Make dir for LTRharvest
			MakeDir('ltrharvest_dir', '{0}/LTRharvest'.format(paths['output_top_dir']))
			paths['LTRharvestGFF'] = '{0}/{1}.ltrharvest.out.gff'.format(paths['ltrharvest_dir'], filenames['inputFasta'])

			# Run LTRharvest
			gt_ltrharvest_call = [ executables['genometools'], 'ltrharvest', '-similar', str(ltrharvest_similar), '-index', paths['inputFastaSuffixArray'], '-gff3', paths['LTRharvestGFF'], '-seqids', 'yes', '-v', 'yes', '-mintsd', str(ltrharvest_mintsd), '-maxtsd', str(ltrharvest_maxtsd), '-xdrop', str(ltrharvest_xdrop), '-mat', str(ltrharvest_mat), '-mis', str(ltrharvest_mis), '-ins', str(ltrharvest_ins), '-del', str(ltrharvest_del), '-minlenltr', str(ltrharvest_minlenltr), '-maxlenltr', str(ltrharvest_maxlenltr), '-mindistltr', str(ltrharvest_mindistltr), '-maxdistltr', str(ltrharvest_maxdistltr), '-vic', str(ltrharvest_vic) ]
			gt_ltrharvest_call_string = '{0} {1}'.format(' '.join(gt_ltrharvest_call),  '1>ltrharvest.stdout 2>ltrharvest.stderr'.format(paths['inputFastaSuffixArray'], paths['LTRharvestGFF'], executables['genometools']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began running LTRharvest on {0} using the call:\n{1}'.format(paths['inputFasta'], gt_ltrharvest_call_string))
			makecall(gt_ltrharvest_call,  '{0}/ltrharvest.stdout'.format(paths['ltrharvest_dir']), '{0}/ltrharvest.stderr'.format(paths['ltrharvest_dir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished gt ltrharvest' )

			# Need to sort GFF3, sometimes it's not sorted like LTRdigest needs it sorted
			gt_sort_call = [ executables['genometools'], 'gff3', '-sort', '-retainids', paths['LTRharvestGFF'] ]
			gt_sort_call_string = 'gt gff3 -sort -retainids {0} > {0}.sorted 2>{0}.gff3sort.err'.format(paths['LTRharvestGFF'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began sorting LTRharvest:\n{0}'.format(gt_ltrharvest_call_string))
			makecall(gt_sort_call,  stdout='{0}.sorted'.format(paths['LTRharvestGFF']), stderr='{0}.gff3sort.err'.format(paths['LTRharvestGFF']))
			os.rename('{0}.sorted'.format(paths['LTRharvestGFF']), paths['LTRharvestGFF'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished sorting GFF3' )

			# Add LTRharvest GFF3 path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('LTRharvestGFF\t{0}\n'.format(paths['LTRharvestGFF']))
			paths['CurrentGFF'] = paths['LTRharvestGFF']


def ltrdigest():
	'''
	Runs LTRdigest. Domains to search for can be be given as a multi HMM with the command line flag: --ltrdigest_hmms
	'''
	global paths
	global filenames

	if LTRDIGEST: # Identify parts of element internal regions with evidence of homology to LTR RT protein coding domains
		os.environ['PATH'] = '{0}:{1}'.format(executables['hmmer'], os.environ['PATH'])
		if not 'LTRdigestGFF' in paths: # If this is in paths this step has been completed. Skip
			if not 'suffixeratorInputFastaCopy' in paths:
				paths['suffixeratorInputFastaCopy'] = paths['inputFasta']

			MakeDir('ltrdigest_dir', '{0}/LTRdigest'.format(paths['output_top_dir']))
			paths['LTRdigestOutputPrefix'] = '{0}/{1}.LTRdigest'.format(paths['ltrdigest_dir'], paths['inputFasta'])
			paths['LTRdigestGFF'] = '{0}.gff'.format(paths['LTRdigestOutputPrefix'])
			filenames['LTRdigestGFF'] = '{0}.LTRdigest.gff'.format(filenames['inputFasta'])

			gt_ltrdigest_call = [ executables['genometools'], '-j', str(procs), 'ltrdigest', '-matchdescstart', '-outfileprefix', paths['LTRdigestOutputPrefix'], '-hmms', '{0}'.format(paths['LTRdigestHMMs']), '-seqfile', paths['suffixeratorInputFastaCopy']]

			gt_ltrdigest_call_string = '{0} -j {1} ltrdigest -matchdescstart -outfileprefix {2} -hmms {3} -seqfile {4} < {5} > {6}'.format(executables['genometools'], procs, paths['LTRdigestOutputPrefix'], paths['LTRdigestHMMs'], paths['suffixeratorInputFastaCopy'], paths['CurrentGFF'], '{0}.gff'.format(paths['LTRdigestOutputPrefix']))

			append2logfile(paths['output_top_dir'], mainlogfile, 'Began running LTRdigest on {0} using the call:\n{1}'.format(paths['inputFasta'], gt_ltrdigest_call_string))
			makecall(gt_ltrdigest_call, '{0}.gff'.format(paths['LTRdigestOutputPrefix']), '{0}/ltrdigest.stderr'.format(paths['ltrdigest_dir']), paths['CurrentGFF'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished gt ltrdigest' )

			# Add LTRdigest GFF3 path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('LTRdigestGFF\t{0}\n'.format(paths['LTRdigestGFF']))
			paths['CurrentGFF'] = paths['LTRdigestGFF']

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


def bestORFs(fasta, outdir, gff, minLen=300):
	'''
	Finds all ORFs in fasta using EMBOSS getorf and writes the best ones to a GFF3 and protein FASTA.
	The best ORFs are the set of non-overlapping ORFs containing the longest ORF out of sets on that
	strand,
	given by records in gff, or, if strand is unknown, the best set out of all.

	The coordinates output by getorf are 1-based

	Only ORFs with nucleotide sequences longer than minLen are kept. (default 300 bp = 80 aa)
	'''
	outgff = '{0}/{1}.orfs.gff'.format(outdir, fasta.split('/')[-1])
	if os.path.isfile(outgff):
		os.remove(outgff)
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
	if os.path.isfile(outseq):
		os.remove(outseq)
	getorf_call = [ executables['getorf'], '-sequence', fasta, '-outseq', outseq ]
	makecall(getorf_call)

	# Read getorf output and put info into dicts
	orfs = list(SeqIO.parse(outseq, 'fasta'))
	nonredundant_orfs = list()
	unique_orfs = set()
	for seqio in orfs:
		if not seqio.description in unique_orfs:
			unique_orfs.add(seqio.description)
			nonredundant_orfs.append(seqio)
	orfs_ordered_lengths = {}
	orfs_ordered_coords = {}
	orfs_seqs_dct = {}
	orfs_coords = {}
	coords2lenkey = {}
	for orf in nonredundant_orfs:
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

		if element in orfs_coords:
			if strand in orfs_coords[element]:
				orfs_coords[element][strand][orfnum] = (start, end)
			else:
				orfs_coords[element][strand] = {orfnum:(start,end)}
		else:
			orfs_coords[element] = {strand:{orfnum:(start,end)}}

		if element in coords2lenkey:
			if strand in coords2lenkey[element]:
				coords2lenkey[element][strand][(start, end)] = (orfnum, length)
			else:
				coords2lenkey[element][strand] = {(start,end):(orfnum, length)}
		else:
			coords2lenkey[element] = {strand:{(start,end):(orfnum, length)}}

		if element in orfs_ordered_lengths:
			if strand in orfs_ordered_lengths[element]:
				orfs_ordered_lengths[element][strand].append((orfnum, length))
			else:
				orfs_ordered_lengths[element][strand] = [(orfnum, length)]
		else:
			orfs_ordered_lengths[element] = {strand:[(orfnum, length)]}

		if element in orfs_ordered_coords:
			if strand in orfs_ordered_coords[element]:
				orfs_ordered_coords[element][strand].append((start, end))
			else:
				orfs_ordered_coords[element][strand] = [(start, end)]
		else:
			orfs_ordered_coords[element] = {strand:[(start, end)]}

		if element in orfs_seqs_dct:
			if strand in orfs_seqs_dct[element]:
				orfs_seqs_dct[element][strand][orfnum] = seq
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
			if s not in orfs_ordered_lengths[element]:
				continue
			orfs_ordered_lengths[element][s].sort(reverse=True, key=lambda x:x[1])
			orfs_ordered_coords[element][s].sort(key=lambda x:x[0])
			while len(orfs_ordered_lengths[element][s]) > i+1: # orfs_ordered_length is a list that is modified. i gets incremented

				orfnum = orfs_ordered_lengths[element][s][i][0]
				coord = orfs_coords[element][s][orfnum] # coords of the current orf

				j = orfs_ordered_coords[element][s].index(coord) # current largest orf
				k = j+1 # check for overlaps with next in proximity
				# Compare j with successively further away orfs until a non-overlap is reached
				if not k > len(orfs_ordered_coords[element][s])-1:
					J = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][j]]) # corresponding occurence in lengths dict for j, the current longest ORF
					K = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][k]]) # corresponding occurence in lengths dict for k, the current ORF closest to j if moving toward position 0
					while Overlaps( orfs_ordered_coords[element][s][j], orfs_ordered_coords[element][s][k] ):
						# k overlaps j. remove k. because it is shorter than j.
						coord_removed = orfs_ordered_coords[element][s][k]
						orfs_ordered_coords[element][s] = [ item for item in orfs_ordered_coords[element][s] if not item == coord_removed ]
						lenkey = coords2lenkey[element][s][coord_removed]
						orfs_ordered_lengths[element][s].remove(lenkey)
						j = orfs_ordered_coords[element][s].index(coord) # current largest orf
						k = j+1 # check for overlaps with next in proximity
						if k > len(orfs_ordered_coords[element][s])-1:
							break
						J = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][j]]) # corresponding occurence in lengths dict for j, the current longest ORF
						K = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][k]]) # corresponding occurence in lengths dict for k, the current ORF closest to j if moving toward position 0
				# Compare j with successively further away orfs until a non-overlap is reached
				j = orfs_ordered_coords[element][s].index(coord) # current largest orf
				m = j-1 # check for overlaps with previous in proximity
				if m > 0:
					J = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][j]]) # corresponding occurence in lengths dict for j, the current longest ORF
					M = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][m]]) # corresponding occurence in lengths dict for m, the current ORF closest to j if moving toward position 0
					while Overlaps( orfs_ordered_coords[element][s][j], orfs_ordered_coords[element][s][m] ):
						coord_removed = orfs_ordered_coords[element][s][m]
						orfs_ordered_coords[element][s] = [ item for item in orfs_ordered_coords[element][s] if not item == coord_removed ]
						lenkey = coords2lenkey[element][s][coord_removed]
						orfs_ordered_lengths[element][s].remove(lenkey)
						j = orfs_ordered_coords[element][s].index(coord) # current largest orf
						m = j-1 # check for overlaps with previous in proximity
						if m < 0:
							break
						J = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][j]]) # corresponding occurence in lengths dict for j, the current longest ORF
						M = orfs_ordered_lengths[element][s].index(coords2lenkey[element][s][orfs_ordered_coords[element][s][m]]) # corresponding occurence in lengths dict for m, the current ORF closest to j if moving toward position 0
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
			with open(outgff, 'a') as outFl:
				outFl.write('{0}\tgetorf\tORF\t{1}\t{2}\t.\t{3}\t.\tParent={4};translated_seq={5}\n'.format(element, start, end, strand_used, element, seq))


def addORFs(maingff, orfgff, newgff):
	'''
	Inserts ORFs into existing LTRdigest/LTRharvest GFF. Expects Orfs were obtained from EMBOSS getorf on output from writeLTRretrotransposonInternalRegions()
	Existing features take precedence, and if ORFs overlap existing features, those ORFs are not included in the final output, newgff.
	'''
	# Read orf gff store lines in lists in dict with parent as key
	orfs = {}
	if newgff.endswith('.gff'):
		orffasta = '{0}/{1}'.format('/'.join(newgff.split('/')[:-1]), '{0}.ORFs.fasta'.format(newgff.split('/')[-1][:-4]))
	else:
		orffasta = '{0}/{1}'.format('/'.join(newgff.split('/')[:-1]), '{0}.ORFs.fasta'.format(newgff.split('/')[-1]))
	if os.path.isfile(orffasta):
		os.remove(orffasta)

	append2logfile(paths['output_top_dir'], mainlogfile, 'Incorporating ORFs in {0} with GFF {1} in to {2}'.format(orfgff, maingff, newgff))
	with open(orfgff, 'r') as inFl:
		for line in inFl:
			if not line.startswith('#'):
				gffLine = GFF3_line(line)
				parent = gffLine.attributes['Parent']
				if parent in orfs:
					orfs[parent].append(gffLine)
				else:
					orfs[parent] = [gffLine]
	# Read in main gff
	GFFLines = []
	NewGFFLines = []
	with open(maingff, 'r') as inFl:
		for line in inFl:
			if line.startswith('#'):
				continue
			else:
				gffLine = GFF3_line(line)
				GFFLines.append(gffLine)
	internalparts = []
	el = None
	firstLTRend = None
	for i in range(len(GFFLines)):
		gl = GFFLines[i]
		if gl.type == 'repeat_region':
			internalparts = []
			el = 'LTR_retrotransposon' + gl.attributes['ID'][13:]
			NewGFFLines.append(gl)
		elif gl.type == 'target_site_duplication':
			NewGFFLines.append(gl)
		elif gl.type == 'LTR_retrotransposon':
			NewGFFLines.append(gl)
		elif gl.type == 'long_terminal_repeat':
			if firstLTRend != None:
				if el in orfs:
					orf_ct = 0
					for orf in orfs[el]:
						orf.start = firstLTRend + int(orf.start)
						orf.end = firstLTRend + int(orf.end)
						orf.seqid = gffLine.seqid # Change the scaffold name
						OVERLAP = False
						for part in internalparts:
						#	print([orf.start, orf.end], [part.start, part.end])
						#	print(Overlaps([orf.start, orf.end], [part.start, part.end]))
							if Overlaps([orf.start, orf.end], [part.start, part.end]):
								OVERLAP = True
								break
						if not OVERLAP:
							orf_ct += 1
							orf.attributes['ID'] = '{0}.ORF.{1:02d}'.format(orf.attributes['Parent'], orf_ct)
							orf.attributes_order.insert(0, 'ID')
							orf.refreshAttrStr()
							with open(orffasta, 'a') as outFl:
								outFl.write('>{0}\n{1}\n'.format(orf.attributes['ID'], orf.attributes['translated_seq']))
							internalparts.append(orf)
				internalparts.sort(key=lambda x:int(x.start))
				NewGFFLines += internalparts
				NewGFFLines.append(gl)
				firstLTRend = None
			elif firstLTRend == None:
				firstLTRend = int(gl.end)
				NewGFFLines.append(gl)
		else:
			internalparts.append(gl)
	
	with open(newgff, 'w') as outFl:
		outFl.write('##gff-version 3\n')
		for gl in NewGFFLines:
			if gl.type == 'repeat_region':
				outFl.write('###\n')
			outFl.write('{0}\n'.format(str(gl)))


def AnnotateORFs(minLen):
	'''
	Uses bestORFs() and addORFs() to add ORFs of length > minLen
	to the GFF3 if they don't overlap existing features.
	'''
	global paths

	if not checkStatusFl('WithORFsGFF'):
		MakeDir('ORFsDir', '{0}/AnnotateORFs'.format(paths['output_top_dir']))
		internalGFF = '{0}/internals.fasta'.format(paths['ORFsDir'])
		internalFASTA = '{0}/internals.gff'.format(paths['ORFsDir'])
		writeLTRretrotransposonInternalRegions(paths['CurrentGFF'], internalGFF, elementSet=None, truncateParent=False)
		getfasta_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', internalGFF ]
		makecall(getfasta_call, internalFASTA)
		ChangeFastaHeaders(internalFASTA, internalGFF, attribute='Parent')
		bestORFs(fasta=internalFASTA, outdir=paths['ORFsDir'], gff=paths['CurrentGFF'], minLen=minLen)
		orfgff = '{0}/{1}.orfs.gff'.format(paths['ORFsDir'], internalFASTA.split('/')[-1])
		withorfsgff='{0}/FullWithORFs_gt_{1}bp.gff'.format(paths['ORFsDir'], minLen)
		addORFs(maingff=paths['CurrentGFF'], orfgff=orfgff, newgff=withorfsgff)
		paths['WithORFsGFF'] = '{0}/{1}.withORFs_gt_{2}bp.gff'.format(paths['GFFOutputDir'], '.'.join(paths['CurrentGFF'].split('/')[-1].split('.')[:-1]), minLen)
		copyfile(withorfsgff, paths['WithORFsGFF'])
		with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
			statusFlAppend.write('WithORFsGFF\t{0}\n'.format(paths['WithORFsGFF']))
		paths['CurrentGFF'] = paths['WithORFsGFF']


def classify_by_homology(KEEPCONFLICTS=False, KEEPNOCLASSIFICATION=False, repbase_tblastx_evalue=1e-5, nhmmer_reporting_evalue=5e-2, nhmmer_inclusion_evalue=1e-2):

	global paths
	global filenames

	if CLASSIFYDFAM or CLASSIFYREPBASE: # Extract LTR_retrotransposon sequences for classification using homology
		if not 'LTRharvest_LTR_retrotransposons_fasta' in paths: # If this is in paths this step has been completed. Skip
			paths['LTRharvest_LTR_retrotransposons_GFF'] = '{0}/LTRharvest_LTR_retrotransposons.gff'.format(paths['GFFOutputDir'])
			paths['LTRharvest_LTR_retrotransposons_fasta'] = '{0}/LTRharvest_LTR_retrotransposons.fasta'.format(paths['FastaOutputDir'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting LTR_retrotransposons from LTRharvest GFF')
			# Write GFF for just LTR_retrotransposon features
			with open(paths['LTRharvestGFF'], 'r') as harvestFl:
				for line in harvestFl:
					if not line.startswith('#'):
						gffLine = GFF3_line(line)
						if gffLine.type == 'LTR_retrotransposon':
							with open(paths['LTRharvest_LTR_retrotransposons_GFF'], 'a') as LTRharvest_LTR_retrotransposons_GFF:
								LTRharvest_LTR_retrotransposons_GFF.write('{0}\n'.format(str(gffLine)))

			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting LTR_retrotransposons from LTRharvest GFF')
			getfasta_ltrretrotransposons_call = [  executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', '{0}'.format(paths['LTRharvest_LTR_retrotransposons_GFF']) ]
			getfasta_ltrretrotransposons_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['LTRharvest_LTR_retrotransposons_GFF'], paths['LTRharvest_LTR_retrotransposons_fasta'], '{0}/bedtools_getfasta.stderr'.format(paths['FastaOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting LTR_retrotransposon sequences from LTRharvest GFF:\n{0}'.format(getfasta_ltrretrotransposons_call_string))
			makecall(getfasta_ltrretrotransposons_call, paths['LTRharvest_LTR_retrotransposons_fasta'], '{0}/bedtools_getfasta.stderr'.format(paths['FastaOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting LTR_retrotransposon sequences from LTRharvest GFF')
			append2logfile(paths['output_top_dir'], mainlogfile, 'Changing FASTA headers from bedtools getfasta-style to LTR_retrotransposon ID')
			ChangeFastaHeaders(paths['LTRharvest_LTR_retrotransposons_fasta'], paths['LTRharvest_LTR_retrotransposons_GFF'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Done changing FASTA headers from bedtools getfasta-style to LTR_retrotransposon ID')
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('LTRharvest_LTR_retrotransposons_fasta\t{0}\n'.format(paths['LTRharvest_LTR_retrotransposons_fasta']))

	if CLASSIFYDFAM: # Find evidence of homology to repeats in Dfam using HMMER. NEED TO CHANGE THIS SO REVERSE COMPLEMENT IS ALSO SEARCHED (nhmmsearch I think)
		
		paths['DfamDB'] = '{0}/RepeatDatabases/Dfam/Dfam.LTRs.hmm'.format(paths['selfDir'])
		paths['DfamTruePosLTRlist'] = '{0}/RepeatDatabases/Dfam/Dfam.annotations.LTR.names.cleaned.txt'.format(paths['selfDir'])
		paths['DfamShortNames'] = '{0}/RepeatDatabases/Dfam/Dfam.annotations.LTR.names.cleaned.map'.format(paths['selfDir'])
		
		if not 'DfamTable' in paths: # If this is in paths this step has been completed. Skip
			# make Dfam classification output dir
			MakeDir('DfamClassificationDir', '{0}/DfamClassification'.format(paths['output_top_dir']))
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('DfamClassificationDir\t{0}\n'.format(paths['DfamClassificationDir']))

			# run hmmsearch of LTR_retrotransposon features from LTRdigest or LTRharvest on Dfam
			paths['nhmmer_DfamHits_table'] = '{0}/{1}.nhmmer_DfamHits.table'.format(paths['DfamClassificationDir'], filenames['inputFasta'])
			nhmmer_dfam_call = [ '{0}/nhmmer'.format(executables['hmmer']), '--tblout', paths['nhmmer_DfamHits_table'], '--incE', str(nhmmer_inclusion_evalue), '-E', str(nhmmer_reporting_evalue), '--cpu', str(procs), paths['DfamDB'], paths['LTRharvest_LTR_retrotransposons_fasta'] ]
			nhmmer_dfam_call_string = '{0} {1}'.format(' '.join(nhmmer_dfam_call), '1>{0}.nhmmer_DfamHits.stdout 2>{0}.nhmmer_DfamHits.stderr'.format('{0}/{1}'.format(paths['DfamClassificationDir'], filenames['inputFasta'])))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began nhmmer of Dfam:\n{0}'.format(nhmmer_dfam_call_string))
			makecall(nhmmer_dfam_call, '{0}.nhmmer_DfamHits.stdout'.format('{0}/{1}'.format(paths['DfamClassificationDir'], filenames['inputFasta'])), '{0}.nhmmer_DfamHits.stderr'.format('{0}/{1}'.format(paths['DfamClassificationDir'], filenames['inputFasta'])))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished nhmmer of Dfam')

			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('DfamTable\t{0}\n'.format(paths['nhmmer_DfamHits_table']))
			paths['DfamTable'] = paths['nhmmer_DfamHits_table']

		# add Dfam classifications to GFF
		if not 'GFFwithDfamClassification' in paths: # If this is in paths this step has been completed. Skip
			# Extract best hits for each query seq
			paths['DfamResultsTableParsed'] = '{0}/{1}.LTR_retrotransposon_DfamBestHits.tab'.format(paths['DfamClassificationDir'], filenames['inputFasta'])
			dfam_results_parse_call_string = '{0}/nhmmer_table2columns.py < {1} > {2} 2>{3}/nhmmer_table2columns.py.stderr'.format(paths['scriptsDir'], paths['DfamTable'], paths['DfamResultsTableParsed'], paths['DfamClassificationDir'])
			dfam_results_parse_call = ['{0}/nhmmer_table2columns.py'.format(paths['scriptsDir']) ]
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting best hits from  nhmmer on Dfam results:\n{0}'.format(dfam_results_parse_call_string))
			makecall(dfam_results_parse_call, paths['DfamResultsTableParsed'], '{0}/nhmmer_table2columns.py.stderr'.format(paths['DfamClassificationDir']), stdin=paths['DfamTable'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting best hits from  nhmmer on Dfam results')

			# Add best hits to GFF
			paths['GFFwithDfamClassification'] = '{0}/{1}.LTRdigest.withDfam.gff'.format(paths['GFFOutputDir'], filenames['inputFasta'])
			add_dfam_hits_to_ltrdigest_gff_call_string = '{0}/gffAddAttr.py -gff {1} -attr dfamClassification -map {2} -mapKey ID -restrictType LTR_retrotransposon -replaceIfNone > {3} 2>{4}/gffAddAttr.py.DfamHits.stderr'.format(paths['scriptsDir'], paths['CurrentGFF'], paths['DfamResultsTableParsed'], paths['GFFwithDfamClassification'], paths['GFFOutputDir'])
			add_dfam_hits_to_ltrdigest_gff_call = [ '{0}/gffAddAttr.py'.format(paths['scriptsDir']), '-gff', paths['CurrentGFF'], '-attr', 'dfamClassification', '-map', paths['DfamResultsTableParsed'], '-mapKey', 'ID', '-restrictType', 'LTR_retrotransposon', '-replaceIfNone' ]
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began adding best hits from nhmmer on Dfam results to LTRdigest GFF:\n{0}'.format(add_dfam_hits_to_ltrdigest_gff_call_string))
			makecall(add_dfam_hits_to_ltrdigest_gff_call, paths['GFFwithDfamClassification'], '{0}/gffAddAttr.py.DfamHits.stderr'.format(paths['GFFOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished adding best hits from nhmmer on Dfam results to LTRdigest GFF')

			append2logfile(paths['output_top_dir'], mainlogfile, 'Update strandedness in GFF3 based on Dfam results')

			# Add LTRdigest GFF3 with Dfam classifications path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('GFFwithDfamClassification\t{0}\n'.format(paths['GFFwithDfamClassification']))
			paths['CurrentGFF'] = paths['GFFwithDfamClassification']

	if CLASSIFYREPBASE: # Find evidence of homology to repeats in Repbase using tblastx

		paths['RepbaseDB'] = '{0}/RepeatDatabases/Repbase/Repbase22.04_LTRs.fasta'.format(paths['selfDir'])
		paths['RepbaseTruePosLTRlist'] = '{0}/RepeatDatabases/Repbase/Repbase.annotations.LTR.names.cleaned.txt'.format(paths['selfDir'])
		paths['RepbaseShortNames'] = '{0}/RepeatDatabases/Repbase/Repbase.annotations.LTR.names.cleaned.map'.format(paths['selfDir'])
		os.environ['BLASTDB'] = '{0}/RepeatDatabases/Repbase/:{0}/{1}'.format(paths['selfDir'],paths['FastaOutputDir'])

		if not 'RepbaseTable' in paths: # If this is in paths this step has been completed. Skip
			# make Repbase annotation dir
			MakeDir('RepbaseClassificationDir', '{0}/RepbaseClassification'.format(paths['output_top_dir']))
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('RepbaseClassificationDir\t{0}\n'.format(paths['RepbaseClassificationDir']))

			# run tblastx of LTR_retrotransposon features from LTRdigest or LTRharvest on Repbase
			paths['tblastx_RepbaseHits_table'] = '{0}/{1}.tblastx_Repbase.tab'.format(paths['RepbaseClassificationDir'], filenames['inputFasta'])
			tblastx_repbase_call = [ '{0}/tblastx'.format(executables['blast']), '-db', 'Repbase22.04_LTRs.fasta', '-query', paths['LTRharvest_LTR_retrotransposons_fasta'], '-evalue', str(repbase_tblastx_evalue), '-outfmt', '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand', '-num_threads', str(procs), '-max_hsps', '25' ]
			tblastx_repbase_call_string = '{0} -db {1} -query {2} -evalue {5} -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -num_threads {3} -max_hsps 25 1>{4} 2>{4}.stderr'.format('{0}/tblastx'.format(executables['blast']), 'Repbase22.04_LTRs.fasta', paths['LTRharvest_LTR_retrotransposons_fasta'], procs, paths['tblastx_RepbaseHits_table'], repbase_tblastx_evalue)
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began tblastx of Repbase:\n{0}'.format(tblastx_repbase_call_string))
			makecall(tblastx_repbase_call, paths['tblastx_RepbaseHits_table'], '{0}.stderr'.format(paths['tblastx_RepbaseHits_table']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished tblastx of Repbase')

			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('RepbaseTable\t{0}\n'.format(paths['tblastx_RepbaseHits_table']))
			paths['RepbaseTable'] = paths['tblastx_RepbaseHits_table']

		# add Repbase classifications to GFF
		if not 'GFFwithRepbaseClassification' in paths: # If this is in paths this step has been completed. Skip

			# Extract best hits for each query seq
			paths['RepbaseResultsTableParsed'] = '{0}/{1}.LTR_retrotransposon_RepbaseBestHits.tab'.format(paths['RepbaseClassificationDir'], filenames['inputFasta'])
			repbase_results_parse_call_string = '{0}/best_blast_hit.py < {1} > {2} 2>{3}/best_blast_hit.py.stderr'.format(paths['scriptsDir'], paths['RepbaseTable'], paths['RepbaseResultsTableParsed'], paths['RepbaseClassificationDir'])
			repbase_results_parse_call = ['{0}/best_blast_hit.py'.format(paths['scriptsDir']) ]
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting best hits from  tblastx on Repbase results:\n{0}'.format(repbase_results_parse_call_string))
			makecall(repbase_results_parse_call, paths['RepbaseResultsTableParsed'], '{0}/best_blast_hit.py.stderr'.format(paths['RepbaseClassificationDir']), stdin=paths['RepbaseTable'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting best hits from  tblastx on Repbase results')

			# Add best hits to GFF
			if CLASSIFYDFAM:
				gff_for_repbase_classification = paths['GFFwithDfamClassification']
				paths['GFFwithRepbaseClassification'] = '{0}/{1}.LTRdigest.withDfam.withRepbase.gff'.format(paths['GFFOutputDir'], filenames['inputFasta'])
			else:
				gff_for_repbase_classification = paths['LTRdigestGFF']
				paths['GFFwithRepbaseClassification'] = '{0}/{1}.LTRdigest.withRepbase.gff'.format(paths['GFFOutputDir'], filenames['inputFasta'])

			add_repbase_hits_to_ltrdigest_gff_call_string = '{0}/gffAddAttr.py -gff {1} -attr repbaseClassification -map {2} -mapKey ID -restrictType LTR_retrotransposon -replaceIfNone > {3} 2>{4}/gffAddAttr.py.RepbaseHits.stderr'.format(paths['scriptsDir'], gff_for_repbase_classification, paths['RepbaseResultsTableParsed'], paths['GFFwithRepbaseClassification'], paths['GFFOutputDir'])
			add_repbase_hits_to_ltrdigest_gff_call = [ '{0}/gffAddAttr.py'.format(paths['scriptsDir']), '-gff', gff_for_repbase_classification, '-attr', 'repbaseClassification', '-map', paths['RepbaseResultsTableParsed'], '-mapKey', 'ID', '-restrictType', 'LTR_retrotransposon', '-replaceIfNone' ]
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began adding best hits from tblastx on Repbase results to LTRdigest GFF:\n{0}'.format(add_repbase_hits_to_ltrdigest_gff_call_string))
			makecall(add_repbase_hits_to_ltrdigest_gff_call, paths['GFFwithRepbaseClassification'], '{0}/gffAddAttr.py.RepbaseHits.stderr'.format(paths['GFFOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished adding best hits from tblastx on Repbase results to LTRdigest GFF')

			# Add LTRdigest GFF3 with Repbase classifications path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('GFFwithRepbaseClassification\t{0}\n'.format(paths['GFFwithRepbaseClassification']))
			paths['CurrentGFF'] = paths['GFFwithRepbaseClassification']

	if CLASSIFYDFAM or CLASSIFYREPBASE: # Remove false positives from LTRdigest GFF3

		if not 'LTRdigestClassifiedNoFP' in paths: # If this is in paths this step has been completed. Skip

			if CLASSIFYREPBASE:
				gff_classified = paths['GFFwithRepbaseClassification'] # Will have both Dfam and Repbase classifications if CLASSIFYDFAM==True also
			else:
				gff_classified = paths['GFFwithDfamClassification']

			paths['LTRdigestClassifiedNoFP'] = '{0}/{1}.LTRdigestClassifiedNoFP.gff'.format(paths['GFFOutputDir'], filenames['inputFasta'])
			TruePositiveLTRclassificationsDct = { 'dfamClassification':paths['DfamTruePosLTRlist'], 'repbaseClassification':paths['RepbaseTruePosLTRlist'] }
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began removing false positives from LTRdigest GFF with classifications.')
			RemoveNonLTRretrotransposons(gff_classified, TruePositiveLTRclassificationsDct, outputFlName=paths['LTRdigestClassifiedNoFP'], REPORTCONFLICTS=True, KEEPCONFLICTS=KEEPCONFLICTS, KEEPNOCLASSIFICATION=KEEPNOCLASSIFICATION, logFilePth='{0}/RemoveNonLTRretrotransposons.log'.format(paths['GFFOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished removing false positives from LTRdigest GFF with classifications. TP file at:\n{0}'.format(paths['LTRdigestClassifiedNoFP']))
			# Add LTRdigest GFF3 with FP removed path to status file (for resuming later)
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('LTRdigestClassifiedNoFP\t{0}\n'.format(paths['LTRdigestClassifiedNoFP']))
			paths['CurrentGFF'] = paths['LTRdigestClassifiedNoFP']


def shortClassif(ElNames=False):
	'''
	- Opens Dfam and Repbase short, class level names for elements, whose paths are hardcoded
	  in this code to positions relative to this file.
	- Opens GFF3 and assigns classification based on attributes. Conflicting attributes leads
	  to the assignment being 'Other'
	#- Returns a dictionary with element name as keys and assignment as values.
	- Returns a dictionary with assignments as keys and a list of elements as values.
	- ElNames=True means this function will return a dict with LTR RT names as keys.
	- ElNames=False means this function will return a dict with classifs as keys and sets of LTR RT names as values.
	'''
	DfamNames = {}
	with open(paths['DfamShortNames']) as DfamFl:
		for line in DfamFl:
			full, short = line.strip().split()
			DfamNames[full] = short
	RepbaseNames = {}
	with open(paths['RepbaseShortNames']) as RepbaseFl:
		for line in RepbaseFl:
			full, short = line.strip().split()
			RepbaseNames[full] = short
	ElementNames = {}
	with open(paths['CurrentGFF']) as gffFl:
		for line in gffFl:
			if not line.startswith('#'):
				if '\tLTR_retrotransposon\t' in line:
					gffLine = GFF3_line(line)
					el = gffLine.attributes['ID']
					for attr in gffLine.attributes:
						if 'dfam' in attr:
							if gffLine.attributes[attr] in DfamNames:
								shortName = DfamNames[gffLine.attributes[attr]]
								if shortName != 'None':
									if el in ElementNames:
										currentName = ElementNames[el]
										if currentName == shortName:
											continue
										else:
											if currentName == 'Other':
												ElementNames[el] = shortName
											else:
												if el != 'Other':
													ElementNames[el] = 'Other' # If classifications disagree, other is assigned as the short classification
												else:
													continue # If competing classifs are 'Other' and something else, 'Other' is not given priority
									else:
										ElementNames[el] = shortName
									
							else:
								if el in ElementNames:
									continue
								else:
									ElementNames[el] = 'Other'

						elif 'repbase' in attr:
							if gffLine.attributes[attr] in RepbaseNames:
								shortName = RepbaseNames[gffLine.attributes[attr]]
								if shortName != 'None':
									if el in ElementNames:
										currentName = ElementNames[el]
										if currentName == shortName:
											continue
										else:
											if currentName == 'Other':
												ElementNames[el] = shortName
											else:
												if el != 'Other':
													ElementNames[el] = 'Other' # If classifications disagree, other is assigned as the short classification
												else:
													continue # If competing classifs are 'Other' and something else, 'Other' is not given priority
									else:
										ElementNames[el] = shortName
							else:
								if el in ElementNames:
									continue
								else:
									ElementNames[el] = 'Other'
	# return 1
	if ElNames:
		return ElementNames

	Classifications = {}
	for el, clasif in ElementNames.items():
		if clasif in Classifications:
			Classifications[clasif].add(el)
		else:
			Classifications[clasif] = set([el])
	# Remove any existing GFF3s to avoid appending to them. May want to change so if a file exists it is read and reclassifying is skipped instead of overwritten.
	classifs = set(list(DfamNames.values()) + list(RepbaseNames.values()))
	for classif in classifs:
		currentGFFpath = '{0}/{1}.gff'.format(paths['GFFOutputDir'], classif)
		if os.path.isfile(currentGFFpath):
			os.remove(currentGFFpath)

	for clasif in Classifications:
		sortedClassifs = sorted(list(Classifications[clasif]))
		classifGFFpath = '{0}/{1}.gff'.format(paths['GFFOutputDir'], clasif)
		with open(classifGFFpath, 'a') as gffFl:
			for el in sortedClassifs:
				gffFl.write('{0}\n'.format(el))
	# return 2
	return Classifications


def checkStatusFl(key):
	'''
	The 'status' file contains paths to files and directories
	to help the program resume at points it has completed.
	'''
	with open('{0}/status'.format(paths['output_top_dir']), 'r') as statusFl:
		for line in statusFl:
			if line.strip().split('\t')[0] == key:
				return True
		return False


def graph2groups(G):
	'''
	Returns lists of each element in each connected
	component in G as dictionary. DFS recursive algorithm.
	Adapted from: https://stackoverflow.com/questions/21078445/find-connected-components-in-a-graph
	'''
	def dfs(node1):

		nonlocal visited
		nonlocal groups
		nonlocal group
		nonlocal G

		for node2 in G[node1]:
			if not node2 in visited:
				visited.add(node2)
				groups[group].append(node2)
				dfs(node2)
	visited = set()
	group = 0
	groups = {}
	for node in G:
		if node not in visited:
			visited.add(node)
			group += 1
			groups[group] = [node]
			dfs(node)
	return groups


def WickerFam(pId=80, percAln=80, minLen=80, use_ltrs=True, use_internal=True):
	'''
	Creates a classification of elements using the Wicker et al. 2007 protocol:

	"80% sequence similarity or more in at least 80% of the aligned sequence.
	Two elements belong to the same family if they:
		1. share 80% (or more) sequence identity
		2. in at least 80% of their coding or internal domain, or within their terminal repeat
	regions, or in both.
	
	To avoid misclassification of short and possibly random stretches of homology,
	we recommend analysing only segments of longer than 80 bp. TEs that are smaller than
	80 bp require specialized analyses.

	The terminal repeat regions and other non-coding regions are the fastest evolving parts of
	TEs. Therefore, they offer the most specificity in defining families. Allowing the 80-80-80
	rule for DNA sequence identity in either the internal domain or in the terminal regions, or
	both, also addresses the problem caused by frequent TE truncations. In some cases only terminal
	repeats and non-coding regions may be present, whereas in other cases only parts of the coding
	region but no terminal repeats may be available for analysis.
	
	In some cases, it may be necessary to add the subfamily taxon, depending on the population
	structure of a TE family. Subfamilies can be populations of non-autonomous deletion derivatives
	or distinct subpopulations in large families that can be clearly segregated. The similarity
	threshold can differ between subfamilies, depending on the number and homogeneity of elements
	described. However, such distinctions are matters for TE specialists and should not be a burden
	for annotators (see the wikiPoson web site for further discussion). Importantly, the term should
	not be used for groupings above the family level."


	LTRs and internal regions will compared separately

	I. LTRs:
		1. Extract all LTRs and run blastn all by all (runblast())
		2. Parse blast output for pairs that satisfy the parameters
	
	II. Internals:
		1. Extract all internals and run blastn all by all (runblast())
		2. Parse blast output for pairs that satisfy the parameters
	
	III. Assign
		1. Write family designations (Like mcl output)
	'''
	global paths
	
	WICKERLTRS = use_ltrs
	WICKERINTERNAL = use_internal
	MakeDir('WickerFamDir', '{0}/WickerFamDir'.format(paths['output_top_dir']))
	time.sleep(.3)
	MakeDir('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, minLen), '{0}/{1}_pId_{2}_percAln_{3}_minLen'.format(paths['WickerFamDir'], pId, percAln, minLen))
	OutDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, minLen)]
	OutDirKey = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, minLen)

	if not checkStatusFl(OutDirKey):
		with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
			statusFlAppend.write("{0}\t{1}\n".format('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, minLen), paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, minLen)]))

	MakeDir('WickerFam_Cluster_dir', '{0}/Clusters'.format(OutDir))
	OutDir = paths['WickerFam_Cluster_dir']
	for classif in clusters_by_classif:

		if not checkStatusFl('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif)):

			MakeDir('WickerFam_{0}_dir'.format(classif), '{0}/{1}'.format(OutDir, classif))
			cOutDir = paths['WickerFam_{0}_dir'.format(classif)]
			G = {} # Dictionary representation of a graph that will hold the blast results

			if WICKERINTERNAL:
				# Extract internal regions of selected elements
				paths['InternalelementsGFF'] = '{0}/internals.gff'.format(cOutDir)
				paths['InternalelementsFasta'] = '{0}/internals.fasta'.format(cOutDir)
				writeLTRretrotransposonInternalRegions(paths['CurrentGFF'], paths['InternalelementsGFF'], elementSet=set(clusters_by_classif[classif]), truncateParent=False)
				getfasta_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', paths['InternalelementsGFF'] ]
				getfasta_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['InternalelementsGFF'], paths['InternalelementsFasta'], '{0}/bedtools_getfasta_InternalRegions.stderr'.format(OutDir))
				makecall(getfasta_call, paths['InternalelementsFasta'], '{0}/bedtools_getfasta_InternalRegions.stderr'.format(cOutDir))
				ChangeFastaHeaders(paths['InternalelementsFasta'], paths['InternalelementsGFF'], attribute='Parent')

				# Get internals seq lengths
				internal_seq_lengths = { s.id:len(s) for s in list(SeqIO.parse(paths['InternalelementsFasta'], 'fasta')) }
				elements = list(internal_seq_lengths.keys())

				# blast all by all internals
				paths['Internals_{0}_selfBlastnOut'.format(classif)] = '{0}/Internals_selfBlastn.tab'.format(cOutDir)
				runblast(query=paths['InternalelementsFasta'], subject=paths['InternalelementsFasta'], out=paths['Internals_{0}_selfBlastnOut'.format(classif)], evalue='1e-5', outfmt='7', percid=pId, blast='blastn', procs=procs)

				# Parse blast hits to define families
				InternalsBlastPth = paths['Internals_{0}_selfBlastnOut'.format(classif)]
				with open(InternalsBlastPth, 'r') as blastFl:
					internal_aln_lens = {}
					for line in blastFl:
						if line.startswith('#'):
							continue

						query, subj, percent_id, alignment_len, mismatches, gap_opens, q_start, q_end, s_start, s_end, E_value, bit_score = line.strip().split('\t')
						# 1. Ignore self-hits
						if query == subj:
							continue
						# 2. Ignore any alignment less than minLen
						if int(alignment_len) < minLen:
							continue
						# 3. Build dict of longest aln len between two elements
						aln_pair = frozenset([query, subj])
						if aln_pair in internal_aln_lens:
							if int(alignment_len) > internal_aln_lens[aln_pair]:
								internal_aln_lens[aln_pair] = int(alignment_len) 
						else:
							internal_aln_lens[aln_pair] = int(alignment_len) 

					# 4. Ignore any alignment with alnLn[i]/seqLen[i][j] < percAln for each seq in each alignment
					for aln_pair in internal_aln_lens:
						aln_pair_lst = list(aln_pair)
						el1 =  aln_pair_lst[0]
						el2 =  aln_pair_lst[1]
						el1_aln_ratio = internal_aln_lens[aln_pair]/internal_seq_lengths[el1]*100
						el2_aln_ratio = internal_aln_lens[aln_pair]/internal_seq_lengths[el2]*100
						if el1_aln_ratio < percAln or el2_aln_ratio < percAln:
							continue
						# 4. Build graph 
						if el1 in G:
							G[el1].add(el2)
						else:
							G[el1] = set([el2])
						if el2 in G:
							G[el2].add(el1)
						else:
							G[el2] = set([el1])

			if WICKERLTRS:
				# Extract all LTRs from CurrentGFF separated by classif
				paths['LTRs_{0}_GFF'.format(classif)] = '{0}/LTRs_{1}.gff'.format(cOutDir, classif)
				paths['LTRs_{0}_Fasta'.format(classif)] = '{0}/LTRs_{1}.fasta'.format(cOutDir, classif)
				writeLTRsGFF(paths['CurrentGFF'], paths['LTRs_{0}_GFF'.format(classif)], elementSet=clusters_by_classif[classif])
				getfasta_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', paths['LTRs_{0}_GFF'.format(classif)] ]
				getfasta_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['LTRs_{0}_GFF'.format(classif)], paths['LTRs_{0}_Fasta'.format(classif)], '{0}/bedtools_getfasta_InternalRegions.stderr'.format(cOutDir))
				makecall(getfasta_call, stdout=paths['LTRs_{0}_Fasta'.format(classif)], stderr='{0}/bedtools_getfasta_InternalRegions.stderr'.format(cOutDir))
				ChangeFastaHeaders(paths['LTRs_{0}_Fasta'.format(classif)], paths['LTRs_{0}_GFF'.format(classif)], attribute='ID')

				# Get LTRs seq lengths
				ltrs_seq_lengths = { s.id:len(s) for s in list(SeqIO.parse(paths['LTRs_{0}_Fasta'.format(classif)], 'fasta')) }
				element_combined_ltr_lengths = {}
				for ltr in ltrs_seq_lengths:
					el = ltr[:-2] # trim off the _L or _R
					if el in element_combined_ltr_lengths:
						element_combined_ltr_lengths[el] += ltrs_seq_lengths[ltr]
					else:
						element_combined_ltr_lengths[el] = ltrs_seq_lengths[ltr]
				elements = list(element_combined_ltr_lengths.keys())
				# blast all by all LTRs
				paths['LTRs_{0}_selfBlastnOut'.format(classif)] = '{0}/LTRs_selfBlastn.tab'.format(cOutDir)
				runblast(query=paths['LTRs_{0}_Fasta'.format(classif)], subject=paths['LTRs_{0}_Fasta'.format(classif)], out=paths['LTRs_{0}_selfBlastnOut'.format(classif)], outfmt='7', evalue='1e-5', percid=pId, blast='blastn', procs=procs)

				# Parse blast hits to define families
				LTRsBlastPth = paths['LTRs_{0}_selfBlastnOut'.format(classif)]
				with open(LTRsBlastPth, 'r') as blastFl:
					ltr_aln_lens = {}
					for line in blastFl:
						if line.startswith('#'):
							continue

						query, subj, percent_id, alignment_len, mismatches, gap_opens, q_start, q_end, s_start, s_end, E_value, bit_score = line.strip().split('\t')

						# 1. Ignore self-hits
						if query == subj:
							continue
						# 2. Ignore any alignment less than minLen
						if int(alignment_len) < minLen:
							continue
						# 3. Ignore any alignment with alnLn[i]/seqLen[i][j] < percAln for each seq in each alignment
						# 3. Build dict of longest aln len between two elements
						aln_pair = frozenset([query, subj])
						q_el = query[:-2]
						s_el = subj[:-2]
						if q_el in ltr_aln_lens:
							if s_el in ltr_aln_lens[q_el]:
								if aln_pair in ltr_aln_lens[q_el][s_el]:
									if int(alignment_len) > ltr_aln_lens[q_el][s_el][aln_pair]:
										ltr_aln_lens[q_el][s_el][aln_pair] = int(alignment_len) 
								else:
									ltr_aln_lens[q_el][s_el][aln_pair] = int(alignment_len) 
							else:
								ltr_aln_lens[q_el][s_el] = {aln_pair:int(alignment_len)}
						else:
							ltr_aln_lens[q_el] = {s_el:{aln_pair:int(alignment_len)}}

						if s_el in ltr_aln_lens:
							if q_el in ltr_aln_lens[s_el]:
								if aln_pair in ltr_aln_lens[s_el][q_el]:
									if int(alignment_len) > ltr_aln_lens[s_el][q_el][aln_pair]:
										ltr_aln_lens[s_el][q_el][aln_pair] = int(alignment_len) 
								else:
									ltr_aln_lens[s_el][q_el][aln_pair] = int(alignment_len) 
							else:
								ltr_aln_lens[s_el][q_el] = {aln_pair:int(alignment_len)}
						else:
							ltr_aln_lens[s_el] = {q_el:{aln_pair:int(alignment_len)}}

					for el1 in ltr_aln_lens:
						for el2 in ltr_aln_lens[el1]:
							pairs = ltr_aln_lens[el1][el2]
							kinds = {}
							for pair in pairs:
								pLst = sorted(list(pair))
								if pLst[0].endswith('_L') and pLst[1].endswith('_L'):
									if 'LL' in kinds:
										if pairs[pair] > kinds['LL']:
											kinds['LL'] = pairs[pair]
									else:
										kinds['LL'] = pairs[pair]

								elif pLst[0].endswith('_R') and pLst[1].endswith('_R'):
									if 'RR' in kinds:
										if pairs[pair] > kinds['RR']:
											kinds['RR'] = pairs[pair]
									else:
										kinds['RR'] = pairs[pair]
								elif pLst[0].endswith('_L') and pLst[1].endswith('_R'):
									if 'LR' in kinds:
										if pairs[pair] > kinds['LR']:
											kinds['LR'] = pairs[pair]
									else:
										kinds['LR'] = pairs[pair]
								elif pLst[0].endswith('_R') and pLst[1].endswith('_L'):
									if 'RL' in kinds:
										if pairs[pair] > kinds['RL']:
											kinds['RL'] = pairs[pair]
									else:
										kinds['RL'] = pairs[pair]
								else:
									sys.exit('Uknown element name LTR suffixes: {0}'.format(pair))

							len_both_ltrs = [element_combined_ltr_lengths[el1], element_combined_ltr_lengths[el2]]
							PASS = False
							for l in len_both_ltrs:
								ratio1 = 0
								ratio2 = 0
								ratio = 0
								if 'LL' in kinds and 'RR' in kinds:
									ratio1 = ((kinds['LL'] + kinds['RR'])/l)*100
								if 'LR' in kinds and 'RL' in kinds:
									ratio2 = ((kinds['LR'] + kinds['RL'])/l)*100
								if not ratio1 == 0 and ratio2 == 0:
									ratio = max([ratio1, ratio2])
								if ratio >= percAln:
									PASS = True
							if PASS:
								# 4. Build graph 
								if el1 in G:
									G[el1].add(el2)
								else:
									G[el1] = set([el2])
								if el2 in G:
									G[el2].add(el1)
								else:
									G[el2] = set([el1])

			# 5. Report nodes in each connected component as each family
			if G != {}:
				wicker_groups = graph2groups(G)
				group_num = max(list(wicker_groups.keys()))
				non_singletons = set()
				for g in wicker_groups:
					for el in wicker_groups[g]:
						non_singletons.add(el)

				for el in elements:
					if el not in non_singletons:
						group_num += 1
						wicker_groups[group_num] = [el]

				#wicker_group_len = sum([ len(wicker_groups[g]) for g in wicker_groups ])
				paths['Wicker_{0}'.format(classif)] = '{0}/wicker_groups_{1}'.format(cOutDir, classif)
				paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif)] = paths['Wicker_{0}'.format(classif)]
				group_lens = sorted([(g,len(wicker_groups[g])) for g in wicker_groups], key=lambda x:x[1], reverse=True)
				# Write families
				with open(paths['Wicker_{0}'.format(classif)], 'w') as outFl:
					for group in group_lens:
						group = group[0]
						outFl.write('\t'.join(wicker_groups[group])+'\n')
					#outFl.write('\n'.join(['\t'.join(wicker_groups[g]) for g in wicker_groups]))
				if not checkStatusFl('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif)):
					with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
						statusFlAppend.write("{0}\t{1}\n".format('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif), paths['Wicker_{0}'.format(classif)]))
			else:
				# Write families: all singletons
				paths['Wicker_{0}'.format(classif)] = '{0}/wicker_groups_{1}'.format(cOutDir, classif)
				paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif)] = paths['Wicker_{0}'.format(classif)]
				if os.path.isfile(paths['Wicker_{0}'.format(classif)]):
					os.remove(paths['Wicker_{0}'.format(classif)])
				for el in elements:
					with open(paths['Wicker_{0}'.format(classif)], 'a') as outFl:
						outFl.write('{0}\n'.format(el))
				
				if not checkStatusFl('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif)):
					with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
						statusFlAppend.write("{0}\t{1}\n".format('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(pId, percAln, minLen, classif), paths['Wicker_{0}'.format(classif)]))


def MCL(I=6, minClustSize=30):
	'''
	If there are less than minClustSize elements for a given classification or for all
	classifications combined, they are put into one cluster and MCL is not used.
	'''
	global paths
	global filenames

	if USEMCL:
		if not 'LTRdigest_LTR_retrotransposons_fasta' in paths:

			MakeDir('MCLdir', '{0}/MCL'.format(paths['output_top_dir'], I))
			paths['LTRdigest_LTR_retrotransposons_GFF'] = '{0}/LTRdigest_LTR_retrotransposons.gff'.format(paths['GFFOutputDir'])
			paths['LTRdigest_LTR_retrotransposons_fasta'] = '{0}/LTRdigest_LTR_retrotransposons.fasta'.format(paths['FastaOutputDir'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting LTR_retrotransposons from {0}'.format(paths['CurrentGFF']))

			# Write GFF for just LTR_retrotransposon features
			with open(paths['CurrentGFF'], 'r') as gffFl:
				for line in gffFl:
					if not line.startswith('#'):
						gffLine = GFF3_line(line)
						if gffLine.type == 'LTR_retrotransposon':
							with open(paths['LTRdigest_LTR_retrotransposons_GFF'], 'a') as LTRdigest_LTR_retrotransposons_GFF:
								LTRdigest_LTR_retrotransposons_GFF.write('{0}\n'.format(str(gffLine)))

			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting LTR_retrotransposons from {0}'.format(paths['CurrentGFF']))

			# Extract sequences for true positive LTR_retrotransposon features
			getfasta_ltrretrotransposons_call = [  executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', '{0}'.format(paths['LTRdigest_LTR_retrotransposons_GFF']) ]
			getfasta_ltrretrotransposons_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['LTRdigest_LTR_retrotransposons_GFF'], paths['LTRdigest_LTR_retrotransposons_fasta'], '{0}/bedtools_getfasta.stderr'.format(paths['FastaOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Began extracting LTR_retrotransposon sequences from LTRdigest_TruePositives GFF:\n{0}'.format(getfasta_ltrretrotransposons_call_string))
			makecall(getfasta_ltrretrotransposons_call, paths['LTRdigest_LTR_retrotransposons_fasta'], '{0}/bedtools_getfasta.stderr'.format(paths['FastaOutputDir']))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Finished extracting LTR_retrotransposon sequences from LTRdigest_TruePositives GFF')
			append2logfile(paths['output_top_dir'], mainlogfile, 'Changing FASTA headers from bedtools getfasta-style to LTR_retrotransposon ID')
			ChangeFastaHeaders(paths['LTRdigest_LTR_retrotransposons_fasta'], paths['LTRdigest_LTR_retrotransposons_GFF'])
			append2logfile(paths['output_top_dir'], mainlogfile, 'Done changing FASTA headers from bedtools getfasta-style to LTR_retrotransposon ID')

			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('LTRdigest_LTR_retrotransposons_fasta\t{0}\n'.format(paths['LTRdigest_LTR_retrotransposons_fasta'])) 

		AllFASTA = list(SeqIO.parse(paths['LTRdigest_LTR_retrotransposons_fasta'], 'fasta'))
		classifFastas = {}

		MakeDir('MCLdir', '{0}/MCL'.format(paths['output_top_dir'], I))
		# If all elements combined are less than the min clust size specified by the user, then all elements are put into one cluster.
		if sum([len(clusters_by_classif[c]) for c in clusters_by_classif]) < minClustSize:
			classif = 'All'
			allClassifs = '{0}/out.allClust'.format(paths['MCLdir'])
			c = []
			for cls in classifs:
				c += clusters_by_classif[cls]
				
			with open(allClassifs, 'w') as outFl:
				outFl.write('{0}\n'.format('\t'.join(c)))
			
			newRecord = 'MCL_{0}_I{1}'.format(classif, I)
			if not checkStatusFl(newRecord):
				paths['MCL_{0}_I{1}'.format(classif, I)] = '{0}/out.allClust'.format(paths['MCLdir'])
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format('MCL_{0}_I{1}'.format(classif, I), paths['MCL_{0}_I{1}'.format(classif, I)])) 
					return

		MakeDir('MCL_I{0}'.format(I), '{0}/I{1}'.format(paths['MCLdir'], I))
		if not checkStatusFl('MCL_I{0}'.format(I)):
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format('MCL_I{0}'.format(I), '{0}/I{1}'.format(paths['MCLdir'], I))) 
		MakeDir('MCL_I{0}_ClustersDir'.format(I), '{0}/Clusters'.format(paths['MCL_I{0}'.format(I)]))
		for classif in classifs:
			if not checkStatusFl('MCL_{0}_I{1}'.format(classif, I)):
				MakeDir('MCL_{0}_I{1}_dir'.format(classif, I), '{0}/{1}'.format(paths['MCL_I{0}_ClustersDir'.format(I)], classif))
				outputPth = paths['MCL_{0}_I{1}_dir'.format(classif, I)]
				if len(clusters_by_classif[classif]) < minClustSize:
					with open('{0}/{1}_MCL_clusters.I{2}'.format(outputPth, classif, I), 'w') as outFl:
						outFl.write('{0}\n'.format('\t'.join(clusters_by_classif[classif])))

					if 'MCL_{0}_I{1}'.format(classif, I) not in paths:
						paths['MCL_{0}_I{1}'.format(classif, I)] = '{0}/{1}_MCL_clusters.I{2}'.format(outputPth, classif, I)
						with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
							statusFlAppend.write('{0}\t{1}\n'.format('MCL_{0}_I{1}'.format(classif, I), paths['MCL_{0}_I{1}'.format(classif, I)])) 

				if not 'MCL_{0}_abc'.format(classif) in paths:
					classifFastas[classif] = '{0}/{1}.fasta'.format(outputPth, classif)
					SeqIO.write([ rec for rec in AllFASTA if rec.id in clusters_by_classif[classif] ], classifFastas[classif], 'fasta')
					# make blast db for all-by-all blast
					makeblastdb_AllByAll_call_string = '{0}/makeblastdb -in {1} -dbtype nucl'.format(executables['blast'], classifFastas[classif])
					makeblastdb_AllByAll_call = [ '{0}/makeblastdb'.format(executables['blast']), '-in', classifFastas[classif], '-dbtype', 'nucl' ]
					append2logfile(paths['output_top_dir'], mainlogfile, 'Began creating blast db for all-by-all blast of LTR_retrotransposon sequences using the call:\n{0}'.format(makeblastdb_AllByAll_call_string))
					makecall(makeblastdb_AllByAll_call, '{0}.makeblastdb.stdout'.format(classifFastas[classif]), '{0}.makeblastdb.stderr'.format(classifFastas[classif]))
					append2logfile(paths['output_top_dir'], mainlogfile, 'Finished creating blast db for all-by-all blast of {0} LTR_retrotransposon sequences.'.format(classif))
					# Perform all-by-all blastn of true positive LTR_retrotransposon sequences
					paths['LTR_retrotransposonAllByAllblastnTable'] = '{0}/{1}.LTRdigest_LTR_retrotransposon_{2}_AllByAll.blastn.tab'.format(outputPth, filenames['inputFasta'], classif)
					blastn_AllByAll_call = [ '{0}/blastn'.format(executables['blast']), '-db', classifFastas[classif], '-query',  classifFastas[classif], '-evalue', '1e-5', '-outfmt', '6', '-num_threads', str(procs), '-max_hsps', '25' ]
					blastn_AllByAll_call_string = '{0} -db {1} -query {1} -evalue 1e-5 -outfmt 6 -num_threads {2} -max_hsps 25 1>{3} 2>{3}.stderr'.format('{0}/blastn'.format(executables['blast']),  classifFastas[classif], procs, paths['LTR_retrotransposonAllByAllblastnTable'])
					append2logfile(paths['output_top_dir'], mainlogfile, 'Began blastn all-by-all blast of LTR_retrotransposon sequences using the call:\n{0}'.format(blastn_AllByAll_call_string))
					makecall(blastn_AllByAll_call, paths['LTR_retrotransposonAllByAllblastnTable'], '{0}.stderr'.format(paths['LTR_retrotransposonAllByAllblastnTable']))
					append2logfile(paths['output_top_dir'], mainlogfile, 'Finished blastn all-by-all blast of {0} LTR_retrotransposon sequences.'.format(classif))

					# Perform MCL clustering based on all-by-all blast results convert blastn output  to abc format for mcl
					paths['MCL_{0}_abc'.format(classif)] = '{0}/LTR_retrotransposon.AllByAllblastn.abc'.format(outputPth)
					blasttable2abc_call = ['cut', '-f', '1,2,11', paths['LTR_retrotransposonAllByAllblastnTable']]
					blasttable2abc_call_string = 'cut -f 1,2,11 {0} > {1} 2>{1}.stderr'.format(paths['LTR_retrotransposonAllByAllblastnTable'], paths['MCL_{0}_abc'.format(classif)])
					append2logfile(paths['output_top_dir'], mainlogfile, 'Began converting blastn tabular output to abc format for MCL using the call:\n{0}'.format(blasttable2abc_call_string))
					makecall(blasttable2abc_call, paths['MCL_{0}_abc'.format(classif)], '{0}.stderr'.format(paths['MCL_{0}_abc'.format(classif)]))
					append2logfile(paths['output_top_dir'], mainlogfile, 'Finished converting blastn tabular output to abc format for MCL')

					with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
						statusFlAppend.write('MCL_{0}_abc\t{1}\n'.format(classif, paths['MCL_{0}_abc'.format(classif)]))

				# create network for mcl
				paths['LTRretrotransposonNetworkmci'] = '{0}/LTR_retrotransposon_{1}.mci'.format(paths['MCL_{0}_I{1}_dir'.format(classif, I)], classif)
				paths['LTRretrotransposonNetworktab'] = '{0}/LTR_retrotransposon_{1}.tab'.format(paths['MCL_{0}_I{1}_dir'.format(classif, I)], classif)
				mcxload_call = ['{0}/mcxload'.format(executables['mcl']), '-abc', paths['MCL_{0}_abc'.format(classif)], '--stream-mirror', '--stream-neg-log10', '-stream-tf', 'ceil(200)', '-o', paths['LTRretrotransposonNetworkmci'], '-write-tab', paths['LTRretrotransposonNetworktab']]
				mcxload_call_string = '{0}/mcxload -abc {1} --stream-mirror --stream-neg-log10 -stream-tf ceil(200) -o {2} -write-tab {3} > {4}/mcxload.stdout 2>{4}mcxload.stderr'.format(executables['mcl'], paths['MCL_{0}_abc'.format(classif)], paths['LTRretrotransposonNetworkmci'], paths['LTRretrotransposonNetworktab'], paths['MCL_{0}_I{1}_dir'.format(classif, I)])
				append2logfile(paths['output_top_dir'], mainlogfile, 'Began creating network for {0} using mcxload:\n{1}'.format(classif, mcxload_call_string))
				makecall(mcxload_call, '{0}/mcxload.stdout'.format(paths['MCL_{0}_I{1}_dir'.format(classif, I)]), '{0}/mcxload.stderr'.format(paths['MCL_{0}_I{1}_dir'.format(classif, I)]))
				append2logfile(paths['output_top_dir'], mainlogfile, 'Finished creating network for {0} LTR RTs using mcxload.'.format(classif))
				
				# cluster using MCL
				current_wd = os.getcwd()
				mcl_call = [ '{0}/mcl'.format(executables['mcl']), 'LTR_retrotransposon_{0}.mci'.format(classif), '-I', str(I), '-use-tab', 'LTR_retrotransposon_{0}.tab'.format(classif), '-te', str(procs) ]
				mcl_call_string = '{0}/mcl {1} -I {2} -use-tab {3} -te {4} > mcl.stdout 2>mcl.stderr'.format(executables['mcl'], 'LTR_retrotransposon_{0}.mci'.format(classif), I, 'LTR_retrotransposons_{0}.tab'.format(classif), str(procs))
				os.chdir(current_wd)
				append2logfile(paths['output_top_dir'], mainlogfile, 'Began clustering {0} using mcl:\n{1}'.format(classif, mcl_call_string))
				os.chdir(paths['MCL_{0}_I{1}_dir'.format(classif, I)])
				makecall(mcl_call, 'mcl.stdout', 'mcl.stderr')
				os.chdir(current_wd)
				append2logfile(paths['output_top_dir'], mainlogfile, 'Finished clustering {0} using mcl'.format(classif))
				os.chdir(paths['MCL_{0}_I{1}_dir'.format(classif, I)])
				os.chdir(current_wd)

				for fl in os.listdir(paths['MCL_{0}_I{1}_dir'.format(classif, I)]):
					if fl.startswith('out'):
						paths['MCL_{0}_I{1}'.format(classif, I)] = '{0}/{1}'.format(paths['MCL_{0}_I{1}_dir'.format(classif, I)], fl)
						
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format('MCL_{0}_I{1}'.format(classif, I), paths['MCL_{0}_I{1}'.format(classif, I)])) 


def full2flankgff(inGFFpth, outGFFpth, bpflank):
	'''
	inGFFpth	A LTRharvest/LTRdigest GFF3
	outGFFpth	The new flanking regions from elements in inGFFpth
	bpflank		The length of each flanking region
	'''
	with open(outGFFpth, 'w') as outFl: # overwrites an old file if present
		outFl.write('##gff-version 3\n')

	with open(inGFFpth, 'r') as inGFFfl:
		outLines = set()
		for line in inGFFfl:
			if not line.startswith('#'):
				gffLine = GFF3_line(line)
				if gffLine.type == 'repeat_region':
					SKIP_L = False
					leftFlankGFFline = GFF3_line()
					leftFlankGFFline.seqid = gffLine.seqid
					leftFlankGFFline.source = 'full2flankgff'
					leftFlankGFFline.type = 'LTR_RT_flank'
					if int(gffLine.start) == 1:
						SKIP_L = True
						leftFlankGFFline.start = 1
					else:
						leftFlankGFFline.start = int(gffLine.start)-bpflank
					if int(leftFlankGFFline.start) < 1:
						leftFlankGFFline.start = 1
					leftFlankGFFline.end = int(gffLine.start)-1
					leftFlankGFFline.score = '.'
					leftFlankGFFline.strand = gffLine.strand
					leftFlankGFFline.phase = '.'
					leftFlankGFFline.attributes['ID'] = '{0}_leftflank_{1}_bp'.format(gffLine.attributes['ID'], bpflank)
					leftFlankGFFline.attributes['Parent'] = gffLine.attributes['ID']
					leftFlankGFFline.attributes_order = ['ID', 'Parent']
					leftFlankGFFline.refreshAttrStr()
					if not SKIP_L:
						outLines.add(str(leftFlankGFFline))
					
					rightFlankGFFline = GFF3_line()
					rightFlankGFFline.seqid = gffLine.seqid
					rightFlankGFFline.source = 'full2flankgff'
					rightFlankGFFline.type = 'LTR_RT_flank'
					rightFlankGFFline.start = int(gffLine.end)+1
					rightFlankGFFline.end = int(gffLine.end)+bpflank
					rightFlankGFFline.score = '.'
					rightFlankGFFline.strand = gffLine.strand
					rightFlankGFFline.phase = '.'
					rightFlankGFFline.attributes['ID'] = '{0}_rightflank_{1}_bp'.format(gffLine.attributes['ID'], bpflank)
					rightFlankGFFline.attributes['Parent'] = gffLine.attributes['ID']
					rightFlankGFFline.attributes_order = ['ID', 'Parent']
					rightFlankGFFline.refreshAttrStr()
					outLines.add(str(rightFlankGFFline))

		with open(outGFFpth, 'w') as outFl:
			for line in outLines:
				outFl.write('{0}\n'.format(line)) 


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


def getfasta(inGFFpth, fastaRefPth, outFastaPth, headerKey='ID'):
	'''
	generic function to run bedtools getfasta
	'''
	call = [  executables['bedtools'], 'getfasta', '-fi', fastaRefPth, '-s', '-bed', inGFFpth ]
	call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2}'.format(fastaRefPth, inGFFpth, outFastaPth)
	makecall(call, stdout=outFastaPth)
	ChangeFastaHeaders(outFastaPth, inGFFpth, attribute=headerKey)
	removeRedundant(outFastaPth)


def runblast(query, subject, out, evalue, outfmt, percid, blast='blastn', procs=1):
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
	blast_call = [ '{0}/{1}'.format(executables['blast'], blast), '-db', subject, '-query', query, '-evalue', str(evalue), '-outfmt', str(outfmt), '-num_threads', str(procs), '-perc_identity', str(percid) ]
	blast_call_string = '{0}/{1} -db {2} -query {3} -evalue {4} -outfmt {5} -num_threads {6} -perc_identity {7} 1>{8} 2>{9}'.format(executables['blast'], blast, subject, query, evalue, outfmt, procs, percid, out, '{0}/runblast.{1}.err'.format('/'.join(out.split('/')[:-1]), blast))
	makecall(blast_call, stdout=out, stderr='{0}/runblast.{1}.err'.format('/'.join(out.split('/')[:-1]), blast))


def reportpairswithhomologousflanks(blastqueryfasta, blastResults, outFlPth, bpflank, perc_len_cutoff):
	'''
	will need to count if elements have 1 or 2 flanks
	'''
	flank_counts = {}
	with open(blastqueryfasta, 'r') as fastafl:
		for line in fastafl:
			if line.startswith('>'):
				el = line[1:].split('_')[1].lstrip('region')
				if el in flank_counts:
					flank_counts[el] += 1
				else:
					flank_counts[el] = 1
	pairs = set()
	with open(blastResults, 'r') as blastfl:
		for line in blastfl:
			if not line.startswith('#'):
				contents = line.strip().split('\t')
				aln_len = int(contents[3])
				rpt1 = contents[0].split('_')[1].lstrip('region')
				rpt2 = contents[1].split('_')[1].lstrip('region')
				if rpt1 == rpt2:
					continue
				else:
					if 100*(aln_len/bpflank) >= float(perc_len_cutoff):
						pairs.add(frozenset([contents[0], contents[1]]))

	pairs = sorted([ sorted(list(pair)) for pair in pairs ], key=lambda x:x[0])
	matches = {}
	for pair in pairs:
		rpt1 = pair[0].split('_')[1].lstrip('region')
		rpt2 = pair[1].split('_')[1].lstrip('region')
		flank1 = pair[0].split('_')[2]
		flank2 = pair[1].split('_')[2]
		if flank1 == 'leftflank':
			f1 = 'L1'
			F1 = 'L2'
		elif flank1 == 'rightflank':
			f1 = 'R1'
			F1 = 'R2'
		if flank2 == 'leftflank':
			f2 = 'L2'
			F2 = 'L1'
		elif flank2 == 'rightflank':
			f2 = 'R2'
			F2 = 'R1'

		match1 = frozenset([f1,f2])
		match2 = frozenset([F1,F2])

		if rpt1 in matches:
			if rpt2 in matches[rpt1]:
				matches[rpt1][rpt2].add(match1)
			else:
				matches[rpt1][rpt2] = set([match1])
		else:
			matches[rpt1] = {rpt2:set([match1])}

		if rpt2 in matches:
			if rpt1 in matches[rpt2]:
				matches[rpt2][rpt1].add(match2)
			else:
				matches[rpt2][rpt1] = set([match2])
		else:
			matches[rpt2] = {rpt1:set([match2])}

	options = [ set([frozenset(['L1', 'L2']), frozenset(['R1','R2'])]), set([frozenset(['L1', 'R2']), frozenset(['L2', 'R1'])]) ]

	matchesset =  set()
	for rpt1 in matches:
		for rpt2 in matches[rpt1]:
			if flank_counts[rpt2] == 1:
				matchesset.add(frozenset([rpt1, rpt2]))
			else:
				if matches[rpt1][rpt2] in options:
					matchesset.add(frozenset([rpt1, rpt2]))
					
	if os.path.isfile(outFlPth):
		os.remove(outFlPth)
	for match in matchesset:
		match = list(match)
		with open(outFlPth, 'a') as outFl:
			outFl.write('{0}\t{1}\n'.format('LTR_retrotransposon{0}'.format(match[0]), 'LTR_retrotransposon{0}'.format(match[1])))


def elementsWithHomologousFlanks(ingff, infasta, outdir, bpflank=None, outfmt='7', percid=None, evalue=None, perc_len_cutoff=None, procs=None):
	'''
	Coordinate discovering element pairs with homologous flanking regions using blastn.
	'''
	MakeDir('FlankDir', outdir)
	blastout = '{0}/blast_all_by_all_flanks'.format(outdir)
	ingffBasename = ingff.split('/')[-1].rstrip('.gff')
	flankfasta = '{0}/{1}.flanks.fasta'.format(outdir, ingffBasename)
	flankgff = '{0}/{1}.flanks.gff'.format(outdir, ingffBasename)
	if not os.path.isfile(blastout):
		if os.path.isfile(flankfasta):
			os.remove(flankfasta)
		if os.path.isfile(flankgff):
			os.remove(flankgff)

		full2flankgff(inGFFpth=ingff, outGFFpth=flankgff, bpflank=bpflank)
		getfasta(inGFFpth=flankgff, fastaRefPth=infasta, outFastaPth=flankfasta, headerKey='ID')
		runblast(query=flankfasta, subject=flankfasta, out=blastout, evalue=str(evalue), outfmt=outfmt, percid=percid, blast='blastn', procs=procs)
	
	reportfl = '{0}/LTR_element_pairs_with_homologous_flanks'.format(outdir)
	reportpairswithhomologousflanks(blastqueryfasta=flankfasta, blastResults=blastout, outFlPth=reportfl, bpflank=bpflank, perc_len_cutoff=perc_len_cutoff)


def aligner(elementList, OutDir, statusFlAlnKey, part):
	'''
	Generic aligner to work with internal and whole elements from LTRharvest GFF3. Called by AutoAlign()
	1. Writes GFF for elements requested (elementList)
	2. Extracts sequences using bedtools getfasta and converts the headers to reflect the GFF3 ID attributes.
	3. Aligns using mafft with given settings.
	4. Trims alignment using trimal -automated1
	5. Adds statsFlAlnKey to the status file.
	'''
	global paths

	if part == 'entire':
		ENTIRE = True
		INTERNAL = False
	elif part == 'internal':
		INTERNAL = True
		ENTIRE = True

	dirKey = '{0}-dir'.format(statusFlAlnKey)
	MakeDir(dirKey, OutDir)
	if not checkStatusFl(statusFlAlnKey):

		if INTERNAL:
			paths['InternalelementsGFF'] = '{0}/elements.gff'.format(OutDir)
			paths['InternalelementsFasta'] = '{1}/elements.fasta'.format(OutDir)
			writeLTRretrotransposonInternalRegions(paths['CurrentGFF'], paths['InternalelementsGFF'], elementSet=set(elementList), truncateParent=True)
			getfasta_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', paths['InternalelementsGFF'] ]
			getfasta_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['InternalelementsGFF'], paths['InternalelementsFasta'], '{0}/bedtools_getfasta_InternalRegions.stderr'.format(OutDir))
			makecall(getfasta_call, paths['InternalelementsFasta'], '{0}/bedtools_getfasta_InternalRegions.stderr'.format(OutDir))
			ChangeFastaHeaders(paths['InternalelementsFasta'], paths['InternalelementsGFF'], attribute='Parent')
			paths['AlnFasta'] = paths['InternalelementsFasta']
			paths['AlnPth'] = '{0}.aln'.format(paths['InternalelementsFasta'])

		elif ENTIRE:

			paths['EntireelementsGFF'] = '{0}/elements.gff'.format(OutDir)
			paths['EntireelementsFasta'] = '{0}/elements.fasta'.format(OutDir)
			writeLTRretrotransposonGFF(inputGFFpth=paths['CurrentGFF'], outputGFFpth=paths['EntireelementsGFF'], elementSet=set(elementList))
			getfasta_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', paths['EntireelementsGFF'] ]
			getfasta_call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(paths['inputFasta'], paths['EntireelementsGFF'], paths['EntireelementsFasta'], '{0}/bedtools_getfasta_Entire.stderr'.format(OutDir))
			makecall(getfasta_call, paths['EntireelementsFasta'], '{0}/bedtools_getfasta_Entire.stderr'.format(OutDir))
			ChangeFastaHeaders(paths['EntireelementsFasta'], paths['EntireelementsGFF'], attribute='ID')
			paths['AlnFasta'] = paths['EntireelementsFasta']
			paths['AlnPth'] = '{0}.aln'.format(paths['EntireelementsFasta'])

		# Align regions from selected elements
		if len(elementList) < mafft_large_minclustsize:
			mafft_call = [ executables['mafft'], '--quiet', '--retree', str(mafft_retree), '--thread', str(procs), '--maxiterate', str(mafft_small_maxiterate), paths['AlnFasta'] ]
			mafft_call_string = '{0} {1}'.format(' '.join(mafft_call),' >{0} 2>{0}.stderr'.format(paths['AlnPth']))
		elif len(elementList) >= mafft_large_minclustsize:
			mafft_call = [ executables['mafft'], '--quiet', '--retree', str(mafft_retree),'--thread', str(procs), '--maxiterate', str(mafft_large_maxiterate), paths['AlnFasta'] ]
			mafft_call_string = '{0} {1}'.format(' '.join(mafft_call),' >{0} 2>{0}.stderr'.format(paths['AlnPth']))

		append2logfile(paths['output_top_dir'], mainlogfile, 'Aligning\n{0}'.format(mafft_call_string))
		makecall(mafft_call, paths['AlnPth'], '{0}.stderr'.format(paths['AlnPth']))
		append2logfile(paths['output_top_dir'], mainlogfile, 'Finished aligning')

		# Trim alignment
		paths['ClusterTrimmedAln'] = '{0}.trimal'.format(paths['AlnPth'])
		trimal_Aln_call =  [ executables['trimal'], '-in', paths['AlnPth'], '-out', paths['ClusterTrimmedAln'], '-automated1' ]
		trimal_Aln_call_string = '{0} -in {1} -out {2} -automated1 >{2}.stdout 2>{2}.stderr'.format(executables['trimal'], paths['AlnPth'], paths['ClusterTrimmedAln'])
		append2logfile(paths['output_top_dir'], mainlogfile, 'Began cleaning alignment using TrimAl:\n{0}'.format(trimal_Aln_call_string))
		makecall(trimal_Aln_call, '{0}.stdout'.format(paths['ClusterTrimmedAln'], '{0}.stdout'.format(paths['ClusterTrimmedAln'])))
		append2logfile(paths['output_top_dir'], mainlogfile, 'Finished cleaning alignment using TrimAl')
		
		paths[dirKey] = OutDir
		if os.path.isfile(paths['ClusterTrimmedAln']): # If the cluster has one element only then there will be no alignment and writing to status could cause problems later
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				paths[statusFlAlnKey] = paths['ClusterTrimmedAln']
				if os.path.isfile(paths[statusFlAlnKey]):
					statusFlAppend.write('{0}\t{1}\n'.format(statusFlAlnKey, paths['ClusterTrimmedAln']))


def AutoAlign(I=6, part='entire', rmgeneconv=False, minClustSize=4, align='clusters', rmhomologflank=False, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, auto_outgroup=False, bpflank=None, flank_pId=None, flank_evalue=None, flank_plencutoff=None, combine_and_do_small_clusters=True):
	'''
	AutoAlign handles aligning whole or internal regions of LTRharvest GFF3 LTR RTs that had been clustered
	using MCL() or WickerFam().

	Options:
	-------------
	I				4		inflation parameter used for MCL
	part				'entire'	'entire' or 'internal', referring to the part of the LTR RTs to align.
	rmgeneconv			False		Remove elements within cluster that show evidence of inter-element gene conversion.
	minClustSize			4		The minimum cluster size to align. If combine_and_do_small_clusters=True, clusters smaller than minClustSize will be aligned together.
	align				'clusters'	'clusters' (align clusters separately) or 'classif' (align all clusters for a given superfamily together)
	rmhomologflank			False		Identify and remove one of every pair of elements from a given cluster whose flanking regions are homologous.
	clustering_method		'WickerFam'	'WickerFam' or 'MCL'. Clustering needs to have been done already.
	WickerParams					Parameters used with the WickerFam() clustering to be aligned.
	auto_outgroup			False		Automatically pick an ougroup and align with cluster. This option is not available for align='classif'
	bpflank				None		If rmhomologflank=True, this specifies how many bp beyond the edge of each element to test for homology.
	flank_pId			None		Minimum percent identity in flank alignment to consider as evidence of homology.
	flank_evalue			None		Minimum E-value of flank alignment to consider as evidence of homology.
	flank_plencutoff		None		Minimum percentage of flanking region required to participate in alignment to consider as evidence of homology.
	combine_and_do_small_clusters	True		Align all elements in clusters smaller than minClustSize together.
	'''
	global paths

	ENTIRE = False
	INTERNAL = False
	AUTO_OUTGROUP = auto_outgroup
	if align == 'clusters':
		ALIGNCLUSTERS = True
		ALIGNCLASSIFS = False
	elif align == 'classifs':
		ALIGNCLASSIFS = True
		ALIGNCLUSTERS = False
	else:
		sys.exit("AutoAlign() parameter 'align' needs to be either 'clusters' or 'classifs', and it is {0}".format(align))
	if rmhomologflank:
		REMOVEHOMOLOGOUSFLANK = True
	else:
		REMOVEHOMOLOGOUSFLANK = False
	if rmgeneconv:
		REMOVEGENECONV = True
		gc = 'GCfiltered'
	else:
		REMOVEGENECONV = False
		gc = 'NoGCfiltering'
	if part == 'entire':
		ENTIRE = True
	elif part == 'internal':
		INTERNAL = True
	else:
		sys.exit('part parameter to AutoAlign() is invalid: {0}. Must be either "entire" or "internal"'.format(part))

	WICKERCLUST=False
	MCLCLUST=False
	if clustering_method == 'WickerFam':
		clusterDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
		MakeDir('WickerAlignments', '{0}/Alignments'.format(clusterDir))
		paths['Alignments'] = paths['WickerAlignments']
		WICKERCLUST=True
	elif clustering_method == 'MCL':
		clusterDir = paths['MCL_I{0}'.format(I)]
		MakeDir('MCLAlignments', '{0}/Alignments'.format(clusterDir))
		paths['Alignments'] = paths['MCLAlignments']
		MCLCLUST=True
	else:
		sys.exit('AutoAlign() parameter clustering_method needs to be either WickerFam or MCL, and it is {0}'.format(clustering_method))

	if INTERNAL:
		MakeDir('InternalRegionsAlignments', '{0}/InternalRegions'.format(paths['Alignments']))
		paths['AlnDir'] = paths['InternalRegionsAlignments']
	elif ENTIRE:
		MakeDir('WholeElementAlignments', '{0}/WholeElements'.format(paths['Alignments']))
		paths['AlnDir'] = paths['WholeElementAlignments']

	if ALIGNCLASSIFS:
		MakeDir('Superfamilies', '{0}/Superfamilies'.format(paths['AlnDir']))
		if REMOVEHOMOLOGOUSFLANK:
			strHomoflank = 'homoflank'
			MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(paths['Superfamilies']))
		else:
			strHomoflank = 'nohomoflank'
			MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['Superfamilies']))
		paths['AlnDir'] = paths['HomoFlankDir']

	elif ALIGNCLUSTERS:
		if REMOVEGENECONV:
			MakeDir('GCDir', '{0}/GeneconversionDisallowed'.format(paths['AlnDir']))
		else:
			MakeDir('GCDir', '{0}/NoGCFiltering'.format(paths['AlnDir']))
		if REMOVEHOMOLOGOUSFLANK:
			strHomoflank = 'homoflank'
			MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(paths['GCDir']))
		else:
			strHomoflank = 'nohomoflank'
			MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['GCDir']))
		if AUTO_OUTGROUP:
			strOutgroup = 'withOutgroup'
			MakeDir('OutgroupDir', '{0}/WithOutgroup'.format(paths['HomoFlankDir']))
		else:
			strOutgroup = 'noOutgroup'
			MakeDir('OutgroupDir', '{0}/NoOutgroup'.format(paths['HomoFlankDir']))
		paths['AlnDir'] = paths['OutgroupDir']

	gcDct = {}
	if REMOVEGENECONV:
		if MCLCLUST:
			gcSummaryPth = 'MCL_I{0}_GENECONV_summary'.format(I)
		elif WICKERCLUST:
			gcSummaryPth = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		if os.path.isfile(paths[gcSummaryPth]):
			with open(paths[gcSummaryPth], 'r') as gcFl:
				for line in gcFl:
					el, clust, classif, g = line.strip().split()
					el = 'LTR_retrotransposon{0}'.format(el)
					clust = int(clust)
					if classif in gcDct:
						if clust in gcDct[classif]:
							if not el in gcDct[classif][clust]:
								gcDct[classif][clust].append(el)
						else:
							gcDct[classif][clust] = [el]
					else:
						gcDct[classif] = {clust:[el]}
		else:
			print('modeltest(removegeneconv=True) used but {0} not found'.format(paths[gcSummaryPth]), file=sys.stderr)

	for classif in classifs:
		if clustering_method == 'WickerFam':
			clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]
		elif clustering_method == 'MCL':
			clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
		smallClusters = []
		MakeDir('AlnDir_{0}'.format(classif), '{0}/{1}'.format(paths['AlnDir'], classif))
		clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]
		if ALIGNCLASSIFS:
			if WICKERCLUST:
				statusFlKey = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_all.{4}.{5}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, strHomoflank)
			if MCLCLUST:
				statusFlKey = 'Aln_{0}_I{1}_cluster_all.{2}.{3}'.format(classif, I, strHomoflank)
			elementlist = [ el for clust in clusters for el in clust ]
			# No gene conversion implemented for this step
			if REMOVEHOMOLOGOUSFLANK:
				outDir = paths['AlnDir_{0}'.format(classif)]
				repeat_region_gff = '{0}/repeat_regions.gff'.format(outDir)
				elements = set()
				for el in elementlist:
					el = el.lstrip('LTR_retrotransposon')
					el = 'repeat_region{0}'.format(el)
					elements.add(el)
				with open(paths['CurrentGFF'], 'r') as inFl:
					with open(repeat_region_gff, 'w') as outFl:
						for line in inFl:
							if '\trepeat_region\t' in line:
								gffLine = GFF3_line(line)
								if gffLine.attributes['ID'] in elements:
									outFl.write(line)

				elementsWithHomologousFlanks(ingff=repeat_region_gff, infasta=paths['inputFasta'], outdir=outDir, bpflank=bpflank, outfmt='7', percid=flank_pId, evalue=flank_evalue, perc_len_cutoff=flank_plencutoff, procs=procs)

				if 'LTR_element_pairs_with_homologous_flanks' in os.listdir(outDir):
					print('Found some homologous flankings', file=sys.stderr)
					print('{0}/LTR_element_pairs_with_homologous_flanks'.format(outDir), file=sys.stderr)
					# parse file and choose LTR RTs to leave out
					with open('{0}/LTR_element_pairs_with_homologous_flanks'.format(outDir), 'r') as pairsFl:
						to_remove = [ line.strip().split('\t')[0] for line in pairsFl ]
						write2summary('These elments were found have homologous flanking sequences with another {0} element and were therefore removed:\n{1}\n'.format(classif, '\n'.join(to_remove)))
						elementlist = [ el for el in elementlist if not el in to_remove ]

				else:
					print('Found no homologous flankings', file=sys.stderr)
					print(outDir, file=sys.stderr)

				aligner(elementList=elementlist, OutDir=outDir, statusFlAlnKey=statusFlKey.format(classif), part=part)

		elif ALIGNCLUSTERS:
			if AUTO_OUTGROUP:
				# Write header for record of outgroups file, overwriting an old file if it exisits.
				if clustering_method == 'WickerFam':
					OutgroupSummaryKey = 'WickerOutgroups_{0}_pId_{1}_percAln_{2}_minLen_{3}_{4}.{5}.{6}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, gc, strHomoflank, strOutgroup)
				elif clustering_method == 'MCL':
					OutgroupSummaryKey = 'MCLOutgroups_{0}_I{1}_{2}.{3}.{4}'.format(classif, I, gc, strHomoflank, strOutgroup)
				outgroupFile = '{0}/Outgroups'.format(paths['AlnDir_{0}'.format(classif)])
				paths[OutgroupSummaryKey] = outgroupFile

				with open(outgroupFile, 'w') as outFl:
					outFl.write('Alignment\toutgroup\toutgroup_cluster\n')

			for j in range(len(clusters)):
				elementlist = clusters[j]
				if REMOVEGENECONV:
					if os.path.isfile(paths[gcSummaryPth]):
						try:
							if j in gcDct[classif]:
								elementlist = [ el for el in clusters[j] if not el in gcDct[classif][j] ]
						except KeyError:
							pass
					else:
						pass # No gene conversion was found at all
				else:
					pass

				if len(elementlist) < minClustSize: # Put small clusters together and align at the end
					if combine_and_do_small_clusters:
						smallClusters += elementlist
				else:
					MakeDir('AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, I, j, gc), '{0}/cluster_{1}'.format(paths['AlnDir_{0}'.format(classif)], j))
					outDir = paths['AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, I, j, gc)]
					if REMOVEHOMOLOGOUSFLANK:
						repeat_region_gff = '{0}/repeat_regions.gff'.format(outDir)
						elements = set()
						for el in elementlist:
							el = el.lstrip('LTR_retrotransposon')
							el = 'repeat_region{0}'.format(el)
							elements.add(el)
						with open(paths['CurrentGFF'], 'r') as inFl:
							with open(repeat_region_gff, 'w') as outFl:
								for line in inFl:
									if '\trepeat_region\t' in line:
										gffLine = GFF3_line(line)
										if gffLine.attributes['ID'] in elements:
											outFl.write(line)

						elementsWithHomologousFlanks(ingff=repeat_region_gff, infasta=paths['inputFasta'], outdir=outDir, bpflank=bpflank, outfmt='7', percid=flank_pId, evalue=flank_evalue, perc_len_cutoff=flank_plencutoff, procs=procs)
						if 'LTR_element_pairs_with_homologous_flanks' in os.listdir(outDir):
							write2summary('Found some homologous flankings among cluster {0}, clustering {1}'.format(j, clustering_method))
							# parse file and choose LTR RTs to leave out
							with open('{0}/LTR_element_pairs_with_homologous_flanks'.format(outDir), 'r') as pairsFl:
								to_remove = [ line.strip().split('\t')[0] for line in pairsFl ]
								write2summary('These elments were found have homologous flanking sequences with another {0} element and were therefore removed:\n{1}\n'.format(classif, '\n'.join(to_remove)))
								elementlist = [ el for el in elementlist if not el in to_remove ]
						else:
							write2summary('Found no homologous flankings among cluster {0}, clustering {1}'.format(j, clustering_method))
					if clustering_method == 'WickerFam':
						statusFlKey = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_{4}_{5}.{6}.{7}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, j, gc, strHomoflank, strOutgroup)
					elif clustering_method == 'MCL':
						statusFlKey = 'Aln_{0}_I{1}_cluster{2}_{3}.{4}.{5}'.format(classif, I, j, gc, strHomoflank, strOutgroup)
					if AUTO_OUTGROUP:
						# the outgroup shall be a random element from cluster k
						# where cluster k is the largest of the clusters that is not j
						# if j is the first cluster then the next smallest cluster is k
						# if there is no other cluster, no outgroup is used
						OUTGROUP_POSSIBLE = True
						if j==0:
							try:
								clusters[j+1]
								k = j+1
							except KeyError:
								OUTGROUP_POSSIBLE = False
						elif j > 0:
							k = 0

						if OUTGROUP_POSSIBLE:
							outgroup = random.choice(clusters[k])
							clustAndOutgroup = clusters[j] + [outgroup]
							with open(outgroupFile, 'a') as outFl:
								outFl.write('{0}\t{1}\t{2}\n'.format(statusFlKey, outgroup, k))
							elementlist = clustAndOutgroup
					if not statusFlKey in paths:
						aligner(elementlist, OutDir=outDir, statusFlAlnKey=statusFlKey, part=part)

			if combine_and_do_small_clusters:

				if len(smallClusters) > 0:
					MakeDir('AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, I, 'small', gc), '{0}/cluster_{1}'.format(paths['AlnDir_{0}'.format(classif)], 'small'))
					outDir = paths['AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, I, 'small', gc)]
					if clustering_method == 'WickerFam':
						statusFlKey = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_{4}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, gc)
					elif clustering_method == 'MCL':
						statusFlKey = 'Aln_{0}_I{1}_clustersmall_{2}'.format(classif, I, gc)
					if not statusFlKey in paths:
						aligner(smallClusters, OutDir=outDir, statusFlAlnKey=statusFlKey, part=part)


def geneconvClusters(trimal=True, g='/g0', force=False, clust=None, I=6, minClustSize=4, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, combine_and_do_small_clusters=True):
	'''
	g can be one of /g0, /g1, or /g2
	g is proportional to the tolerance for mismatches in fragments by geneconv
	'''
	# GENECONV output:
	##   Names                                                                                                    Pvalue  Pvalue   Begin  End   Len  Poly Dif  Difs Pen.
	#GI      LTR_retrotransposon4189;LTR_retrotransposon4189 0.0001  0.01017 26      100     75      75      0       10      None
	##
	global paths
	global filenames
	TRIMAL = trimal
	
	if GENECONVCLUSTERS:

		WICKERCLUST = False
		MCLCLUST = False
		if clustering_method == 'WickerFam':
			WICKERCLUST = True
		elif clustering_method == 'MCL':
			MCLCLUST = True
		else:
			sys.exit('geneconvClusters() parameter clustering_method needs to be either WickerFam or MCL, and it is: {0}'.format(clustering_method))

		if WICKERCLUST:

			WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
			paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])] = '{0}/GENECONV'.format(WickerDir)
			WickerGCdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
			if not checkStatusFl(WickerGCdirkey):
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format(WickerGCdirkey, paths[WickerGCdirkey]))
			geneconvOutputDir = WickerGCdirkey
			gcSummaryFl = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_{3}_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], g[1:])

		elif MCLCLUST:

			MCLdir = paths['MCL_I{0}'.format(I)]
			paths['MCL_I{0}_GENECONVdir'.format(I)] = '{0}/GENECONV'.format(MCLdir)
			if not checkStatusFl('MCL_I{0}_GENECONVdir'.format(I)):
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format('MCL_I{0}_GENECONVdir'.format(I), paths['MCL_I{0}_GENECONVdir'.format(I)]))
			geneconvOutputDir = 'MCL_I{0}_GENECONVdir'.format(I)
			gcSummaryFl = 'MCL_I{0}_GENECONV_{1}_summary'.format(I, g[1:])

		if not gcSummaryFl in paths:

			geneconv_calls = []
			append2logfile(paths['output_top_dir'], mainlogfile, 'Checking directory structure for GENECONV using {0}'.format(g) )

			for classif in classifs:
				MakeDir('GENECONV_{0}_dir'.format(classif), '{0}/{1}'.format(paths[geneconvOutputDir], classif))
				MakeDir('GENECONV_{0}_{1}_dir'.format(classif, g[1:]), '{0}/{1}'.format(paths['GENECONV_{0}_dir'.format(classif)], g[1:]))

				if MCLCLUST:
					clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
				elif WICKERCLUST:
					clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]

				clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]
				smalls = []
				for j in range(len(clusters)):
					if len(clusters[j]) >= minClustSize:
						if MCLCLUST:
							alnPth = paths['Aln_{0}_I{1}_cluster{2}_NoGCfiltering.nohomoflank.noOutgroup'.format(classif, I, j)]
						elif WICKERCLUST:
							alnPth = paths['WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_{4}_NoGCfiltering.nohomoflank.noOutgroup'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, j)]
						if os.path.isfile(alnPth):
							if not os.stat(alnPth).st_size == 0: # if alignment file is non-empty
								append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing to run GENECONV on:\n{0}'.format(alnPth))
								call = [ executables['geneconv'], alnPth, '/w124', g, '-include_monosites', '-nolog', '-Dumptab', '-Fancy' ]
								geneconv_calls.append((call, '/dev/null', None, None))
						else:
							continue
					else:
						if combine_and_do_small_clusters:
							smalls += clusters[j]

				if combine_and_do_small_clusters:
					if smalls != []:
						if MCLCLUST:
							if 'Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(classif, I) in paths:
								alnPth = paths['Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(classif, I)]
							else:
								continue
						elif WICKERCLUST:
							if 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_NoGCfiltering'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif) in paths:
								alnPth = paths['WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_NoGCfiltering'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]
							else:
								continue
						if os.path.isfile(alnPth):
							if not os.stat(alnPth).st_size == 0: # if alignment file is non-empty
								append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing to run GENECONV on:\n{0}'.format(alnPth))
								call = [ executables['geneconv'], alnPth, '/w124', g, '-include_monosites', '-nolog', '-Dumptab', '-Fancy' ]
								geneconv_calls.append((call, '/dev/null', '{0}.geneconv.err'.format(alnPth), None))
						else:
							continue
			if geneconv_calls == []:
				return
			append2logfile(paths['output_top_dir'], mainlogfile, 'Running GENECONV, example call:\n{0}'.format(' '.join(geneconv_calls[0][0] )))
			chunk_size = ceil(len(geneconv_calls)/procs)
			with Pool(processes=procs) as p:
				p.map(makecallMultiprocessing, geneconv_calls, chunksize=chunk_size)
			p.join()

			hasEvidenceOfGC = set()
			for classif in classifs:
				if MCLCLUST:
					clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
				elif WICKERCLUST:
					clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]

				clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]
				for j in range(len(clusters)):
					if len(clusters[j]) >= minClustSize: 
						if MCLCLUST:
							trimalOutput = paths['Aln_{0}_I{1}_cluster{2}_NoGCfiltering.nohomoflank.noOutgroup'.format(classif, I, j)].split('.')
						elif WICKERCLUST:
							trimalOutput = paths['WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_{4}_NoGCfiltering.nohomoflank.noOutgroup'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, j)].split('.')
						trimalOutput[-1] = 'tab'
						geneconvOutputPth = '.'.join(trimalOutput)

						with open(geneconvOutputPth, 'r') as gcFl:
							for line in gcFl:
								if line.startswith('GI'):
									## Add parsing of new format here and make Circos plots.
								#	#   Seq       Sim     BC KA    Aligned Offsets         In Seq1            In Seq2        Num  Num  Tot  MisM
								#	#   Names    Pvalue   Pvalue   Begin  End   Len    Begin  End   Len   Begin  End   Len   Poly Dif  Difs Pen.
								#	GI  S18;S37  0.0000  4.73e-41   235   1605 1371     235   1605 1371    229   1599 1371    213   0  229  None
								#	GI  S21;S37  0.0000  3.32e-38   235   1383 1149     235   1383 1149    229   1377 1149    179   0  249  None
								#	GI  S34;S37  0.0000  6.18e-37   322   1605 1284     322   1605 1284    316   1599 1284    206   0  218  None
								#	GI  S26;S37  0.0000  1.33e-33   426   1383  958     426   1383  958    420   1377  958    157   0  252  None
									GI, els, sim_p_val, bc_ka_p_val, aln_start, aln_end, aln_len, el1_start, el1_end, el1_len, el2_start, el2_end, el2_len, num_polys, num_difs, tot_difs, mism_pen = line.strip().split('\t')
									element1, element2 = [ i[1:] for i in els.split(';') ]
								#	totDifs = int(line.strip().split()[9])
									if int(tot_difs) < 3:
										continue
									#element1 = line.strip().split()[1].split(';')[0]
									#element2 = line.strip().split()[1].split(';')[1]
									hasEvidenceOfGC.add((element1, j, classif, g[1:]))
									hasEvidenceOfGC.add((element2, j, classif, g[1:]))
									line='{0}\t{1}\t{2}\n'.format(line.strip(), classif, j)
									with open('{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, g[1:]), 'a') as outFl:
										outFl.write(line)

						fName = geneconvOutputPth.split('/')[-1]
						os.rename(geneconvOutputPth, '{0}/clust{1}_{2}'.format(paths['GENECONV_{0}_{1}_dir'.format(classif, g[1:])], j, fName))

				trimalOutput = None
				if combine_and_do_small_clusters:
					if MCLCLUST:
						if 'Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(classif, I) in paths:
							if os.path.isfile(paths['Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(classif, I)]):
								trimalOutput = paths['Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(classif, I)].split('.')
					elif WICKERCLUST:
						if 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_NoGCfiltering'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif) in paths:
							if os.path.isfile(paths['WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_NoGCfiltering'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]): # 1 sequence, no alignment
								trimalOutput = paths['WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_NoGCfiltering'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)].split('.')

				if not trimalOutput == None:

					trimalOutput[-1] = 'tab'
					geneconvOutputPth = '.'.join(trimalOutput)

					with open(geneconvOutputPth, 'r') as gcFl:
						for line in gcFl:
							if line.startswith('GI'):
								totDifs = int(line.strip().split()[9])
								if totDifs < 5:
									continue
								element1 = line.strip().split()[1].split(';')[0]
								element2 = line.strip().split()[1].split(';')[1]
								hasEvidenceOfGC.add((element1.lstrip('S'), j, classif, g[1:]))
								hasEvidenceOfGC.add((element2.lstrip('S'), j, classif, g[1:]))
								line='{0}\t{1}\t{2}\n'.format(line.strip(), classif, j)
								with open('{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, g[1:]), 'a') as outFl:
									outFl.write(line)

					fName = geneconvOutputPth.split('/')[-1]
					os.rename(geneconvOutputPth, '{0}/clust{1}_{2}'.format(paths['GENECONV_{0}_{1}_dir'.format(classif, g[1:])], j, fName))

			for el in sorted(list(hasEvidenceOfGC), key=lambda x:(x[2], x[1], x[3])):
				with open('{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir]), 'a') as outFl:
					outFl.write('{0}\t{1}\t{2}\t{3}\n'.format(el[0], el[1], el[2], g[1:]))
			paths['GENECONV_clusters_I{0}_summary'.format(I)] = '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])

			if MCLCLUST:
				paths['MCL_I{0}_GENECONV_summary'.format(I)] = '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])
			elif WICKERCLUST:
				paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])] = '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])

			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				if MCLCLUST:
					if not checkStatusFl('MCL_I{0}_GENECONV_{1}_summary'.format(I, g[1:])):
						statusFlAppend.write('{0}\t{1}\n'.format('MCL_I{0}_GENECONV_{1}_summary'.format(I, g[1:]), '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])))
					if not checkStatusFl('MCL_I{0}_GENECONV_summary'.format(I)):
						statusFlAppend.write('{0}\t{1}\n'.format('MCL_I{0}_GENECONV_summary'.format(I), '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])))
				elif WICKERCLUST:
					if not checkStatusFl('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_{3}_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], g[1:])):
						statusFlAppend.write('{0}\t{1}\n'.format('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_{3}_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], g[1:]), '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])))
					if not checkStatusFl('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])):
						statusFlAppend.write('{0}\t{1}\n'.format('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen']), '{0}/ElementsWithEvidenceOfGeneConversion'.format(paths[geneconvOutputDir])))
		
def modeltest(iters=1, I=6, removegeneconv=True, part='entire', clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, minClustSize=4, bpflank=None, combine_and_do_small_clusters=True):

	'''
	iters	number of times to run jModeltest2. default 10
	part refers to what portion of the element to align. possible values are:
	entire or internal
	'''
	global paths
	global fileames

	REMOVEGENECONV = removegeneconv
	WICKERCLUST = False
	MCLCLUST = False
	if clustering_method == 'WickerFam':
		WICKERCLUST = True
	elif clustering_method == 'MCL':
		MCLCLUST = True
	else:
		sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL, and it is: {0}'.format(clustering_method))
	
	if REMOVEGENECONV:
		gc = 'GCfiltered'
	else:
		gc = 'NoGCfiltering'
	keySuffix = '{0}.nohomoflank.noOutgroup'.format(gc)

	if LTRDIVERGENCE:
		if WICKERCLUST:
			WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
			WickerMTdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_modeltestDir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
			paths[WickerMTdirkey] = '{0}/Modeltest'.format(WickerDir)
			paths['ModelTestDir'] = paths[WickerMTdirkey]
			if not checkStatusFl(WickerMTdirkey):
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format(WickerMTdirkey, paths[WickerMTdirkey]))
			if REMOVEGENECONV:
				AutoAlign(I=None, part='entire', rmgeneconv=True, minClustSize=minClustSize, align='clusters', rmhomologflank=False, clustering_method='WickerFam', WickerParams={'pId':WickerParams['pId'], 'percAln':WickerParams['percAln'], 'minLen':WickerParams['minLen']}, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
				MakeDir('ModelTestDir_{0}'.format(keySuffix), '{0}/GeneconversionDisallowed'.format(paths['ModelTestDir']))
			else:
				AutoAlign(I=None, part='entire', rmgeneconv=False, minClustSize=minClustSize, align='clusters', rmhomologflank=False, clustering_method='WickerFam', WickerParams={'pId':WickerParams['pId'], 'percAln':WickerParams['percAln'], 'minLen':WickerParams['minLen']}, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
				MakeDir('ModelTestDir_{0}'.format(keySuffix), '{0}/NoGCFiltering'.format(paths['ModelTestDir']))

			alignmentsForModeltesting = [ pth for pth in paths if pth.startswith('WickerAln_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])) and pth.endswith('{0}.nohomoflank.noOutgroup'.format(gc)) ] # model testing not done for small clusters

		elif MCLCLUST:
			MCLdir = paths['MCL_I{0}'.format(I)]
			MCL_MT_dirkey = 'MCL_I{0}_modeltestDir'.format(I)
			paths[MCL_MT_dirkey] = '{0}/Modeltest'.format(MCLdir)
			paths['ModelTestDir'] = paths[MCL_MT_dirkey]
			if not checkStatusFl(MCL_MT_dirkey):
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format(MCL_MT_dirkey, paths[MCL_MT_dirkey]))
			if REMOVEGENECONV:
				AutoAlign(I=I, part='entire', rmgeneconv=True, minClustSize=minClustSize, align='clusters', rmhomologflank=False, clustering_method='MCL', WickerParams=None, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
				MakeDir('ModelTestDir_{0}'.format(keySuffix), '{0}/GeneconversionDisallowed'.format(paths['ModelTestDir']))
			else:
				AutoAlign(I=I, part='entire', rmgeneconv=False, minClustSize=minClustSize, align='clusters', rmhomologflank=False, clustering_method='MCL', WickerParams=None, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
				MakeDir('ModelTestDir_{0}'.format(keySuffix), '{0}/NoGCFiltering'.format(paths['ModelTestDir']))
			alignmentsForModeltesting = [ pth for pth in paths if pth.startswith('Aln_') and pth.endswith('{0}.nohomoflank.noOutgroup'.format(gc)) and 'I{0}'.format(I) in pth ]

		if alignmentsForModeltesting == []:
			with open('{0}/NO_MODEL_TESTING_DONE'.format(paths['ModelTestDir_{0}'.format(keySuffix)]), 'w') as outFl:
				outFl.write('Perhaps there are too few LTR RTs')
			
		for aln in alignmentsForModeltesting:
			if not os.path.isfile(paths[aln]):
				continue
			if MCLCLUST:
				# e.g. aln = Aln_Other_I8_cluster0_GCfiltered.nohomoflank.noOutgroup
				a = aln.split('_')
				classif = a[1]
				if '-' in a[2]:
					j = a[2].split('-')[1][7:]
					suffix = a[3] 

				else:
					j = a[3][7:]
					suffix = a[4]
				OutDirKey = 'MCLModelTestDir_{0}_iters_I{1}_{2}_cluster_{3}_{4}'.format(iters,I, classif, j, suffix)
				OutSummaryKey = 'MCLModelTestSummary_{0}_iters_I{1}_{2}'.format(iters,I, classif)
				if j == 'small':
					if not combine_and_do_small_clusters:	
						continue
					j = 'clustersmall'
			elif WICKERCLUST:
				# e.g. aln = WickerAln_80_pId_80_percAln_80_minLen_Copia_cluster_0_GCfiltered
				a = aln.split('_')
				classif = a[7]
				j = a[-2]
				suffix = a[-1]
				OutDirKey = 'WickerModelTestDir_{0}_iters_{1}_pId_{2}_percAln_{3}_minLen_{4}_cluster_{5}_{6}'.format(iters, WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, j, suffix)
				OutSummaryKey = 'WickerModelTestSummary_{0}_iters_{1}_pId_{2}_percAln_{3}_minLen_{4}'.format(iters, WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)
					
			MakeDir('ModelTest_{0}_{1}_dir'.format(suffix, classif), '{0}/{1}'.format(paths['ModelTestDir_{0}'.format(suffix)], classif))
			sessionDir = paths['ModelTest_{0}_{1}_dir'.format(suffix, classif)]

			if not OutDirKey in paths: # If true it means the model testing already finished for that group

				MakeDir(OutDirKey, '{0}/{1}_iters'.format(sessionDir, iters))
				MakeDir('ModelTestClustDir', '{0}/cluster_{1}'.format(paths[OutDirKey], j))

				for i in range(int(iters)):

					MakeDir('ModelTestIterationDir', '{0}/iter_{1}'.format(paths['ModelTestClustDir'], str(i+1)))

					# FastTree
					paths['Tree'] = '{0}/{1}_I{2}_{3}.tree'.format(paths['ModelTestIterationDir'], classif, I, j)
					filenames['Tree'] = '{0}_I{1}_{2}.tree'.format(classif, I, j)
					fasttree_call = [ executables['fasttree'], '-nt', '-gtr' ]
					fasttree_call_string =  '{0} -nt -gtr <{1} >{2} 2>{2}.stderr'.format(executables['fasttree'], paths[aln],paths['Tree'])
					append2logfile(paths['output_top_dir'], mainlogfile, 'Began inferring phylogeny using FastTree:\n{0}'.format(fasttree_call_string))
					makecall(fasttree_call, stdout=paths['Tree'], stderr='{0}.stderr'.format(paths['Tree']), stdin=paths[aln])
					append2logfile(paths['output_top_dir'], mainlogfile, 'Finished inferring phylogeny using FastTree')
					paths['jModeltest2out_{0}'.format(classif)] = '{0}/{1}.jModelTest2.out'.format(paths['ModelTestIterationDir'], filenames['Tree'])
					jmodeltestCallString = 'java -jar {0} -d {1} -w -g 4 -f -BIC -a -u {2} -o {3} -tr {4} -s 11'.format(executables['jmodeltest2'], paths[aln], paths['Tree'], str(procs), paths['jModeltest2out_{0}'.format(classif)])
					append2logfile(paths['output_top_dir'], mainlogfile, 'Starting jModeltest2\n{0}'.format(jmodeltestCallString))
					subprocess.call([ 'java', '-jar', executables['jmodeltest2'], '-d', paths[aln], '-w', '-g', '4', '-f', '-BIC', '-a', '-u', paths['Tree'], '-o', paths['jModeltest2out_{0}'.format(classif)], '-tr', str(procs), '-s', '11'])
									
					topDir = paths[OutDirKey]
					jmt2summaryFlPth = '{0}/Summary.txt'.format(topDir)
					paths['jModeltest2summary_{0}'.format(classif)] = jmt2summaryFlPth
					#with open('{0}/{1}/{2}'.format(topDir, d, f), 'r') as jmt2outFl:
					clustSize = len([ 1 for line in open(paths[aln],'r').read().split('\n') if line.startswith('>') ])

					with open(paths['jModeltest2out_{0}'.format(classif)], 'r') as jmt2outFl:
						END = False
						PAUP = False
						getPAUP = False
						COMPLETE_RUN = False
						paupLines = ''
						for line in jmt2outFl:
							if line.startswith('PAUP* Commands Block:'):
								PAUP = True
							elif line.startswith('::Best Models::'):
								END = True
					
								#method, model, a, c, g, t, kappa, titv, Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma = line.strip().split()
								with open(jmt2summaryFlPth, 'a') as summaryFl:
									summaryFl.write('{0}\t{1}\t{2}\n'.format(line.strip(), j, clustSize))
									summaryFl.write(paupLines)
									COMPLETE_RUN = True
							if line.startswith('END;') and getPAUP:
								paupLines += 'Dset distance=ML;\n'
								paupLines += 'SaveDist format=oneColumn file=dist;\n'
								paupLines += line
								getPAUP = False
							elif line.startswith('[!') and PAUP:
								paupLines += line
								PAUP = False
								getPAUP = True
							elif getPAUP:
								paupLines += line

			if not checkStatusFl(OutDirKey):
				with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
					statusFlAppend.write('{0}\t{1}\n'.format(OutDirKey, paths[OutDirKey]))
					#statusFlAppend.write('{0}\t{1}\n'.format('jModeltest2out_{0}'.format(classif), paths['jModeltest2out_{0}'.format(classif)]))
					if 'jModeltest2summary_{0}'.format(classif) in paths and COMPLETE_RUN:
						if not checkStatusFl(OutSummaryKey):
							statusFlAppend.write('{0}\t{1}\n'.format(OutSummaryKey, paths['jModeltest2summary_{0}'.format(classif)]))
					else:
						pass

def align_ltrs(trimal=True, I=6, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}):

	'''
	Run through GFF3
	Make a GFF3 for every LTR pair
	Create a bedtools getfasta call for every LTR pair GFF3
	Change sequence headers
	Create a mafft call for every LTR pair FASTA
	Execute calls using multiprocessing
	'''
	global paths
	global filenames

	TRIMAL = trimal
	WICKERCLUST = False
	MCLCLUST = False
	# And set up directory structure for output (LTR pairs GFFs, FASTAs, and alignments)
	if clustering_method == 'WickerFam':
		WICKERCLUST = True
		key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		WickerDir = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		MakeDir('Alignments', '{0}/Alignments'.format(paths[WickerDir]))
		AlnKey = '{0}.LTRAlnDir'.format(key_base)
		paths[AlnKey] = '{0}/LTRs'.format(paths['Alignments'])
		MakeDir(AlnKey, paths[AlnKey])
		GFFKey = '{0}.GFFDir'.format(key_base)
		paths[GFFKey] = '{0}/GFFs'.format(paths[WickerDir])
		MakeDir(GFFKey, paths[GFFKey])
		FASTAKey = '{0}.FASTADir'.format(key_base)
		paths[FASTAKey] = '{0}/FASTAs'.format(paths[WickerDir])
		MakeDir(FASTAKey, paths[FASTAKey])
		TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
	elif clustering_method == 'MCL':
		MCLCLUST = True
		key_base = 'MCL_I{0}'.format(I)
		MCLdir = 'MCL_I{0}'.format(I)
		MakeDir('Alignments', '{0}/Alignments'.format(paths[MCLdir]))
		AlnKey = '{0}.LTRAlnDir'.format(key_base)
		paths[AlnKey] = '{0}/LTRs'.format(paths['Alignments'])
		MakeDir(AlnKey, paths[AlnKey])
		GFFKey = '{0}.GFFDir'.format(key_base)
		paths[GFFKey] = '{0}/GFFs'.format(paths[MCLdir])
		MakeDir(GFFKey, paths[GFFKey])
		FASTAKey = '{0}.FASTADir'.format(key_base)
		paths[FASTAKey] = '{0}/FASTAs'.format(paths[MCLdir])
		MakeDir(FASTAKey, paths[FASTAKey])
		TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
	else:
		sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL and it is: {0}'.format(clustering_method))

	if checkStatusFl('{0}.LTR_divergence_complete'.format(key_base)):
		append2logfile(paths['output_top_dir'], mainlogfile, 'ltr_divergence() already completed: {0}'.format(paths['{0}.LTR_divergence_complete'.format(key_base)]))
		return

	if GENECONVLTRS or DIVERGENCE:

		MakeDir('LTRsGFFoutputDir', '{0}/LTRs'.format(paths['GFFOutputDir']))
		MakeDir('LTRsFASTAoutputDir', '{0}/LTRs'.format(paths['FastaOutputDir']))
		MakeDir('Alignments', '{0}/Alignments'.format(paths['output_top_dir']))
		MakeDir('LTRsAlignments_dir', '{0}/LTRs'.format(paths['Alignments']))

		ltrs = {} # this will have the GFF lines for each LTR RT's LTRs
		ltrs_getfasta_calls = {}
		ltrs_changefastaheaders_calls = {}
		ltrs_mafft_calls = {}
		ltrs_trimal_calls = {}
		num_pairs = 0

		append2logfile(paths['output_top_dir'], mainlogfile, 'Parsing LTRs from GFF3:\n{0}'.format(paths['CurrentGFF']))
		with open(paths['CurrentGFF'], 'r') as GFF_fl:
			for line in GFF_fl:
				if '\tlong_terminal_repeat\t' in line:
					gffLine = GFF3_line(line)
					elementName = gffLine.attributes['Parent']

					if elementName in ltrs:
						ltrs[elementName].append(line)
						if len(ltrs[elementName]) == 2:
							LTRsGFFfilepath = '{0}/{1}_LTRs.gff'.format(paths[GFFKey], elementName)
							LTRsFASTAfilepath = '{0}/{1}_LTRs.fasta'.format(paths[FASTAKey], elementName)
							LTRsAlignmentFilepath = '{0}/{1}_LTRs.fasta.aln'.format(paths[AlnKey], elementName)
							LTRsTrimmedAlnFilepath = '{0}.trimmed'.format(LTRsAlignmentFilepath)
							ltrs[elementName].append(LTRsGFFfilepath) # For writing GFF3
							getfasta_LTRs_call = [ executables['bedtools'], 'getfasta', '-fi', paths['inputFasta'], '-s', '-bed', LTRsGFFfilepath ]  
							ltrs_getfasta_calls[elementName] = (getfasta_LTRs_call, LTRsFASTAfilepath, None, None)
							ltrs_changefastaheaders_calls[elementName] = (LTRsFASTAfilepath, LTRsGFFfilepath, 'Parent')
							mafft_LTRs_call = [ executables['mafft'], '--quiet', '--globalpair', '--maxiterate', '1000', LTRsFASTAfilepath ]
							ltrs_mafft_calls[elementName] = (mafft_LTRs_call, LTRsAlignmentFilepath, None, None)

							if TRIMAL:
								trimal_LTRs_call = [ executables['trimal'], '-in', LTRsAlignmentFilepath, '-out', LTRsTrimmedAlnFilepath, '-automated1' ]
								ltrs_trimal_calls[elementName] = (trimal_LTRs_call, LTRsTrimmedAlnFilepath, None, None)

							num_pairs += 1
					else:
						ltrs[elementName] = [line]

		chunk_size = ceil(num_pairs/procs)
		# Write to log file here about chunk size and processors used
		append2logfile(paths['output_top_dir'], mainlogfile, 'For align_ltrs(): procs={0} chunk_size={1}'.format(procs,chunk_size))
		if not checkStatusFl(GFFKey):
			append2logfile(paths['output_top_dir'], mainlogfile, 'Writing GFF3s for each LTR pair:\n{0}'.format(paths[GFFKey]) )
			with Pool(processes=procs) as p: # Write GFF3s for each LTR pair
				p.map(write_ltrs_gff3, ltrs.values(), chunksize=chunk_size)
			p.join()
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(GFFKey, paths[GFFKey])) # Add LTRs GFF path to status file (for resuming later)

		if not checkStatusFl(FASTAKey):
			append2logfile(paths['output_top_dir'], mainlogfile, 'Extracting sequences for LTR pairs:\n{0}'.format(list(ltrs_getfasta_calls.values())[0]))
			with Pool(processes=procs) as p: # Write FASTA for each LTR pair
				p.map(makecallMultiprocessing, ltrs_getfasta_calls.values(), chunksize=chunk_size)
			p.join()
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(FASTAKey, paths[FASTAKey])) # Add LTRs FASTA path to status file (for resuming later)

		NewHeaderKey = '{0}.LTRsFASTAnewheaders'.format(key_base)
		if not checkStatusFl(NewHeaderKey):
			append2logfile(paths['output_top_dir'], mainlogfile, 'Changing bedtools getfasta default headers to LTR RT names:\n{0}'.format(paths[FASTAKey] ))
			with Pool(processes=procs) as p: # Write FASTA for each LTR pair
				p.map(ChangeFastaHeadersMultiprocessing, ltrs_changefastaheaders_calls.values(), chunksize=chunk_size)
			p.join()
			paths[NewHeaderKey] = paths[FASTAKey]
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(NewHeaderKey, paths[FASTAKey])) # Add LTRs FASTA path to status file (for resuming later)

		if not checkStatusFl(AlnKey):
			append2logfile(paths['output_top_dir'], mainlogfile, 'Making MAFFT alignments for each LTR pair:\n{0}'.format(list(ltrs_mafft_calls.values())[0]))
			with Pool(processes=procs) as p: # Do alignment for each LTR pair
				p.map(makecallMultiprocessing, ltrs_mafft_calls.values(), chunksize=chunk_size)
			p.join()
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(AlnKey, paths[AlnKey])) # Add LTRs FASTA path to status file (for resuming later)

		if not checkStatusFl(TrimalKey):
			append2logfile(paths['output_top_dir'], mainlogfile, 'Trimming MAFFT alignment using TrimAl -automated1:\n{0}'.format(list(ltrs_trimal_calls.values())[0]) )
			with Pool(processes=procs) as p: # Do alignment for each LTR pair
				p.map(makecallMultiprocessing, ltrs_trimal_calls.values(), chunksize=chunk_size)
			p.join()
			paths[TrimalKey] = paths[AlnKey]
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(TrimalKey, paths[TrimalKey])) # Add LTRs FASTA path to status file (for resuming later)
				statusFlAppend.write('{0}.LTR_divergence_complete\t{0}'.format(key_base))

def geneconvLTRs(trimal=True, g='/g0', force=False, I=6, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}):
	'''
	g can be one of /g0, /g1, or /g2
	g is proportional to the tolerance for mismatches in fragments by geneconv
	trimal should be true if LTRs were aligned then trimmed with trimal
	'''
	# GENECONV output:
	##   Names                                                                                                    Pvalue  Pvalue   Begin  End   Len  Poly Dif  Difs Pen.
	#GI      LTR_retrotransposon4189;LTR_retrotransposon4189 0.0001  0.01017 26      100     75      75      0       10      None
	##
	global paths
	global filenames

	TRIMAL = trimal
	WICKERCLUST = False
	MCLCLUST = False
	# And set up directory structure for output (LTR pairs GFFs, FASTAs, and alignments)
	if clustering_method == 'WickerFam':
		WICKERCLUST = True
		key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		WickerDir = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])

		AlnKey = '{0}.LTRAlnDir'.format(key_base)
		if not checkStatusFl(AlnKey):
			sys.exit('LTR alignment not in status file: {0}'.format(AlnKey))

		GENECONVTopDirKey = '{0}.GENECONV'.format(key_base)
		paths[GENECONVTopDirKey] = '{0}/GENECONV'.format(paths[WickerDir])
		MakeDir(GENECONVTopDirKey, paths[GENECONVTopDirKey])

		GENECONVDirKey = '{0}.GENECONVLTRs'.format(key_base)
		paths[GENECONVDirKey] = '{0}/LTRs'.format(paths[GENECONVTopDirKey])
		MakeDir(GENECONVDirKey, paths[GENECONVDirKey])

		GENECONVgDirKey = 'GENECONVgDir'
		paths[GENECONVgDirKey] = '{0}/{1}'.format(paths[GENECONVDirKey], g[1:])
		MakeDir(GENECONVgDirKey, paths[GENECONVgDirKey])

		SummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)
		
	elif clustering_method == 'MCL':

		MCLCLUST = True
		key_base = 'MCL_I{0}'.format(I)
		MCLdir = 'MCL_I{0}'.format(I)

		AlnKey = '{0}.LTRAlnDir'.format(key_base)
		if not checkStatusFl(AlnKey):
			sys.exit('LTR alignment not in status file: {0}'.format(AlnKey))

		GENECONVTopDirKey = '{0}.GENECONV'.format(key_base)
		paths[GENECONVTopDirKey] = '{0}/GENECONV'.format(paths[MCLdir])
		MakeDir(GENECONVTopDirKey, paths[GENECONVTopDirKey])

		GENECONVDirKey = '{0}.GENECONVLTRs'.format(key_base)
		paths[GENECONVDirKey] = '{0}/LTRs'.format(paths[GENECONVTopDirKey])
		MakeDir(GENECONVDirKey, paths[GENECONVDirKey])

		GENECONVgDirKey = 'GENECONVgDir'
		paths[GENECONVgDirKey] = '{0}/{1}'.format(paths[GENECONVDirKey], g[1:])
		MakeDir(GENECONVgDirKey, paths[GENECONVgDirKey])

		SummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)

	else:
		sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL and it is: {0}'.format(clustering_method))

	if checkStatusFl(SummaryKey):
		append2logfile(paths['output_top_dir'], mainlogfile, 'ltr_divergence() already completed: {0}'.format(paths['{0}.GENECONVLTRs'.format(key_base)]))
		return

	if not checkStatusFl('{0}.GENECONVLTRs.{1}'.format(key_base, g[1:])):

		geneconv_calls = []
		append2logfile(paths['output_top_dir'], mainlogfile, 'Checking directory structure for GENECONV using {0}'.format(g) )

		if TRIMAL:
			files = [ f for f in os.listdir(paths[AlnKey]) if f.endswith('trimmed') ]
		else:
			files = [ f for f in os.listdir(paths[AlnKey]) if f.endswith('aln') ]

		num_elements = len(files)
		alnLens = {}
		show = False
		append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing to run GENECONV for finding intraelement gene conversion between LTRs')
		for f in files:
			flpth = '{0}/{1}'.format(paths[AlnKey], f)
			call = [ executables['geneconv'], flpth, '/w123', g, '-include_monosites', '-nolog', '-Dumptab' ]
			geneconv_calls.append((call, '/dev/null', None, None))
			elementName = '_'.join(f.split('_')[:2])
			seqs = list(SeqIO.parse(flpth, 'fasta'))

			try:
				alnLen = len(seqs[0].seq)
			except IndexError:
				continue

			elementName = seqs[0].id
			alnLens[elementName] = alnLen
			show=False
		# Make calls for each ltr pair
		if not checkStatusFl('{0}.GENECONVLTRs.{1}'.format(key_base, g[1:])):
			chunk_size = ceil(len(files)/procs)
			with Pool(processes=procs) as p:
				p.map(makecallMultiprocessing, geneconv_calls, chunksize=chunk_size)
			p.join()

			output = [ f for f in os.listdir(paths[AlnKey]) if f.endswith('tab') ]
			# Move geneconv files to geneconv dir
			for f in output:
				os.rename('{0}/{1}'.format(paths[AlnKey], f), '{0}/{1}'.format(paths[GENECONVgDirKey], f))

			paths['{0}.GENECONVLTRs.{1}'.format(key_base, g[1:])] = paths[GENECONVgDirKey]
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}.GENECONVLTRs.{1}\t{2}\n'.format(key_base, g[1:], paths[GENECONVgDirKey])) 
		
		append2logfile(paths['output_top_dir'], mainlogfile, 'Parsing GENECONV output')
		# Parse geneconv files
		sig = [] # sig is short for "significant". To contain significant global inner fragments identified by GENECONV
		for f in os.listdir(paths[GENECONVgDirKey]):
			if f.endswith('tab'):
				with open('{0}/{1}'.format(paths[GENECONVgDirKey], f)) as fl:
					sig += [ line for line in fl.read().split('\n') if line.startswith('GI') ]

		sig = sorted(sig, key=lambda x:x[1]) # Sort by element name, which at sig[i][1] are as: LTR_retrotransposon1;LTR_retrotransposon1
		
		paths['GENECONVsummary'] = '{0}/GENECONV_{1}.summary'.format(paths[GENECONVDirKey], g[1:])
		paths['GENECONV_output'] = '{0}/GENECONVoutput_{1}.tab'.format(paths[GENECONVDirKey], g[1:])
		IAGCpositive = set()
		with open(paths['GENECONV_output'], 'w') as outputFl:
			with open(paths['GENECONVsummary'], 'w') as summaryFl:
				summaryFl.write('# ratio is the ratio of the alignment length to the alignment length minus the gene conversion tract\n')
				try:
					summaryFl.write('# {0} elements out of {1}, or {2:.1f}% with possible evidence of gene conversion\n'.format(len(sig), num_elements, ((len(sig)/num_elements)*100)))
				except ZeroDivisionError:
					summaryFl.write('# {0} elements out of {1}, or 0% with possible evidence of gene conversion\n'.format(len(sig), num_elements))
				summaryFl.write('# 5% are expected by chance\n')
				summaryFl.write('#elementName\tsim_p-val\talnLen\tstart\tend\ttractLen\tratio\n')
				for line in sig:
					totDifs = int(line.strip().split()[9])
					if totDifs < 3:
						continue
					outputFl.write(line + '\n')
					contents = line.split('\t')
					element = contents[1].split(';')[0]
					IAGCpositive.add(element)
					sim_p = float(contents[2])
					#KA_p = float(contents[3])
					start=int(contents[4])
					end=int(contents[5])
					tractLen = int(contents[6])
					alnLen = alnLens[element]
					ratio = alnLen / (alnLen-tractLen)
					summaryFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(element, sim_p, alnLen, start, end, tractLen, ratio))

		append2logfile(paths['output_top_dir'], mainlogfile, 'GENECONV output and a summary written to:\n{0}\n{1}'.format(paths['GENECONV_output'], paths['GENECONVsummary'] ))
		paths[SummaryKey] = paths['GENECONVsummary']
		with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
			statusFlAppend.write('{0}\t{1}\n'.format(SummaryKey, paths['GENECONVsummary']))

def ltr_divergence(I=6, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, iters=1):
	'''
	Runs PAUP
	iters used with modeltest(). deprecated pretty much
	'''

	global paths

	WICKERCLUST = False
	MCLCLUST = False
	if clustering_method == 'WickerFam':
		key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		WICKERCLUST = True
		WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
		WickerLTRdivDirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_LTRdivDir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		statusFlKey = '{0}_PAUP_divergence_dir'.format(WickerLTRdivDirkey)
		paths[WickerLTRdivDirkey] = '{0}/LTR_divergence'.format(WickerDir)
		paths['DivergenceTopDir'] = paths[WickerLTRdivDirkey]
		SummaryKey = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_LTR_divergence_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
	elif clustering_method == 'MCL':
		key_base = 'MCL_I{0}'.format(I)
		MCLCLUST = True
		MCLdir = paths['MCL_I{0}'.format(I)]
		MCL_LTRdivDirkey = 'MCL_I{0}_LTRdivDir'.format(I)
		statusFlKey = '{0}_PAUP_divergence_dir'.format(MCL_LTRdivDirkey)
		paths[MCL_LTRdivDirkey] = '{0}/LTR_divergence'.format(MCLdir)
		paths['DivergenceTopDir'] = paths[MCL_LTRdivDirkey]
		SummaryKey = 'MCL_I{0}_LTR_divergence_summary'.format(I)
		TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
	else:
		sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL and it is: {0}'.format(clustering_method))
	MakeDir('PAUPdivergenceDir', '{0}/PAUP'.format(paths['DivergenceTopDir']))
	MakeDir('PAUPNexusInputDir', '{0}/nexus'.format(paths['PAUPdivergenceDir'.format(I)]))
	MakeDir('PAUPDivOutDir', '{0}/divergences'.format(paths['PAUPdivergenceDir'.format(I)]))
	if checkStatusFl(SummaryKey):
		append2logfile(paths['output_top_dir'], mainlogfile, 'ltr_divergence() already completed: {0}'.format(paths[SummaryKey]))
		return
	if LTRDIVERGENCE:
		paupCalls = []
		modeltestResults = {}
		if not checkStatusFl(statusFlKey):
			for classif in classifs:
				# parse model test results for paup block and summary for best model
				if MCLCLUST:
					ModeltestSummaryKey = 'MCLModelTestSummary_{0}_iters_I{1}_{2}'.format(iters,I, classif)
				elif WICKERCLUST:
					ModeltestSummaryKey = 'WickerModelTestSummary_{0}_iters_{1}_pId_{2}_percAln_{3}_minLen_{4}'.format(iters, WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)
				modeltestResults[classif] = {}
				# Create PAUP calls for each cluster
				if MCLCLUST:
					clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
				elif WICKERCLUST:
					clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]
				clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]

				LTRalnPths = { '_'.join(f.split('_')[:2]):'{0}/{1}'.format(paths[TrimalKey], f) for f in os.listdir(paths[TrimalKey]) if f.endswith('trimmed') } # LTR_retrotransposon148_LTRs.fasta.aln.trimmed 
				if ModeltestSummaryKey in paths:
					with open(paths[ModeltestSummaryKey]) as testOutputFl:
						paupLines = ''
						method = None
						model = None
						clust = None
						for line in testOutputFl:
							if line.startswith('BIC'):
								if not method == None:
									modeltestResults[classif][clust] = (paupLines, model)
								if len(line.strip().split()) == 17:
									method, model, pa, pc, pg, pt, kappa, titv, Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma, clust = line.strip().split()
								else:
									method, model, pa, pc, pg, pt, kappa, titv, Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma, clust, clustSize = line.strip().split()
								paupLines = ''
							else:
								paupLines += line

						if not method == None:
							modeltestResults[classif][clust] = (paupLines, model)
					for j in range(len(clusters)):
						for el in clusters[j]:
							paupBlock = ''
							seqRec = list(SeqIO.parse(LTRalnPths[el], 'fasta'))
							paupBlock = '''#NEXUS
begin DATA;
DIMENSIONS ntax=2 nchar={0};
FORMAT datatype=dna missing=? gap=-;
MATRIX
'{1}' {3}
'{2}' {4}
;
END;
'''.format(len(seqRec[0].seq), seqRec[0].id+'_L', seqRec[1].id+'_R', str(seqRec[0].seq), str(seqRec[1].seq))

							if str(j) in modeltestResults[classif]: # Have modeltest result for this cluster
								model = modeltestResults[classif][str(j)][1]
								for line in modeltestResults[classif][str(j)][0].split('\n'):
									if line.startswith('SaveDist'):
										paupBlock += 'SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};\n'.format(model, el)
									else:
										paupBlock += line + '\n'
							else: # No model test result, perhaps the cluster is too small (default may be no modeltesting for clusters <5 elements)
								if len(set(str(seqRec[0].seq))) == 2 and len(set(str(seqRec[1].seq))) == 2: # HKY+85 does not work when there are only 2 character states
									model = 'JC'
									append2logfile(paths['output_top_dir'], mainlogfile, 'Model testing not done for {0}. Using JC'.format('jModeltest2summary_{0}'.format(classif)))
									paupBlock += '''[!
Likelihood settings from best-fit model (JC) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold and there are only 2 character
states.]

BEGIN PAUP;
Dset distance=jc;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
'''.format(model, el)
								else:
									model = 'HKY85'
									append2logfile(paths['output_top_dir'], mainlogfile, 'Model testing not done for {0}. Using HKY85'.format('jModeltest2summary_{0}'.format(classif)))
									paupBlock += '''[!
Likelihood settings from best-fit model (HKY) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold.]

BEGIN PAUP;
Dset distance=hky85;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
'''.format(model, el)
							with open('{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(paths['PAUPNexusInputDir'], seqRec[0].id, model, classif, j), 'w') as nexusFl:
								nexusFl.write(paupBlock)
							paup_call = [ executables['paup'], '-n', '{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(paths['PAUPNexusInputDir'], seqRec[0].id, model, classif, j) ]
							packet = (paup_call, None, None, None)
							paupCalls.append(packet)

				else:
					for j in range(len(clusters)):
						for el in clusters[j]:
							paupBlock = ''
							seqRec = list(SeqIO.parse(LTRalnPths[el], 'fasta'))
							paupBlock = '''#NEXUS
begin DATA;
DIMENSIONS ntax=2 nchar={0};
FORMAT datatype=dna missing=? gap=-;
MATRIX
'{1}' {3}
'{2}' {4}
;
END;
'''.format(len(seqRec[0].seq), seqRec[0].id+'_L', seqRec[1].id+'_R', str(seqRec[0].seq), str(seqRec[1].seq))
							if len(set(str(seqRec[0].seq))) == 2 and len(set(str(seqRec[1].seq))) == 2: # HKY+85 does not work when there are only 2 character states
								append2logfile(paths['output_top_dir'], mainlogfile, 'Model testing not done for {0}. Using JC'.format('jModeltest2summary_{0}'.format(classif)))
								model = 'JC'
								paupBlock += '''[!
Likelihood settings from best-fit model (JC) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold and there are only 2 character
states.]

BEGIN PAUP;
Dset distance=jc;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
'''.format(model, el)
							else:
								append2logfile(paths['output_top_dir'], mainlogfile, 'Model testing not done for {0}. Using HKY85'.format('jModeltest2summary_{0}'.format(classif)))
								model = 'HKY85'
								paupBlock += '''[!
Likelihood settings from best-fit model (HKY) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold.]

BEGIN PAUP;
Dset distance=hky85;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
'''.format(model, el)
							with open('{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(paths['PAUPNexusInputDir'], seqRec[0].id, model, classif, j), 'w') as nexusFl:
								nexusFl.write(paupBlock)
							paup_call = [ executables['paup'], '-n', '{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(paths['PAUPNexusInputDir'], seqRec[0].id, model, classif, j) ]
							packet = (paup_call, None, None, None)
							paupCalls.append(packet)

			chunk_size = ceil(len(paupCalls)/procs)
			with Pool(processes=procs) as p:
				p.map(makecallMultiprocessing, paupCalls, chunksize=chunk_size)
			p.join()
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(statusFlKey, paths['PAUPdivergenceDir']))
	gcDct = {}
	if GENECONVLTRS:
		GENECONVSummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)
		with open(paths[GENECONVSummaryKey], 'r') as gcFl:
			for line in gcFl:
				if not line.startswith('#'):
					el, p, alnLen, start, end, tractLen, ratio = line.strip().split()
					ratio = float(ratio)
					alnLen = int(alnLen)
					start = int(start)
					end = int(end)

					if not el in gcDct:
						gcDct[el] = [(start, end, alnLen)]
					else:
						gcDct[el].append((start, end, alnLen))

		# Average ratio for those elements with multiple predicted GC tracts.
		for el in gcDct:
			if len(gcDct[el]) == 1:
				tractLen = gcDct[el][0][1] - gcDct[el][0][0]
				alnLen = gcDct[el][0][2]
				gcDct[el] = alnLen / (alnLen - tractLen)
				assert gcDct[el] != float(1), "GENECONV ratio of alnLen/(alnLen-tractLen) > 1. Shouldn't be the case since GENECONV found evidence of gene conversion. see ltr_divergece()"
			else:
				alnLen = gcDct[el][0][2]
				coords = []
				for tract in gcDct[el]:
					coords.append(tract[:2])
				coords = sorted(coords, key=lambda x:x[0])
				c1 = coords[0]
				tracts = []
				for c2 in coords[1:]:
					c3 = mergeCoords(c1, c2) 
					if c3[0][0] == c3[0][1]:
						c1 = c3[0][0]
					else:
						tracts.append(c3[0][0])
						c1 = c3[0][1]
				tracts.append(c1)
				tractLen = 0
				for t in tracts:
					tractLen += t[1] - t[0] + 1
				gcDct[el] = alnLen / (alnLen - tractLen)
				assert gcDct[el] != float(1), "GENECONV ratio of alnLen/(alnLen-tractLen) > 1. Shouldn't be the case since GENECONV found evidence of gene conversion. see ltr_divergece()"
	# Read PAUP output
	clustLenDcts = {}
	clustDct = {}
	for classif in classifs:
		clustDct[classif] = {}
		if MCLCLUST:
			clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
		elif WICKERCLUST:
			clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]
		clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]
		for j in range(len(clusters)):
			for el in clusters[j]:
				clustDct[classif][el] = j

		clustLenDct = {i:len(clusters[i]) for i in range(len(clusters))}
		clustLenDcts[classif] = clustLenDct

	paths['DivergenceSummary'] = '{0}/LTR_divergences.tab'.format(paths['DivergenceTopDir'])

	with open(paths['DivergenceSummary'.format(I)],'w') as divSummaryFl:
		divSummaryFl.write('elementName\tclassification\tMCLinflationValue\tcluster\tclusterSize\tmodel\tdivergence\tcorrectedDivergence\tIntraelementGeneConversion\n')
	for fname in os.listdir(paths['PAUPDivOutDir']):
		fPth = '{0}/{1}'.format(paths['PAUPDivOutDir'], fname)
		el1, el2, div = open(fPth, 'r').read().strip().split('\t')
		el = el1[:-2]
		div = float(div)
		model = fname.split('.')[1]
		if list(clustDct.keys()) == ['All']:
			classif == 'All'
		else:
			classif = classifs_by_element[el]

		if classif in clustDct:
			if el in clustDct[classif]:
				clust = clustDct[classif][el]
			else:
				continue
		else:
			continue
		clustSize = clustLenDcts[classif][clust]
		if el in gcDct:
			divc = div*gcDct[el]
			GC = 'Yes'
		else:
			divc = div
			GC = 'No'
	
		# Possibly exclude proportion of invariable sites in model, it was part of the models that estimated too high
		if div > 5 or divc > 5:
			print('{0} LTRs have estimated substitutions per site of {1} and {2} (gene conversion-corrected) using {3}. It was excluded from the summary table at: {4}'.format(el, div, divc, model, paths['DivergenceSummary']), file=sys.stderr)
			continue

		# Combine as new output
		with open(paths['DivergenceSummary'],'a') as divSummaryFl:
			divSummaryFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(el, classif, I, clust, clustSize, model, div, divc, GC))
							
	if not SummaryKey in paths:
		with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
			statusFlAppend.write('{0}\t{1}\n'.format(SummaryKey, paths['DivergenceSummary']))


def phylo(removegeneconv=True, BOOTSTRAP=True, I=6, align='cluster', removehomologouspair=True, part='entire', clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, auto_outgroup=False, bootstrap_reps=100, minClustSize=4, convert_to_ultrametric=False, bpflank=None, combine_and_do_small_clusters=True):
	'''
	align: one of cluster or classifs
	I only applies to cluster

	removegeneconv	If True GENECONV output will be parsed and elements with evidence of gene conversion will be removed from cluster prior to aligning for tree
	BOOTSTRAP=False is not implemented yet.
	I		inflation parameter used for MCL clusters desired to use for this function
	align		'cluster'  align and make trees for each WickerFam or MCL cluster, depending on setting of clustering_method
			'classif'  align and make trees for each superfamily
	removehomologouspair	Not implementd for align='cluster' yet.
				If True then the flanking X bp on both sides of the element are extracted and blastn'd together
				(constraints can be changed, see ). One of each pair of hits is removed prior to alignment.
	part		'entire'  Use whole sequence of element for trees. Not including TSD. From: start(LTR1) -> end(LTR2)
			'inernal' Use sequence between LTRs for trees. From: end(LTR1) -> start(LTR2)
	clustering_method	'WickerFam'
				'MCL'
	WickerParams	A dictionary of the format {'pId':A,'percAln':B,'minLen':C} where A, B, C are numbers specifying the WickerFam() params.
	auto_outgroup	Not available for align='classif'.
			# Automatically roots each cluster's tree with an outgroup by this process:
			#   the outgroup shall be a random element from cluster k
			#   where cluster k is the largest of the clusters that is not j
			#   if j is the first cluster then the next smallest cluster is k
			#   if there is no other cluster, no outgroup is used
			
	'''
	global paths

	AUTO_OUTGROUP = auto_outgroup
	REMOVEHOMOLOGOUSFLANK = removehomologouspair
	REMOVEGENECONV = removegeneconv
	ULTRAMETRIC = convert_to_ultrametric

	if AUTO_OUTGROUP:
		strOutgroup = 'withOutgroup'
	else:
		strOutgroup = 'noOutgroup'
	if REMOVEHOMOLOGOUSFLANK:
		strHomoflank = 'homoflank'
	else:
		strHomoflank = 'nohomoflank'

	WICKERCLUST = False
	MCLCLUST = False
	if clustering_method == 'WickerFam':
		WICKERCLUST = True
		WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
		WickerTreesDirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_TreesDir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
		paths[WickerTreesDirkey] = '{0}/Trees'.format(WickerDir)
		paths['TreesDir'] = paths[WickerTreesDirkey]
		MakeDir(WickerTreesDirkey, paths['TreesDir'])

	elif clustering_method == 'MCL':
		MCLCLUST = True
		MCLdir = paths['MCL_I{0}'.format(I)]
		MCLTreesDirkey = 'MCL_I{0}_TreesDir'.format(I)
		paths[MCLTreesDirkey] = '{0}/Trees'.format(MCLdir)
		paths['TreesDir'] = paths[MCLTreesDirkey]
		MakeDir(MCLTreesDirkey, paths['TreesDir'])
	else:
		sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL, and it is: {0}'.format(clustering_method))

	OutPth = None
	alnPthKeys = []

	ENTIRE = False
	INTERNAL = False
	if part == 'entire':
		ENTIRE = True
	elif part == 'internal':
		INTERNAL = True
	else:
		sys.exit('part parameter to AutoAlign() is invalid: {0}. Must be either "entire" or "internal"'.format(part))

	if INTERNAL:
		MakeDir('InternalRegionsAlignments', '{0}/InternalRegions'.format(paths['TreesDir']))
		paths['RegionDir'] = paths['InternalRegionsAlignments']
	elif ENTIRE:
		MakeDir('WholeElementAlignments', '{0}/WholeElements'.format(paths['TreesDir']))
		paths['RegionDir'] = paths['WholeElementAlignments']
	if align == 'classif': # OUTGROUP is not possible here
		MakeDir('Superfamilies', '{0}/Superfamilies'.format(paths['RegionDir']))
		if REMOVEHOMOLOGOUSFLANK:
			strHomoflank = 'homoflank'
			MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(paths['Superfamilies']))
		else:
			strHomoflank = 'nohomoflank'
			MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['Superfamilies']))
		OutPth = paths['HomoFlankDir']
		if WICKERCLUST:
			AutoAlign(I=None, part=part, rmgeneconv=removegeneconv, minClustSize=minClustSize, align='classif', rmhomologflank=REMOVEHOMOLOGOUSFLANK, clustering_method='WickerFam', WickerParams={'pId':WickerParams['pId'],'percAln':WickerParams['percAln'],'minLen':['minLen']}, auto_outgroup=AUTO_OUTGROUP, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
			alnPthKeys.append('WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_all.{4}.{5}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, strHomoflank))
		elif MCLCLUST:
			AutoAlign(I=I, part=part, rmgeneconv=removegeneconv, minClustSize=minClustSize, align='classif', rmhomologflank=REMOVEHOMOLOGOUSFLANK, clustering_method='MCL', WickerParams=None, auto_outgroup=AUTO_OUTGROUP, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)
			alnPthKeys.append('Aln_{0}_I{1}_cluster_all.{2}.{3}'.format(classif, I, strHomoflank))
	elif align == 'cluster':
		if REMOVEGENECONV:
			MakeDir('GCDir', '{0}/GeneconversionDisallowed'.format(paths['RegionDir']))
		else:
			MakeDir('GCDir', '{0}/NoGCFiltering'.format(paths['RegionDir']))
		if REMOVEHOMOLOGOUSFLANK:
			strHomoflank = 'homoflank'
			MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(paths['GCDir']))
		else:
			strHomoflank = 'nohomoflank'
			MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['GCDir']))
		if AUTO_OUTGROUP:
			strOutgroup = 'withOutgroup'
			MakeDir('OutgroupDir', '{0}/WithOutgroup'.format(paths['HomoFlankDir']))
		else:
			strOutgroup = 'noOutgroup'
			MakeDir('OutgroupDir', '{0}/NoOutgroup'.format(paths['HomoFlankDir']))
		OutPth= paths['OutgroupDir']

		if WICKERCLUST:
			AutoAlign(I=None, part=part, rmgeneconv=removegeneconv, minClustSize=minClustSize, align='clusters', rmhomologflank=REMOVEHOMOLOGOUSFLANK, clustering_method='WickerFam', WickerParams={'pId':WickerParams['pId'],'percAln':WickerParams['percAln'],'minLen':WickerParams['minLen']}, auto_outgroup=AUTO_OUTGROUP, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)

		elif MCLCLUST:
			AutoAlign(I=I, part=part, rmgeneconv=removegeneconv, minClustSize=minClustSize, align='clusters', rmhomologflank=REMOVEHOMOLOGOUSFLANK, clustering_method='MCL', WickerParams=None, auto_outgroup=AUTO_OUTGROUP, bpflank=bpflank, combine_and_do_small_clusters=combine_and_do_small_clusters, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)

	for classif in classifs:
		# Read clusters
		if MCLCLUST:
			clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
		elif WICKERCLUST:
			clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]
		clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]
		# Read GENECONV output
		if REMOVEGENECONV:
			gc = 'GCfiltered'
			append2logfile(paths['output_top_dir'], mainlogfile, 'Excluding {0} elements with evidence of gene conversion'.format(classif))
			if MCLCLUST:
				gcSummaryPth = 'MCL_I{0}_GENECONV_summary'.format(I)
			elif WICKERCLUST:
				gcSummaryPth = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
			if gcSummaryPth in paths:
				if os.path.isfile(gcSummaryPth):
					gcDct = {}
					with open(paths[gcSummaryPth], 'r') as gcFl:
						for line in gcFl:
							el, clust, classif, g = line.strip().split()
							el = 'LTR_retrotransposon{0}'.format(el)
							clust = int(clust)
							if classif in gcDct:
								if clust in gcDct[classif]:
									if not el in gcDct[classif][clust]:
										gcDct[classif][clust].append(el)
								else:
									gcDct[classif][clust] = [el]
							else:
								gcDct[classif] = {clust:[el]}
			elif gcSummaryPth not in paths:
				print('phylo(): modeltest(removegeneconv=True) used but {0} not in paths or status file. This happens at least when there is only 1 element with a classification and therefore nothing is aligned'.format(gcSummaryPth), file=sys.stderr)

			if os.path.isfile(gcSummaryPth):
				for j in range(len(clusters)):
					if j in gcDct[classif]:
						clusters[j] = [ el for el in clusters[j] if not el in gcDct[classif][j] ]
		else:
			gc = 'NoGCfiltering'
		
		# Begin cluster routine
		if align == 'cluster':
			for j in range(len(clusters)):
				if WICKERCLUST:
					alnPth = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_{4}_{5}.{6}.{7}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, j, gc, strHomoflank, strOutgroup)
				elif MCLCLUST:
					alnPth = 'Aln_{0}_I{1}_cluster{2}_{3}.{4}.{5}'.format(classif, I, j, gc, strHomoflank, strOutgroup)
				if  alnPth in paths:
					alnPthKeys.append(alnPth)

			# clustersmall will be left out of OUTGROUP. No DTT for clustersmall
			if combine_and_do_small_clusters:
				if WICKERCLUST:
					alnPth = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_{4}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, gc)
					if alnPth in paths:
						alnPthKeys.append(alnPth)
				elif MCLCLUST:
					alnPth = 'Aln_{0}_I{1}_clustersmall_{2}'.format(classif, I, gc)
					if alnPth in paths:
						alnPthKeys.append(alnPth)

	if BOOTSTRAP:
		append2logfile(paths['output_top_dir'], mainlogfile, 'Ready to begin bootstrapping for:\n{0}'.format('\n'.join(alnPthKeys)))
		if ULTRAMETRIC:
			if WICKERCLUST:
				bootstrap(alnPthsLst=alnPthKeys, reps=bootstrap_reps, OutPth=OutPth, convert_to_ultrametric=True, WickerParams={'pId':WickerParams['pId'],'percAln':WickerParams['percAln'],'minLen':WickerParams['minLen']}, gc=gc, strHomoflank=strHomoflank, strOutgroup=strOutgroup, I=None)
			elif MCLCLUST:
				bootstrap(alnPthsLst=alnPthKeys, reps=bootstrap_reps, OutPth=OutPth, convert_to_ultrametric=True, WickerParams=None,  gc=gc, strHomoflank=strHomoflank, strOutgroup=strOutgroup, I=I)
		else:
			if WICKERCLUST:
				bootstrap(alnPthsLst=alnPthKeys, reps=bootstrap_reps, OutPth=OutPth, convert_to_ultrametric=False, WickerParams={'pId':WickerParams['pId'],'percAln':WickerParams['percAln'],'minLen':WickerParams['minLen']}, gc=gc, strHomoflank=strHomoflank, strOutgroup=strOutgroup, I=None)
			elif MCLCLUST:
				bootstrap(alnPthsLst=alnPthKeys, reps=bootstrap_reps, OutPth=OutPth, convert_to_ultrametric=False, WickerParams=None, gc=gc, strHomoflank=strHomoflank, strOutgroup=strOutgroup, I=I)


def SeqbootCall(Call):
	baseDir = os.getcwd()
	os.chdir(Call[1])
	subprocess.call(Call[0], stdin=open('seqboot.conf', 'r'))
	os.chdir(baseDir)


def bootstrap(alnPthsLst, reps, OutPth=None, convert_to_ultrametric=False, WickerParams=None, gc=None, strHomoflank=None, strOutgroup=None, I=None):

	global paths

	append2logfile(paths['output_top_dir'], mainlogfile, 'Began bootstrap()')
	OutgroupSummaryKey = None
	seqbootCalls = []
	ULTRAMETRIC = convert_to_ultrametric

	for AlnFasta in alnPthsLst:
		if not os.path.isfile(paths[AlnFasta]):
			continue
		append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing to run SEQBOOT for:\n{0}'.format(paths[AlnFasta]))
		AlignIO.convert(in_file=paths[AlnFasta], out_file='{0}.phylip'.format(paths[AlnFasta]), in_format='fasta', out_format='phylip')
		SeqBootInstructionsFlPth = '{0}/seqboot.conf'.format('/'.join(paths[AlnFasta].split('/')[:-1]))
		with open(SeqBootInstructionsFlPth, 'w') as SeqBootInstructionsFl:
			Phylip = '{0}.phylip'.format(paths[AlnFasta]).split('/')[-1]
			SeqBootInstructionsFl.write('{0}\nR\n{1}\nY\n1\nR\n'.format(Phylip, reps))
		
		seqbootCalls.append([['{0}/seqboot'.format(executables['phylip'])], '/'.join(paths[AlnFasta].split('/')[:-1])])

	chunk_size = ceil(len(seqbootCalls)/procs)
	with Pool(processes=procs) as p:
		p.map(SeqbootCall, seqbootCalls, chunksize=chunk_size)
	p.join()

	append2logfile(paths['output_top_dir'], mainlogfile, 'Finished running SEQBOOT')

	fasttreeCalls = []
	for AlnFasta in alnPthsLst:
		alignment_length = len(list(SeqIO.parse(paths[AlnFasta], 'fasta'))[0])
		Alns = list(SeqIO.parse(paths[AlnFasta], 'fasta'))
		numAlns = len(Alns)
		alnLen = len(Alns[0].seq)
		pth = '/'.join(paths[AlnFasta].split('/')[:-1])
		clustMethod = paths[AlnFasta].split('/')[1]
		if clustMethod.startswith('Wicker'):
			clustMethod = 'WickerFam'
		if clustMethod == 'MCL':
			#LTRan_output/MCL/I8/Alignments/WholeElements/NoFiltering/Copia/cluster_0/elements.fasta.aln.trimal
			settings  = paths[AlnFasta].split('/')[2]
			classif = paths[AlnFasta].split('/')[-3]
			if classif == 'whole_classif':
				classif = paths[AlnFasta].split('/')[-2]
			clust = paths[AlnFasta].split('/')[-2].split('_')[-1]
			if clust.endswith('.fasta.aln.trimal'):
				clust = 'wholeClassif'
		elif clustMethod == 'WickerFam':
			#['LTRan_output', 'WickerFamDir', '80_pId_80_percAln_80_minLen', 'Alignments', 'WholeElements', 'NoFiltering', 'Other', 'cluster_small', 'elements.fasta.aln.trimal']
			classif = paths[AlnFasta].split('/')[-3]
			clust = paths[AlnFasta].split('/')[-2]
			settings = paths[AlnFasta].split('/')[2]
		else:
			sys.exit('modeltest() parameter clustering_method needs to be either WickerFam or MCL, and it is, from within bootstrap(): {0}'.format(clustMethod))

		if ULTRAMETRIC:
			if clustMethod == 'WickerFam':
				OutgroupSummaryKey = 'WickerOutgroups_{0}_pId_{1}_percAln_{2}_minLen_{3}_{4}.{5}.{6}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, gc, strHomoflank, strOutgroup)
			elif clustMethod == 'MCL':
				OutgroupSummaryKey = 'MCLOutgroups_{0}_I{1}_{2}.{3}.{4}'.format(classif, I, gc, strHomoflank, strOutgroup)
		#LTRan_output/Alignments/WholeElements/GeneconversionDisallowed/whole_classif/NoHomologousFlank/All/elements.fasta.aln.trimal
		seqbootOutputPhylip = '{0}/outfile'.format(pth)
		# Split seqboot multi-phylip output
		MakeDir('classifDir', '{0}/{1}'.format(OutPth, classif))
		MakeDir('clustDir', '{0}/{1}'.format(paths['classifDir'], clust))
		append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing to run FastTree for:\n{0}'.format(paths[AlnFasta]))
		# Integrate both Wicker and MCL
		runID = '{0}_{1}_Bootstrap_{2}_{3}'.format(clustMethod, settings, classif, clust)
		# Run FastTree for main tree
		if not 'mainTree_{0}'.format(runID) in paths:
			mainTree = '{0}/{1}.main_tree'.format(paths['clustDir'], paths[AlnFasta].split('/')[-1])
			makecallMultiprocessing(([ executables['fasttree'], '-nt', '-gtr' ], mainTree, '{0}/fasttree.err'.format(pth), paths[AlnFasta] ))
			paths['mainTree_{0}'.format(runID)] = mainTree
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format('mainTree_{0}'.format(runID), mainTree))
			append2logfile(paths['output_top_dir'], mainlogfile, 'Generated main tree:\n{0}'.format(mainTree))

		if not 'Trees_{0}'.format(runID) in paths:
			with open(seqbootOutputPhylip, 'r') as alnIters:
				replicate = 0
				repFlPth = None
				append2logfile(paths['output_top_dir'], mainlogfile, 'Preparing for bootstrapping')
				for line in alnIters:
					if line.strip().replace(' ', '') == '{0}{1}'.format(numAlns, alnLen): # new alignment
						
						if repFlPth != None:
							iterDir = '{0}/{1}'.format(paths['clustDir'], replicate)
							fasttree_call = ( [ executables['fasttree'], '-nt', '-gtr' ], '{0}/tree.newick'.format(iterDir), '{0}/fasttree.err'.format(iterDir), repFlPth )
							#fasttree_call_string =  '{0} -nt -gtr <{1} >{2} 2>{3}'.format(executables['fasttree'], repFlPth, '{0}/tree.newick'.format(iterDir), '{0}/fasttree.err'.format(iterDir))
							fasttreeCalls.append(fasttree_call)
						replicate += 1
						# MAKE DIR FOR BOOTSTRAP REPLICATES OUTPUT
						MakeDir('iterDir', '{0}/{1}'.format(paths['clustDir'], replicate))
						repFlPth = '{0}/alnReplicate_{1}.fasta'.format(paths['iterDir'], replicate)
						with open(repFlPth, 'a') as repFl:
							repFl.write(line)
					else:
						with open(repFlPth, 'a') as repFl:
							repFl.write(line)
						
				fasttree_call = ( [ executables['fasttree'], '-nt', '-gtr' ], '{0}/tree.newick'.format(paths['iterDir']), '{0}/fasttree.err'.format(paths['iterDir']), repFlPth )
				#fasttree_call_string =  '{0} -nt -gtr <{1} >{2} 2>{3}'.format(executables['fasttree'], repFlPth, '{0}/tree.newick'.format(iterDir), '{0}/fasttree.err'.format(iterDir))
				fasttreeCalls.append(fasttree_call)

			chunk_size = ceil(len(fasttreeCalls)/procs)
			with Pool(processes=procs) as p:
				p.map(makecallMultiprocessing, fasttreeCalls, chunksize=chunk_size)
			p.join()
			
			paths['Trees_{0}'.format(runID)] = '{0}/allReplicates.newick'.format(paths['clustDir'])
			for d in [D for D in os.listdir(paths['clustDir']) if os.path.isdir('{0}/{1}'.format(paths['clustDir'], D))]:
				for f in os.listdir('{0}/{1}'.format(paths['clustDir'], d)):
					if f.endswith('newick'):
						tree = open('{0}/{1}/{2}'.format(paths['clustDir'], d, f), 'r').read().strip()
						with open(paths['Trees_{0}'.format(runID)], 'a') as treeFl:
							treeFl.write('{0}\n'.format(tree))

			append2logfile(paths['output_top_dir'], mainlogfile, 'All bootstrap replicate trees in this file:\n{0}'.format(paths['Trees_{0}'.format(runID)]))
			paths[runID] = paths['clustDir']
			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format(runID,  paths[runID]))

			with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
				statusFlAppend.write('{0}\t{1}\n'.format('Trees_{0}'.format(runID), paths['Trees_{0}'.format(runID)]))

		# Get bootstrap values
		compare2bootstrap_call = [ executables['perl'], '{0}/{1}'.format( paths['scriptsDir'], 'CompareToBootstrap.pl'), '-tree', paths['mainTree_{0}'.format(runID)], '-boot', paths['Trees_{0}'.format(runID)]]
		makecall(compare2bootstrap_call, stdout = '{0}.bootstrapped'.format(paths['mainTree_{0}'.format(runID)]))

		bootstrapped = '{0}.bootstrapped'.format(paths['mainTree_{0}'.format(runID)])

		# run PATHd8 to convert tree to ultrametric
		if ULTRAMETRIC:
			append2logfile(paths['output_top_dir'], mainlogfile, 'Beginning Ultrametric transformation for {0}'.format(runID))
			with open(bootstrapped, 'r') as bstreefl:
				tree = bstreefl.read().strip()
			outgroup = None
			with open(paths[OutgroupSummaryKey], 'r') as outgroup_file:
				for line in outgroup_file:
					if line == 'Alignment\toutgroup\toutgroup_cluster\n':
						continue
					key, outG, outgroup_clust = line.strip().split('\t')
					if key == AlnFasta:
						outgroup = outG[19:]

			if outgroup == None:
				print('Unable to root {0}. No ultrametric tree will be made'.format(AlnFasta))
			else:
				other_taxon = random.choice(list(SeqIO.parse(paths[AlnFasta], 'fasta'))).id # a random element from the current cluster to tell PATHd8 where the mrca is (outgroupXother_taxon) for determining relative branch lengths
				append2logfile(paths['output_top_dir'], mainlogfile, 'Beginning PATHd8 for {0}'.format(runID))
				pathd8_file_str = '''Sequence length={0};
mrca: {1}, {2}, fixage=1;
{3}
	'''.format(alignment_length, outgroup, other_taxon, tree)
				pathd8flpth = '{0}/pathd8.in'.format(paths[runID])
				with open(pathd8flpth, 'w') as outFl:
					outFl.write(pathd8_file_str)
				paths['PATHd8_output_{0}'.format(runID)] = '{0}.pathd8_ultrametric'.format(bootstrapped)
				pathd8_call = [ executables['pathd8'], pathd8flpth, paths['PATHd8_output_{0}'.format(runID)] ]
				makecall(pathd8_call)
				if os.path.isfile(paths['PATHd8_output_{0}'.format(runID)]):
					with open(paths['PATHd8_output_{0}'.format(runID)], 'r') as d8Fl:
						for line in d8Fl:
							if line.startswith('d8 tree '):
								d8tree = line[line.index('('):]
								newLocation = '{0}/{1}_{2}.bootstrapped.pathd8_ultrametric.outgroup_{3}.newick'.format(paths['classifDir'], classif, clust, outgroup)
								if clustMethod == 'WickerFam':
									newLocationKey = 'WickerOutgroups_{0}_pId_{1}_percAln_{2}_minLen_{3}_{4}.{5}.outgroup_{6}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif, gc, strHomoflank, outgroup)
								elif clustMethod == 'MCL':
									newLocationKey = 'MCLOutgroups_{0}_I{1}_{2}.{3}.outgroup_{4}'.format(classif, I, gc, strHomoflank, outgroup)
								with open(newLocation, 'w') as outFl:
									outFl.write(d8tree)
								if not checkStatusFl(newLocationKey):
									with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFl:
										statusFl.write('{0}\t{1}\n'.format(newLocationKey, newLocation))

							
				append2logfile(paths['output_top_dir'], mainlogfile, 'Finished PATHd8 for {0}'.format(runID))

		bootstrappedLocation = '{0}/{1}_{2}.bootstrapped.newick'.format(paths['classifDir'], classif, clust)
		copyfile(bootstrapped, bootstrappedLocation)
		append2logfile(paths['output_top_dir'], mainlogfile, 'Finished bootstrapping. Tree with bootstrap support values here:\n{0}'.format(bootstrappedLocation))


def div2Rplots(I=6):

	R_wd = os.path.abspath('/'.join(paths['DivergenceSummary_I{0}'.format(I)].split('/')[:-1]))
	fl = os.path.abspath(paths['DivergenceSummary_I{0}'.format(I)])
	call = [ executables['rscript'], '{0}/div2density.R'.format(paths['scriptsDir']), R_wd, fl ]
	call_str = '{0} {1} {2} {3}'.format(executables['rscript'], '{0}/div2density.R'.format(paths['scriptsDir']), R_wd, fl)
	append2logfile(paths['output_top_dir'], mainlogfile, 'Making R plot\n{0}'.format(call_str))
	with open('{0}/RscriptCalls.txt'.format(paths['output_top_dir']), 'a') as RcallFl:
		RcallFl.write(call_str+'\n')
	subprocess.call(call)


def geneconv2circoslinks(geneconvfile, ltrharvestgff, outfile, append=False, output='file', linksdct=None, transposeLinks=True):
	'''
	Converts GI tract pairs from geneconvClusters() output and writes a links file for Circos
	The GFF3 is needed to get the scaffold name.
	seqlengths needs to be a dictionary with the lengths of the sequences whose names correspond
	to the sequence names in the gff for the features with gene conversion tracts.
	Assumes LTR_retrotransposon features were used.
	append=True will append to the outfile if it exists.

	output	'file', or 'return'. If file, the links will be written to a file at outfile. If 'return'
				     then the links will be returned as a dictionary like:

				     links = { orf1:[linkline1, linkline2, linkline3], ... }

				     Only the first orf in each link pair needs to be considered because gene conversion
				     was evaluated within clusters only.
	
	linksdct	If 'output'=return and a links dct is provided here it will be updated and returned.

	transposeLinks=True will write links for elements in their position on the scaffolds.
	transposeLinks=False will write links for elements with a start position of 0 at one end of the element.
	'''
	global paths

	#lengths = { l[0]:int(l[1]) for  l in open(seqlengths).read().strip().split('\n')}

	seqs = {}
	starts = {}
	if append:
		mode = 'a'
	else:
		mode = 'w'
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
	
	if output == 'file':
		with open(outfile, mode) as outFl:
			with open(geneconvfile, 'r') as inFl:
				for line in inFl:
					if line.startswith('GI'):
						rec = line.strip().split('\t')
						el1, el2 = [ 'LTR_retrotransposon{0}'.format(e[1:]) for e in rec[1].split(';') ]
						if transposeLinks:
							el1start  = int(rec[7]) + starts[el1] - 1
							el1end  = int(rec[8]) + starts[el1] - 1
							el2start  = int(rec[10]) + starts[el2] - 1
							el2end  = int(rec[11]) + starts[el2] - 1
							el1seq = copy(seqs[el1])
							el2seq =copy(seqs[el2])
						else:
							el1start  = int(rec[7])
							el1end  = int(rec[8])
							el2start  = int(rec[10])
							el2end  = int(rec[11])
							el1seq = copy(el1)
							el2seq = copy(el2)
						if 'g0.summary' in geneconvfile:
							outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end))
						elif 'g1.summary' in geneconvfile:
							outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end))
						elif 'g2.summary' in geneconvfile:
							outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end))
                     
	elif output == 'return':
		if linksdct != None:
			links = linksdct
			links_untransposed = linksdct
		else:
			links = {}
			links_untransposed = {}
		with open(geneconvfile, 'r') as inFl:
			
			el1start, el1end, el2start, el2end, el1seq, el2seq = [None]*6

			for line in inFl:
				if line.startswith('GI'):
					rec = line.strip().split('\t')
					el1, el2 = [ 'LTR_retrotransposon{0}'.format(e[1:]) for e in rec[1].split(';') ]
					if transposeLinks:
						el1start  = int(rec[7]) + starts[el1] - 1
						el1end  = int(rec[8]) + starts[el1] - 1
						el2start  = int(rec[10]) + starts[el2] - 1
						el2end  = int(rec[11]) + starts[el2] - 1
						el1seq = copy(seqs[el1])
						el2seq = copy(seqs[el2])
						# Different colored links for different gscale parameters. g values > 2 are possible but not implemented.
						if 'g0.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end)
							g = 'g0'
						elif 'g1.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end)
							g = 'g1'
						elif 'g2.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(el1seq, el1start, el1end, el2seq, el2start, el2end)
							g = 'g2'
						if g in links:
							if el1 in links[g]:
								links[g][el1].append(outline)
							else:
								links[g][el1] = [outline]
						else:
							links[g] = {el1:[outline]}
					else:
						el1start  = int(rec[7])
						el1end  = int(rec[8])
						el2start  = int(rec[10])
						el2end  = int(rec[11])
						# Different colored links for different gscale parameters. g values > 2 are possible but not implemented.
						if 'g0.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(el1, el1start, el1end, el2, el2start, el2end)
							g = 'g0'
						elif 'g1.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(el1, el1start, el1end, el2, el2start, el2end)
							g = 'g1'
						elif 'g2.summary' in geneconvfile:
							outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(el1, el1start, el1end, el2, el2start, el2end)
							g = 'g2'
						if g in links:
							if el1 in links[g]:
								links[g][el1].append(outline)
							else:
								links[g][el1] = [outline]
						else:
							links[g] = {el1:[outline]}
		return links


def circosMultiprocessing(packet):
	'''
	For using with multiprocessing.Pool()
	'''
	global paths

	circosdir = packet[0]
	circos_call = packet[1]
	classif = packet[2]
	i = packet[3]
	outpath = packet[4]
	G = packet[5]
	current_wd = os.getcwd()
	os.chdir(circosdir)
	subprocess.call(circos_call, stdout=open('out', 'w'), stderr=open('err','w'))
	os.chdir(current_wd)
	png = '{0}/circos.png'.format(circosdir)
	svg = '{0}/circos.svg'.format(circosdir)
	newpng = '{0}/{1}.cluster_{2}.geneconv_{3}.png'.format(outpath, classif, i, '_'.join(G))
	newsvg = '{0}/{1}.cluster_{2}.geneconv_{3}.svg'.format(outpath, classif, i, '_'.join(G))
	if not os.path.isfile(newpng):
		if os.path.isfile(png):
			copyfile(png, newpng)
	if not os.path.isfile(newsvg):
		if os.path.isfile(svg):
			copyfile(svg, newsvg)


def Circos(window='1000000', plots='clusters', I=6, clustering_method='WickerFam', WickerParams={'pId':80,'percAln':80,'minLen':80}, g='g0,g1,g2'):
	'''
	Generate a Circos plot for each cluster, showing gene interelement gene conversion tracts
	by using links.
	plots needs to be either 'clusters' or 'classifs'
	'''
	global paths

	G = g.split(',')
	CLASSIFS = False
	CLUSTERS = False
	WICKERCLUST = False
	MCLCLUST = False
	if plots == 'classifs':
		CLASSIFS = True
	elif plots == 'clusters':
		CLUSTERS = True
	if clustering_method == 'WickerFam':
		WICKERCLUST = True
	elif clustering_method == 'MCL':
		MCLCLUST = True
	append2logfile(paths['output_top_dir'], mainlogfile, 'Beginning making Circos plots')
	MakeDir('CircosTopDir', '{0}/Circos'.format(paths['output_top_dir']))

	# Get seq lengths
	scafLengthsFlPth = '{0}/seqLengths.tab'.format(paths['CircosTopDir'])
	scafLengths = {}
	with open(paths['inputFasta'], 'r') as inFl:
		currentSeqName = None
		for line in inFl:
			if line.startswith('>'):
				currentSeqName = line.strip()[1:].split(' ')[0]
				scafLengths[currentSeqName] = 0
			else:
				scafLengths[currentSeqName] += len(line.strip())
	
	with open(scafLengthsFlPth, 'w') as outFl:
		for scaf in sorted(list(scafLengths.keys())):
			outFl.write('{0}\t{1}\n'.format(scaf, str(scafLengths[scaf])))

	# Create ideogram file!
	allscafs = '{0}/seqs.track'.format(paths['CircosTopDir'])
	ideogramCall = [ '{0}/ideogramFromLengths.py'.format(paths['scriptsDir']) ]
	makecall(ideogramCall, stdin=scafLengthsFlPth, stdout=allscafs)

	# Make track for elements!
	if CLASSIFS:
		# Separate out GFFs by classif
		allGFFoutPth = '{0}/all.gff'.format(paths['CircosTopDir'])
		with open(paths['CurrentGFF']) as gffFl:
			for line in gffFl:
				if '\tLTR_retrotransposon\t' in line:
					gffLine = GFF3_line(line)
					classif = classifs_by_element[gffLine.attributes['ID']]
					MakeDir('classifDir', '{0}/{1}'.format(paths['CircosTopDir'], classif))
					GFFoutPth = '{0}/{1}.gff'.format(paths['classifDir'], classif)
					with open(GFFoutPth, 'a') as GFFoutFl:
						GFFoutFl.write(line)
					with open(allGFFoutPth, 'a') as GFFoutFl:
						GFFoutFl.write(line)
		append2logfile(paths['output_top_dir'], mainlogfile, 'Created GFF files for each classification for converting to Circos heatmap tracks.')

		for classif in classifs:
			classifDir = '{0}/{1}'.format(paths['CircosTopDir'], classif)
			GFFoutPth = '{0}/{1}.gff'.format(classifDir, classif)
			gff2heatmapCall = [ '{0}/gff2circos-heatmap.py'.format(paths['scriptsDir']), '-gff', GFFoutPth, '-window', window, '-scafLens', scafLengthsFlPth ]
			makecall(gff2heatmapCall, stdout='{0}/{1}.heatmap.track'.format(classifDir, classif))

			# Geneconv output to circos links track
			paths['GENECONV_{0}_dir'.format(classif)] = '{0}/{1}'.format(paths[geneconvOutputDir], classif)
			outfile = '{0}/{1}.testlinks'.format(paths['CircosTopDir'], classif)
			g0fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g0')
			g1fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g1')
			g2fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g2')
			if os.path.isfile(g0fl):
				# Convert GENECONV output to Circos links track
				geneconv2circoslinks(g0fl, paths['CurrentGFF'], outfile)
			if os.path.isfile(g1fl):
				# Convert GENECONV output to Circos links track
				geneconv2circoslinks(g1fl, paths['CurrentGFF'], outfile, append=True)
			if os.path.isfile(g2fl):
				# Convert GENECONV output to Circos links track
				geneconv2circoslinks(g2fl, paths['CurrentGFF'], outfile, append=True)
			append2logfile(paths['output_top_dir'], mainlogfile, 'Created links tracks for Circos from intra-cluster inter-element GENECONV output')

		gff2heatmapCall = [ '{0}/gff2circos-heatmap.py'.format(paths['scriptsDir']), '-gff', allGFFoutPth, '-window', window, '-scafLens', scafLengthsFlPth ]
		makecall(gff2heatmapCall, stdout='{0}/all.heatmap.track'.format(paths['CircosTopDir']))
		append2logfile(paths['output_top_dir'], mainlogfile, 'Converted GFFs to heatmap tracks for Circos')

	elif CLUSTERS:

		if WICKERCLUST:
			WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])]
			paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])] = '{0}/GENECONV'.format(WickerDir)
			WickerGCdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'])
			if not checkStatusFl(WickerGCdirkey):
				sys.exit('Circos() not possible: geneconvClusters() not done yet.')
			geneconvOutputDir = WickerGCdirkey
			MakeDir('CurrentTopDir', '{0}/WickerFam_{1}_pId_{2}_percAln_{3}_minLen'.format(paths['CircosTopDir'], WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen']))
		elif MCLCLUST:
			MCLdir = paths['MCL_I{0}'.format(I)]
			paths['MCL_I{0}_GENECONVdir'.format(I)] = '{0}/GENECONV'.format(MCLdir)
			if not checkStatusFl('MCL_I{0}_GENECONVdir'.format(I)):
				sys.exit('Circos() not possible: geneconvClusters() not done yet.')
			geneconvOutputDir = 'MCL_I{0}_GENECONVdir'.format(I)
			MakeDir('CurrentTopDir', '{0}/MCL_I{1}'.format(paths['CircosTopDir'], I))

		# Create a Circos plot for each cluster
		heatmapcalls = []
		tilecalls = []
		circoscalls = []
		for classif in classifs:

			# Geneconv output to circos links track
			paths['GENECONV_{0}_dir'.format(classif)] = '{0}/{1}'.format(paths[geneconvOutputDir], classif)
			outfile = '{0}/{1}.testlinks'.format(paths['CurrentTopDir'], classif)
			outfile_untransposed = '{0}/{1}.testlinks.untransposed'.format(paths['CurrentTopDir'], classif)
			g0fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g0')
			g1fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g1')
			g2fl = '{0}/{1}_{2}.summary'.format(paths['GENECONV_{0}_dir'.format(classif)],classif, 'g2')
			links = {}
			links_untransposed = {}
			if os.path.isfile(g0fl):
				# Convert GENECONV output to Circos links track
				links = geneconv2circoslinks(g0fl, paths['CurrentGFF'], outfile, append=False, output='return', linksdct=None)
				links_untransposed = geneconv2circoslinks(g0fl, paths['CurrentGFF'], outfile, append=False, output='return', linksdct=None, transposeLinks=False)
			if os.path.isfile(g1fl):
				# Convert GENECONV output to Circos links track
				links = geneconv2circoslinks(g1fl, paths['CurrentGFF'], outfile, append=True, output='return', linksdct=links)
				links_untransposed = geneconv2circoslinks(g1fl, paths['CurrentGFF'], outfile, append=True, output='return', linksdct=links_untransposed, transposeLinks=False)
			if os.path.isfile(g2fl):
				# Convert GENECONV output to Circos links track
				links = geneconv2circoslinks(g2fl, paths['CurrentGFF'], outfile, append=True, output='return', linksdct=links)
				links_untransposed = geneconv2circoslinks(g2fl, paths['CurrentGFF'], outfile, append=True, output='return', linksdct=links_untransposed, transposeLinks=False)
			append2logfile(paths['output_top_dir'], mainlogfile, 'Created links tracks for Circos from intra-cluster inter-element GENECONV output')
			# Modify geneconv2circoslinks to include an option to return the links information instead of writing to file.
			# Then use the returned infromation in the cluster loop to write a links track just for the cluster.
			#
			#
			# Modify links and then do four runs of the following (g0,g1,g2,all_3)
			#
			# Do not make small element figures if the link file is empty or nonexistent
			#
			#
			G_incl = [ ['g0'], ['g1'], ['g2'], ['g0', 'g1', 'g2'] ]

			for G in G_incl:
				if MCLCLUST:
					clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
				elif WICKERCLUST:
					clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(WickerParams['pId'], WickerParams['percAln'], WickerParams['minLen'], classif)]

				clusters = [ clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n') ]

				totallengths = {}
				element_coords = {}
				for i in range(len(clusters)):
					if len(clusters[i]) < 2:
						continue
					clusterscafs = set()
					outputlinks = []
					outputlinks_untransposed = []
					GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(paths['CurrentTopDir'], classif, i)
					if os.path.isfile(GFFoutPth):
						os.remove(GFFoutPth)
					highlights_ltrs_fl = '{0}/{1}.cluster_{2}.LTR_highlights.track'.format(paths['CurrentTopDir'], classif, i)
					if os.path.isfile(highlights_ltrs_fl):
						os.remove(highlights_ltrs_fl)
					with open(paths['CurrentGFF']) as gffFl:
						for line in gffFl:
							if '\tLTR_retrotransposon\t' in line:
								gffLine = GFF3_line(line)
								start = int(gffLine.start)
								end = int(gffLine.end)
								name = gffLine.attributes['ID']
								element_coords[name] = (start, end)
								el = gffLine.attributes['ID']
								# Only add links from elements in i
								if el not in clusters[i]:
									continue
								# Add link to output links
								for g in G:
									if g in links:
										if el in links[g]:
											outputlinks += links[g][el]
											outputlinks_untransposed += links_untransposed[g][el]
								scaf = gffLine.seqid
								if scaf not in clusterscafs:
									clusterscafs.add(scaf)
								with open(GFFoutPth, 'a') as GFFoutFl:
									GFFoutFl.write(line)
								gff2heatmapCallPacket = ([ '{0}/gff2circos-heatmap.py'.format(paths['scriptsDir']), '-gff', GFFoutPth, '-window', window, '-scafLens', scafLengthsFlPth ], '{0}.heatmap.track'.format(GFFoutPth), None, None)
								gff2tileCallPacket = ([ '{0}/gff2circos-tile.py'.format(paths['scriptsDir']), '-valueDef', 'LTR', '-gff', GFFoutPth ], '{0}.tile.track'.format(GFFoutPth), None, None)
								#gff2textLabelcallPacket = ([ '{0}/gff2circos-tile.py'.format(paths['scriptsDir']), '-valueDef', 'text', '-gff', GFFoutPth ], '{0}.tile.text.track'.format(GFFoutPth), None, None)
								tilecalls.append(gff2tileCallPacket)
								#tilecalls.append(gff2textLabelcallPacket)
								append2logfile(paths['output_top_dir'], mainlogfile, 'gff2circos-heatmap.py:\n{0}'.format(' '.join(gff2heatmapCallPacket[0])))
								append2logfile(paths['output_top_dir'], mainlogfile, 'gff2circos-tile.py:\n{0}'.format(' '.join(gff2tileCallPacket[0])))
								heatmapcalls.append(gff2heatmapCallPacket)
							# Write highlights track
							elif '\tlong_terminal_repeat\t' in line:
								GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(paths['CurrentTopDir'], classif, i)
								gffLine = GFF3_line(line)
								start = int(gffLine.start)
								end = int(gffLine.end)
								name = gffLine.attributes['Parent']
								newstart = start - element_coords[name][0] + 1
								newend = end - element_coords[name][0] + 1
								with open(highlights_ltrs_fl, 'a') as outFl:
									outFl.write('{0}\t{1}\t{2}\tfill_color=black\n'.format(name, newstart, newend))



								# Append to current LTR highlight track
								# Subtract the start and end of the LTR from the 
					with open('{0}/{1}.cluster_{2}.geneconv_{3}.links.track'.format(paths['CurrentTopDir'], classif, i, '_'.join(G)), 'w') as outFl:
						outFl.write('\n'.join(outputlinks))
					with open('{0}/{1}.cluster_{2}.geneconv_{3}.links_untransposed.track'.format(paths['CurrentTopDir'], classif, i, '_'.join(G)), 'w') as outFl:
						outFl.write('\n'.join(outputlinks_untransposed))

					# Write ideogram file for just scafs for this cluster
					totalseq = 0
					with open(allscafs, 'r') as inFl:
						ideoOut  = '{0}/{1}.cluster_{2}.seq.track'.format(paths['CurrentTopDir'], classif, i)
						with open(ideoOut, 'w') as outFl:
							for line in inFl:
								scaf = line.split()[2]
								if scaf in clusterscafs:
									contents = line.split()
									totalseq += int(contents[-2]) - int(contents[-3])
									outFl.write(line)
					totallengths[i] = totalseq
	#chr - Sacu_v1.1_s0011	11	0	2262239	greys-6-seq-4

					append2logfile(paths['output_top_dir'], mainlogfile, 'Created GFF files for each classification for converting to Circos heatmap tracks.')
						
				chunk_size = ceil(len(heatmapcalls)/procs)
				with Pool(processes=procs) as p:
					p.map(makecallMultiprocessing, heatmapcalls, chunksize=chunk_size)
				p.join()
				append2logfile(paths['output_top_dir'], mainlogfile, 'Converted GFFs to heatmap tracks for Circos.')

				chunk_size = ceil(len(tilecalls)/procs)
				with Pool(processes=procs) as p:
					p.map(makecallMultiprocessing, tilecalls, chunksize=chunk_size)
				p.join()
				append2logfile(paths['output_top_dir'], mainlogfile, 'Converted GFFs to tile tracks for Circos.')
				
				# Circos plot 1: ideograms are scaffolds
				for i in range(len(clusters)):
					if len(clusters[i]) < 2:
						continue
					GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(paths['CurrentTopDir'], classif, i)
					tilefl = '{0}.tile.track'.format(GFFoutPth)
					textfl = '{0}.tile.text.track'.format(GFFoutPth)
					linksfl = '{0}/{1}.cluster_{2}.geneconv_{3}.links.track'.format(paths['CurrentTopDir'], classif, i, '_'.join(G))
					seqfl = '{0}/{1}.cluster_{2}.seq.track'.format(paths['CurrentTopDir'], classif, i)
					if os.stat(linksfl).st_size == 0 and not G == ['g0', 'g1', 'g2']: # If no links for this cluster, don't draw a Circos plot for the elements without scaffolds.unless this is the composite with all gscale values, to ensure at least one plot is drawn for each cluster.
						continue
					#if os.path.isfile(tilefl) and os.path.isfile(linksfl) and os.path.isfile(seqfl) and os.path.isfile(textfl):
					if os.path.isfile(tilefl) and os.path.isfile(linksfl) and os.path.isfile(seqfl):
						# Files exist. copy and run Circos.

						# Plot with scaffolds
						circosdir = '{0}/circos.{1}.cluster_{2}.geneconv_{3}'.format(paths['CurrentTopDir'], classif, i, '_'.join(G))
						if not os.path.exists(circosdir):
							copytree('{0}/circos'.format(paths['scriptsDir']), circosdir) # copy circos conf files and dir structure
						newtilefl = '{0}/data/{1}'.format(circosdir, tilefl.split('/')[-1])
						if not os.path.isfile(newtilefl):
							copyfile(tilefl, newtilefl)
						newlinksfl = '{0}/data/{1}'.format(circosdir, linksfl.split('/')[-1])
						if not os.path.isfile(newlinksfl):
							copyfile(linksfl, newlinksfl)
						newseqfl = '{0}/data/{1}'.format(circosdir, seqfl.split('/')[-1])
						if not os.path.isfile(newseqfl):
							copyfile(seqfl, newseqfl)
						newtextfl = '{0}/data/{1}'.format(circosdir, textfl.split('/')[-1])
						if not os.path.isfile(newtextfl):
							copyfile(tilefl, newtextfl)
						conffl = '{0}/etc/circos.conf'.format(circosdir)
						confbasename = conffl.split('/')[-1]
						tileblock = '''
<plot>
type	=	tile
thickness	=	30
file	=	data/{0}
color	=	dorange
r1	=	0.84r
r0	=	0.78r
</plot>
'''.format(newtilefl.split('/')[-1])
						glyphblock = '''
<plot>
type	=	scatter
glyph	=	circle
glyph_size = 60
file	=	data/{0}
color	=	vdorange
orientation = out
r1	=	0.80r
r0	=	0.80r
</plot>
'''.format(newtilefl.split('/')[-1])
						if totallengths[i] > 5000000:
							plotblock = glyphblock
						else:
							plotblock = tileblock
						circos_conf_str = '''<<include colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype   = data/{0}

chromosomes_units           = 1000000

<plots>
<plot>
type             = text
color            = black
file             = data/{1}

r0 = 0.84r
r1 = 0.99r

show_links     = no
link_dims      = 0p,10p,60p,10p,0p
link_thickness = 10p
link_color     = red

label_snuggle         = yes
max_snuggle_distance  = 1r
snuggle_tolerance     = 0.40r
snuggle_sampling      = 2
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 2p
snuggle_refine        = yes
label_rotate = yes
label_size   = 100p
label_font   = condensed

padding  = 1p
rpadding = 1p

</plot>

{2}


</plots>

<links>

radius = 0.78r
crest  = 1
ribbon           = yes
flat             = yes
stroke_color     = vdgrey
stroke_thickness = 2
color            = grey_a3

bezier_radius        = 0r
bezier_radius_purity = 0.5

<link>
file       = data/{3}
</link>

</links>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim'''.format(newseqfl.split('/')[-1], newtextfl.split('/')[-1],  plotblock, newlinksfl.split('/')[-1])

						with open(conffl, 'w') as outFl:
							outFl.write(circos_conf_str)
						
						confbasename = conffl.split('/')[-1]
						imagesize = totallengths[i]/10
						if imagesize > 6000:
							imagesize = 6000
						conffl = '{0}/etc/image.generic.conf'.format(circosdir)
						circosimageconfstr = '''dir   = . 
file  = circos.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius = {0}

# by default angle=0 is at 3 o'clock position
angle_offset      = -96

#angle_orientation = counterclockwise

auto_alpha_colors = yes
auto_alpha_steps  = 5'''.format(imagesize)
						with open(conffl, 'w') as outFl:
							outFl.write(circosimageconfstr)
						#circos_call = [executables['circos'], '-silent', '-conf', confbasename]
						circos_call = [executables['perl'], executables['circos']]
						circoscalls.append([circosdir, circos_call, classif, i, '{0}/plots.scaffolds'.format(paths['CurrentTopDir']), G])
					

				MakeDir('plotdir', '{0}/plots.scaffolds'.format(paths['CurrentTopDir']))
				chunk_size = ceil(len(circoscalls)/procs)
				with Pool(processes=procs) as p:
					p.map(circosMultiprocessing, circoscalls, chunksize=chunk_size)
				p.join()
				append2logfile(paths['output_top_dir'], mainlogfile, 'Made Circos plots.')


				# Circos plot 2: ideograms are elements
				for i in range(len(clusters)):
					if len(clusters[i]) < 2:
						continue
					GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(paths['CurrentTopDir'], classif, i)
					tilefl = '{0}.tile.track'.format(GFFoutPth)
					seqfl = '{0}/{1}.cluster_{2}.seq.track'.format(paths['CurrentTopDir'], classif, i)
					highlights_ltrs_fl = '{0}/{1}.cluster_{2}.LTR_highlights.track'.format(paths['CurrentTopDir'], classif, i)
					links_untransposedfl = '{0}/{1}.cluster_{2}.geneconv_{3}.links_untransposed.track'.format(paths['CurrentTopDir'], classif, i, '_'.join(G))
					if os.stat(links_untransposedfl).st_size == 0: # If no links for this cluster, don't draw a Circos plot for the elements without scaffolds.
						continue
					if os.path.isfile(highlights_ltrs_fl) and os.path.isfile(tilefl) and os.path.isfile(links_untransposedfl) and os.path.isfile(seqfl):
						# Files exist. copy and run Circos.
						circosdir = '{0}/circos.{1}.cluster_{2}.geneconv_{3}.justelements'.format(paths['CurrentTopDir'], classif, i, '_'.join(G))
						if not os.path.exists(circosdir):
							copytree('{0}/circos'.format(paths['scriptsDir']), circosdir) # copy circos conf files and dir structure

						totallengthsLTRs = 0
						newseqfl = '{0}/data/{1}'.format(circosdir, '{0}.seq.track'.format('.'.join(tilefl.split('/')[-1].split('.')[:-2])))
						newhlfl = '{0}/data/{1}'.format(circosdir, '{0}.LTR_highlights.track'.format('.'.join(highlights_ltrs_fl.split('/')[-1].split('.')[:-2])))
						# Copy hl fl to circos dir
						copyfile(highlights_ltrs_fl, newhlfl)
						if os.path.isfile(newseqfl):
							os.remove(newseqfl)
						# Convert tile file to ideogram track
						with open(tilefl, 'r') as inFl:
							with open(newseqfl, 'w') as outFl:
								for line in inFl:
									scaf, start, end, val = line.strip().split()
									length = int(end) - int(start) + 1
									totallengthsLTRs += length
									color = 'dorange'
									outline = 'chr - {0} {1} 0 {2} {3}\n'.format('LTR_retrotransposon{0}'.format(val), val, length, color)
									outFl.write(outline)
								
								#==> Copia.cluster_0.gff.tile.track <==
								#Sacu_v1.1_s0016	1385103	1391765 123
								#
								#888> Copia.cluster_0.seq.track <==
								#chr - Sacu_v1.1_s0001	1	0	4132625	greys-6-seq-4

						newlinksuntransposedfl = '{0}/data/{1}'.format(circosdir, links_untransposedfl.split('/')[-1])
						copyfile(links_untransposedfl, newlinksuntransposedfl)

						# Don't use this ideogram track
						#newseqfl = '{0}/data/{1}'.format(circosdir, seqfl.split('/')[-1])
						#if not os.path.isfile(newseqfl):
						#	copyfile(seqfl, newseqfl)
						conffl = '{0}/etc/circos.conf'.format(circosdir)
						circos_conf_str = '''<<include colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype   = data/{0}

chromosomes_units           = 1000000

<highlights>
 <highlight>
 file       = data/{1}
 ideogram   = yes
 color = black
 </highlight>
</highlights>


<links>
radius = 0.999r
#radius = 1r
crest  = 1
ribbon           = yes
flat             = yes
#stroke_color     = vdgrey
stroke_thickness = 2
color            = grey_a3

bezier_radius        = 0r
bezier_radius_purity = 0.5

<link>
file       = data/{2}
</link>

</links>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim'''.format(newseqfl.split('/')[-1], newhlfl.split('/')[-1], newlinksuntransposedfl.split('/')[-1])
						with open(conffl, 'w') as outFl:
							outFl.write(circos_conf_str)
						
						confbasename = conffl.split('/')[-1]
						imagesize = totallengthsLTRs/10
						if imagesize > 8000:
							imagesize = 8000
						if imagesize < 1000:
							imagesize = 1000
						conffl = '{0}/etc/image.generic.conf'.format(circosdir)
						circosimageconfstr = '''dir   = . 
file  = circos.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius = {0}

# by default angle=0 is at 3 o'clock position
angle_offset      = -96

#angle_orientation = counterclockwise

auto_alpha_colors = yes
auto_alpha_steps  = 5'''.format(imagesize)
						with open(conffl, 'w') as outFl:
							outFl.write(circosimageconfstr)
						#circos_call = [executables['circos'], '-silent', '-conf', confbasename]
						circos_call = [executables['perl'], executables['circos']]
						circoscalls.append([circosdir, circos_call, classif, i, '{0}/plots.elements'.format(paths['CurrentTopDir']), G])
					
				MakeDir('plotdir', '{0}/plots.elements'.format(paths['CurrentTopDir']))
				chunk_size = ceil(len(circoscalls)/procs)
				with Pool(processes=procs) as p:
					p.map(circosMultiprocessing, circoscalls, chunksize=chunk_size)
				p.join()
				append2logfile(paths['output_top_dir'], mainlogfile, 'Made Circos plots.')




def shortHelp():

	print('''
	  Usage:
	  ------------
	  phaltr -f|--fasta <path> [--logfile <path>] [-p|--procs <int>] [-c|--clean]
	  [-o|--output_dir <path>] [-lh|--ltrharvest] [--minlenltr <int>] [--maxlenltr <int>]
	  [--mindistltr <int>] [--maxdistltr <int>] [--similar <int|float>] [--vic <int>]
	  [--mintsd <int>] [--maxtsd <int>] [--xdrop <int>] [--mat <int>] [--mis <int>]
	  [--ins <int>] [--del <int>] [-ld|--ltrdigest] [--ltrdigest_hmms <path>] [--classify]
	  [--classify_dfam] [--classify_repbase] [--nhmmer_reporting_evalue <int|float>]
	  [--nhmmer_inclusion_evalue <int|float>] [--repbase_tblastx_evalue <int|float>] 
	  [--keep_conflicting_classificaitons] [--keep_no_classifications] [--min_clust_size <int>]
	  [--wicker] [--wicker_pId <int|float>] [--wicker_minLen <int>] [--wicker_pAln <int|float>]
	  [--wicker_no_ltrs] [--wicker_no_internals] [--mcl] [-I <int|float>] [--nosmalls]
	  [--geneconvltrs] [--geneconv_g <[g0[,g1[,g2]]]>] [--ltrdivergence] [--remove_GC_from_modeltest_aln]
	  [--modeltest_criterion <str>] [--gc_ltr_div_scaling <int>] [--maxiterate_small_clusters <int>]
	  [--maxiterate_large_clusters <int>] [--min_clustsize_for_faster_aln <int>] [--mafft_retree <int>]
	  [--geneconvclusters] [--DTT] [--phylo] [--bootstrap_reps] [--bpflank <int>]
	  [--flank_evalue <int|float>] [--flank_pId <int|float>][--flank_plencutoff <int|float>]
	  [--min_orf_len <int>]
	  ''', file=sys.stderr)

def help2():
	print('''
	  Option			    ArgType	       Default
	----------------------------------------------------------------
	--fasta				    <path>		NONE
	--procs				    <int>		1
	--clean				    BINARY		OFF
	--output_dir			    <path>		LTRan_output
	--ltrharvest			    BINARY		OFF
	--minlenltr			    <int>		100
	--maxlenltr			    <int>		1000
	--mindistltr			    <int>		1000
	--maxdistltr			    <int>		15000
	--similar			    <int|float>		85.0
	--vic				    <int>		60
	--mintsd			    <int>		4
	--maxtsd			    <int>		20
	--xdrop				    <int>		5
	--mat				    <int>		2
	--mis				    <int>		-2
	--ins				    <int>		-3
	--del				    <int>		-3
	--ltrdigest			    BINARY		OFF
	--ltrdigest_hmms		    <path>		{0}
	--classify			    BINARY		OFF
	--classify_dfam			    BINARY		OFF
	--classify_repbase		    BINARY		OFF
	--nhmmer_reporting_evalue	    <int|float>		10
	--nhmmer_inclusion_evalue	    <int|float>		1e-2
	--repbase_tblastx_evalue	    <int|float>		1e-5
	--keep_conflicting_classificaitons  BINARY		OFF
	--keep_no_classifications	    BINARY		OFF
	--min_clust_size		    <int>		7
	--wicker			    BINARY		OFF
	--wicker_pId			    <int|float>		80
	--wicker_minLen			    <int>		80
	--wicker_pAln			    <int|float>		80
	--wicker_no_internals		    BINARY		OFF
	--wicker_no_ltrs		    BINARY		OFF
	--mcl				    BINARY		MCL
	-I				    <int|float>		6
	--nosmalls			    BINARY		OFF
	--geneconvltrs			    BINARY		OFF
	--geneconv_g			    <str>		g0,g1,g2
	--ltrdivergence			    BINARY		OFF
	--remove_GC_from_modeltest_aln	    BINARY		OFF
	--modeltest_criterion		    <str>		BIC
	--gc_ltr_div_scaling		    <int>		1
	--maxiterate_small_clusters	    <int>		30
	--maxiterate_large_clusters	    <int>		3
	--min_clustsize_for_faster_aln	    <int>		40
	--mafft_retree			    <int>		2
	--geneconvclusters		    BINARY		OFF
	--DTT				    BINARY		OFF
	--phylo				    BINARY		OFF
	--bootstrap_reps		    <int>		100
	--bpflank			    <int>		500
	--flank_evalue			    <int|float>		1e-5
	--flank_pId			    <int|float>		70
	--flank_plencutoff		    <int|float>		70
	--min_orf_len			    <int>		300

	'''.format('{0}/LTRdigest_HMMs/hmms'.format(paths['selfDir']), file=sys.stderr))

def help():

	print('''

	  Usage:
	  ------------
	  phaltr [options] -fasta <input.fasta>
	  
	  Description:
	  ------------
	  The main options would be:
	  phaltr --fasta fasta.fa \\
	  	 --procs 40 \\
		 --ltrharvest \\
		 --ltrdigest \\				Use: phaltr -h
		 --classify \\				   to see the defaults
		 --mcl \\
		 --wicker \\
		 --geneconvltrs \\
		 --geneconvclusters \\
		 --ltrdivergence \\
		 --DTT

	  
	  ------------------------------
	  Global Options:
	  ------------------------------
	  -f | --fasta			<path>	Sequences to analyze. Mandatory.
	  -p | --procs			<int>	Number of processors (default 1)
	  -c | --clean				Remove all temporary files (NOT IMPLEMENTED YET)
	  -o | --output_dir		<path>	Output directory. Default is "LTRan_output
	  --logfile			<path>  Path to where log file is written (default <output_dir>/log.txt)
	  -h					Defaults
	  -help					Long help
	  

	  -------------------------
	  Program-specific Options:
	  -------------------------

	  LTRharvest
	  ----------
	  -lh | --ltrharvest		Run LTRharvest on file given by --fasta (default ON)
  	  --minlenltr	<int>		minimum length allowed for LTRs for element calling (default 100 bp)
  	  --maxlenltr	<int>		maximum length allowed for LTRs for element caling (default 1000 bp)
  	  --mindistltr	<int>		minimum distance allowed between LTRs for element calling (default 1000 bp)
  	  --maxdistltr	<int>		maximum distance allowed between LTRs for element calling (default 15000 bp)
  	  --similar	<int|float>	minimum % similarity for LTR calling (default 85.0)
  	  --vic		<int>		# of nucleotides to left and right to search for TSDs (default 60)
  	  --mintsd	<int>		minimum length allowed for TSDs (use with --maxtsd) (default 4)
  	  --maxtsd	<int>		maximum length allowed for TSDs (use with --mintsd) (default 20)
  	  --xdrop	<int>		xdropbelow score for extension-alignment (default 5)
  	  --mat		<int>		matchscore score for extension-alignment (default 2)
  	  --mis		<int>		mismatchscore score for extension-alignment (default -2)
  	  --ins		<int>		insertionscore for extension-alignment (default -3)
  	  --del		<int>		deletionscore for extension-alignment (default -3)



	  --min_orf_len			<int>	(default 300)


	  LTRdigest
	  ---------
	  -ld | --ltrdigest			Run LTRdigest on file given by --fasta and --ltrharvest results (GFF) (default ON)
	  --ltrdigest_hmms		<path>	Path to a file with one or more protein profile HMMs for LTRdigest (HMMER)
	  					(default {0})


	  Classification of LTR RTs to superfamily using homology to annotated sequences in Repbase and/or Dfam
	  -----------------------------------------------------------------------------------------------------
	  --classify				Do both Dfam and Repbase classifications
	  --classify_dfam			Run hmmsearch on Dfam database (default ON)
	  --classify_repbase			Run tblastx on Repbase database (default ON)
	  --keep_conflicting_classifications		If an element has two classifications that disagree (i.e. Repbase and Dfam) and one classification
	  						is an LTR RT hit and one is a non-LTR RT hit, and this flag is set, the element will be kept.
							(default OFF; the element is discarded as a false positive)
	  --keep_no_classification			If an element does not have evidence of homology to any sequence in Repbase or Dfam, keep it.
	  						(default OFF; elements without homology to LTR RT in a database are discarded as false positives)
	  --repbase_tblastx_evalue	<int|float>	Max allowed E-value for tblastx hits (default 1e-5)
	  --nhmmer_reporting_evalue	<int|float> 	(default 10)
	  --nhmmer_inclusion_evalue	<int|float>	(default 1e-2)


	  Clustering
	  -----------
	  --min_clust_size			Minimum allowed cluster size. Clusters with < minclustsize elements get assembled together
	  					but no model testing, gene conversion, or outgroup/DTT analysis is done.
  	  --mcl					Cluster using MCL
	  -I					Inflation/granularity parameter for MCL (default 6)
  	  --wicker				Cluster using '80-80-80' rule, or custom values specified below.
  	  --wicker_pId				Minimum % ID in pairwise alignment between any two elements (default 80)
  	  --wicker_minLen			Minimum alignment length to considered. (default 80)
  	  --wicker_pAln				Minimum percentage of whole sequence (LTRs or internal region) for alignment
	  					to be considered. (default 80)
	  --wicker_no_internals			(default OFF)
	  --wicker_no_ltrs			(default OFF)


	  MAFFT (for cluster, not LTR alignment)
	  -----------------------------
	  --maxiterate_small_clusters		<int>	Max number of iterations for MAFFT algorithm. 1 is fastest; greater numbers will improve the
	  						alignment (default 30)
	  --maxiterate_large_clusters		<int>	Number of iterations for MAFFT algorithm for large clusters (default 3)
	  --min_clustsize_for_faster_aln	<int>	Clusters smaller than this will be aligned using the settings from --maxiterate_small_clusters,
	  						while those larger will get aligned using the settings from --maxiterate_large_clusters (default 40)
	  --retree				<int>	Guide tree is built <int> times in the progressive stage. Valid with 6mer distance. (default 2)
	  --nosmalls					Do not combine and assemble clusters smaller than --min_clust_size (see clustering options)

	  GENECONV
	  --------
	  --geneconv_g	<comma-sep-list>	A comma-separated list of the values for the 'g' parameter of GENECONV. Valid values are g0, g1, and g2.
	  					(default g0,g1,g2)
	  --geneconvltrs
	  --geneconvclusters


	  LTR divergence estimation
	  -------------------------
	  --ltrdivergence			Find statistially best supported (BIC) substitution model for each cluster (default ON)
	  					and estimate substitutions per site between LTRs for each element. 
	  --modeltest_criterion		<str>	AIC, AICc, or BIC. (default BIC)
						MULTIPLE CHOICES NOT YET IMPLEMENTED
	  --gc_ltr_div_scaling		<int>	For reporting scaled divergence estimates to account for effects of gene conversion, if observed. (default 1)
						MULTIPLE CHOICES NOT YET IMPLEMENTED
	  					1. divC = div * len(aln)/(len(aln)-len(gc_trac))
	  --default_model		<str>	When model testing is not possible (i.e. cluster is too small) HKY85 or JC,
	  					if HKY85 is not possible due to dinucleotide LTRs.
	  					MULTIPLE CHOICES NOT YET IMPLEMENTED
	  --remove_GC_from_modeltest_aln	Remove elements with suspected gene conversion tracts.


	  Finding pairs of elements within clusters that have homologous flanking regions
	  -------------------------------------------------------------------------------
	  --rmhomoflank				Remove one of each pair of elements within each alignment (and therefore, each tree).
	  					(default OFF; Fixed ON when using --DTT)
	  --bpflank		<int>		Number of bases on either side of each element to search for homology. (default 500 bp)
	  --flank_evalue	<int|float>	E-value ceiling for considering blastn hits as evidence of homology. (default 1e-5)
	  --flank_pId		<int|float>	Minimum percent identity in blastn alignment to consider hit as evidence of homology. (default 70)
	  --flank_plencutoff	<int|float>	Minimum percentage of flanking region required to participate in alignment to consider
	  					blastn hit as evidence of homology. (default 70)


	  Phylogenetic analysis
	  ---------------------
	  --phylo				##### not implemented yet
	  --nosmalls				Do not combine and perform phylogentic analyses on clusters smaller than --min_clust_size.
	  --DTT					Turns on --rmhomoflank, --convert_to_ultrametric, and --auto_outgroup. Generates and attempts to run
	  					Rscript that generates a DTT (LTT) plot for each cluster for which rooting is possible.
	  --bootstrap_reps		<int>	Number of replicates to perform when bootstrapping (default 100)
	  --convert_to_ultrametric		Convert trees to ultrametric using PATHd8. (default OFF; ON when using --DTT)
	  --auto_outgroup			Pick an outgroup automatically:
						  The outgroup shall be a random element from cluster k where cluster k is the largest of the clusters
						  that is not j if j is the first cluster then the next smallest cluster is k if there is no other
						  cluster, no outgroup is used.

		'''.format('{0}/LTRdigest_HMMs/hmms'.format(paths['selfDir']), file=sys.stderr))
args=sys.argv

# get executable paths from CONFIG file, which should be in the same directory as this script
executables = {}
commentPattern = re.compile('#.*$')
with open('{0}/CONFIG'.format(os.path.dirname(os.path.realpath(__file__)))) as config_file:
	paths = [ re.sub(commentPattern, '', line) for line in config_file.read().strip().split('\n') ]
	for path in paths:
		if not path == '':
			p = path.split('=')
			executables[p[0]] = p[1]
	
filenames = {}
paths = {}
paths_toClean = {}
params = {}
# paths may contain:
# inputFasta			Input fasta
# inputFastaSuffixArray		Basename for suffixerator output
# LTRharvestGFF			GFF3 output from LTRharvest
# LTRdigestOutputPrefix		Prefix for LTRdigest output files
# LTRdigestHMMs			Directory for individual Pfam HMMs for LTR retrotransposon domains

paths['selfDir'] = '/'.join(os.path.realpath(__file__).split('/')[:-1])
paths['scriptsDir'] = '{0}/scripts'.format(paths['selfDir'])

if '-h' in args or '--h' in args:
	help2()
	sys.exit(0)

if '-help' in args or '--help' in args:
	help()
	sys.exit(0)
if len(args) < 3:
	shortHelp()
	sys.exit(0)

if '--fasta' in args:
	paths['inputFasta'] = args[args.index('--fasta') + 1]
elif '-f' in args:
	paths['inputFasta'] = args[args.index('-f') + 1]
else:
	help()
	print('''

		MUST SPECIFY INPUT FASTA WITH -f or --fasta

		''', file=sys.stderr)
	sys.exit(0)

filenames['inputFasta'] = paths['inputFasta'].split('/')[-1]

if '--nosmalls' in args:
	SMALLS = False
else:
	SMALLS = True
if '--del' in args:
	ltrharvest_del = int(args[args.index('--del')+1])
else:
	ltrharvest_del = -3
if '--ins' in args:
	ltrharvest_ins = int(args[args.index('--ins')+1])
else:
	ltrharvest_ins = -3
if '--mis' in args:
	ltrharvest_mis = int(args[args.index('--mis')+1])
else:
	ltrharvest_mis = -2
if '--mat' in args:
	ltrharvest_mat = int(args[args.index('--mat')+1])
else:
	ltrharvest_mat = 2
if '--xdrop' in args:
	ltrharvest_xdrop = int(args[args.index('--xdrop')+1])
else:
	ltrharvest_xdrop = 5 
if '--minlenltr' in args:
	ltrharvest_minlenltr = int(args[args.index('--minlenltr')+1])
else:
	ltrharvest_minlenltr = 100
if '--maxlenltr' in args:
	ltrharvest_maxlenltr = int(args[args.index('--maxlenltr')+1])
else:
	ltrharvest_maxlenltr = 1000
if '--mindistltr' in args:
	ltrharvest_mindistltr = int(args[args.index('--mindistltr')+1])
else:
	ltrharvest_mindistltr = 1000
if '--maxdistltr' in args:
	ltrharvest_maxdistltr = int(args[args.index('--maxdistltr')+1])
else:
	ltrharvest_maxdistltr = 15000
if '--similar' in args:
	ltrharvest_similar = float(args[args.index('--similar')+1])
else:
	ltrharvest_similar = 85
if '--vic' in args:
	ltrharvest_vic = int(args[args.index('--vic')+1])
else:
	ltrharvest_vic = 60
if '--mintsd' in args:
	ltrharvest_mintsd = int(args[args.index('--mintsd')+1])
else:
	ltrharvest_mintsd = 4
if '--maxtsd' in args:
	ltrharvest_maxtsd = int(args[args.index('--maxtsd')+1])
else:
	ltrharvest_maxtsd = 20
if '--logfile' in args:
	mainlogfile = args[args.index('--logfile')+1]
else:
	mainlogfile = 'log.txt'
if '--procs' in args:
	procs = int(args[args.index('--procs') + 1])
if '-p' in args:
	procs = int(args[args.index('-p') + 1])
else:
	procs = 1
if '--output_dir' in args:
	MakeDir('output_top_dir', args[args.index('--output_dir') + 1])
if '-o' in args:
	MakeDir('output_top_dir', args[args.index('-o') + 1])
else:
	MakeDir('output_top_dir', 'LTRan_output')

if '--ltrharvest' in args or '-lh' in args: # Turn on LTRharvest for file given by --fasta
	LTRHARVEST = True
else:
	LTRHARVEST = False

if '--ltrdigest' in args or '-ld' in args: # Turn on LTRdigest for LTRharvest results
	LTRDIGEST = True
else:
	LTRDIGEST = False

if '--ltrdigest_hmms' in args: # Check for user-supplied location of HMMs for LTRdigest, set default otherwise (comes with package)
	paths['LTRdigestHMMs'] = args[args.index('--ltrdigest_hmms') + 1]
else:
	paths['LTRdigestHMMs'] = '{0}/LTRdigest_HMMs/hmms'.format(paths['selfDir'])
if '--min_orf_len' in args:
	min_orf_len = int(args[args.index('--min_orf_len')+1])
else:
	min_orf_len = 300

if '--classify_dfam' or '--classify' in args: # Classification Parameters
	CLASSIFYDFAM = True
else:
	CLASSIFYDFAM = False

if '--repbase_tblastx_evalue' in args:
	repbase_tblastx_evalue = float(args[args.index('--repbase_tblastx_evalue')+1])
else:
	repbase_tblastx_evalue = 1e-5

if '--classify_repbase' in args or '--classify' in args:
	CLASSIFYREPBASE = True
else:
	CLASSIFYREPBASE = False
	os.environ['BLASTDB'] = paths['FastaOutputDir']
if '--keep_conflicting_classifications':
	KEEPCONFLICTS=True
else:
	KEEPCONFLICTS=False
if '--keep_no_classification':
	KEEPNOCLASSIFICATION=True,
else:
	KEEPNOCLASSIFICATION=False,
if '--nhmmer_reporting_evalue' in args:
	nhmmer_reporting_evalue = float(args[args.index('--nhmmer_reporting_evalue')+1])
else:
	nhmmer_reporting_evalue = 10
if '--nhmmer_inclusion_evalue' in args:
	nhmmer_inclusion_evalue = float(args[args.index('--nhmmer_inclusion_evalue')+1])
else:
	nhmmer_inclusion_evalue = 1e-5
if '--wicker' in args:
	WICKER = True
else:
	WICKER = False

if '--mcl' in args:
	USEMCL = True
	if '-I' in args:
		MCL_I = args[args.index('-I')+1]
	else:
		MCL_I = '6'

else:
	USEMCL = False

if '--min_clust_size' in args:
	MinClustSize = int(args[args.index('--min_clust_size')+1])
else:
	MinClustSize = 7

if '--ltrdivergence' in args:
	LTRDIVERGENCE = True
else:
	LTRDIVERGENCE = False

if '--geneconvltrs' in args:
	GENECONVLTRS = True
else:
	GENECONVLTRS = False

if '--geneconvclusters' in args:
	GENECONVCLUSTERS = True
else:
	GENECONVCLUSTERS = False
if '--geneconv_g' in args:
	gcparams = args[args.index('--geneconv_g')+1].split(',')
	if 'g0' in gcparams:
		GENECONV_G0 = True
	else:
		GENECONV_G0 = False
	if 'g1' in gcparams:
		GENECONV_G1 = True
	else:
		GENECONV_G1 = False
	if 'g2' in gcparams:
		GENECONV_G2 = True
	else:
		GENECONV_G2 = False
else:
	GENECONV_G0 = True
	GENECONV_G1 = True
	GENECONV_G2 = True

if '--estimate_divergence' in args:
	DIVERGENCE = True
else:
	DIVERGENCE = False

if '--wicker_pId' in args:
	wicker_pId = int(args[args.index('--wicker_pId')+1])
else:
	wicker_pId = 80
if '--wicker_pAln' in args:
	wicker_pAln = int(args[args.index('--wicker_pAln')+1])
else:
	wicker_pAln = 80

if '--wicker_minLen' in args:
	wicker_minLen = int(args[args.index('--wicker_minLen')+1])
else:
	wicker_minLen = 80
if '--remove_GC_from_modeltest_aln' in args:
	remove_GC_from_modeltest_aln = True
else:
	remove_GC_from_modeltest_aln = False
if '--bootstrap_reps' in args:
	bootstrap_reps = int(args[args.index('--bootstrap_reps')+1])
else:
	bootstrap_reps = 100
if '--rmhomoflank' in args:
	RMHOMOFLANK = True
else:
	RMHOMOFLANK = False
if '--convert_to_ultrametric' in args:
	ULTRAMETRIC = True
else:
	ULTRAMETRIC = False
if '--auto_outgroup' in args:
	AUTO_OUTGROUP = True
else:
	AUTO_OUTGROUP = False
if '--DTT' in args:
	AUTO_OUTGROUP = True
	RMHOMOFLANK = True
	DTT = True
	ULTRAMETRIC = True
else:
	AUTO_OUTGROUP = False
	RMHOMOFLANK = False
	DTT = False
	ULTRAMETRIC = False

# MAFFT parameters
if '--maxiterate_small_clusters' in args:
	mafft_small_maxiterate = int(args[args.index('--maxiterate_small_clusters')+1])
else:
	mafft_small_maxiterate = 30
if '--maxiterate_large_clusters' in args:
	mafft_large_maxiterate = int(args[args.index('--maxiterate_large_clusters')+1])
else:
	mafft_large_maxiterate = 3
if '--min_clustsize_for_faster_aln' in args:
	mafft_large_minclustsize = int(args[args.index('--min_clustsize_for_faster_aln')+1])
else:
	mafft_large_minclustsize = 40
if '--retree' in args:
	mafft_retree = int(args[args.index('--retree')+1])
else:
	mafft_retree = 2
if '--bpflank' in args:
	bpflank = int(args[args.index('--bpflank')+1])
else:
	bpflank = 500
if '--flank_evalue' in args:
	flank_evalue = int(args[args.index('--flank_evalue')+1])
else:
	flank_evalue = 1e-5
if '--flank_pId' in args:
	flank_pId = float(args[args.index('--flank_pId')+1])
else:
	flank_pId = 70.0
if '--flank_plencutoff' in args:
	flank_plencutoff = float(args[args.index('--flank_plencutoff')+1])
else:
	flank_plencutoff = 70.0
if '--wicker_no_ltrs' in args:
	wicker_use_ltrs = False
else:
	wicker_use_ltrs = True
if '--wicker_no_internals' in args:
	wicker_use_internal = False
else:
	wicker_use_internal = True
	
paths['RepbaseShortNames'] = '{0}/RepeatDatabases/Repbase/Repbase.annotations.LTR.names.cleaned.map'.format(paths['selfDir'])
paths['DfamShortNames'] = '{0}/RepeatDatabases/Dfam/Dfam.annotations.LTR.names.cleaned.map'.format(paths['selfDir'])
MakeDir('FastaOutputDir', '{0}/FASTA_output'.format(paths['output_top_dir']))
MakeDir('GFFOutputDir', '{0}/GFF_output'.format(paths['output_top_dir']))
paths['CurrentGFF'] = None # This path will have the path to the best GFF3 to use.
			    # Div > GC > TP > Ld > Lh
try:
	statusFlRead = open('{0}/status'.format(paths['output_top_dir']), 'r')
except:
	pass

def write2summary(text):
	with open('{0}/summary'.format(paths['output_top_dir']), 'a') as summaryFl:
		summaryFl.write('{0}\n'.format(text))

# Check for status file. if exists parse it and skip the sections that are done.
try:
	statusContents = statusFlRead.read()
	if not statusContents == '':
		for pair in statusContents.strip().split('\n'):
			pair = pair.split('\t')
			paths[pair[0]] = pair[1]
except:
	pass

if 'LTRdigestClassifiedNoFP' in paths:
	paths['CurrentGFF'] = paths['LTRdigestClassifiedNoFP']
elif 'GFFwithRepbaseClassification' in paths:
	paths['CurrentGFF'] = paths['GFFwithRepbaseClassification']
elif 'GFFwithDfamClassification' in paths:
	paths['CurrentGFF'] = paths['GFFwithDfamClassification']
elif 'WithORFsGFF' in paths:
	paths['CurrentGFF'] = paths['WithORFsGFF']
elif 'LTRdigestGFF' in paths:
	paths['CurrentGFF'] = paths['LTRdigestGFF']
elif 'LTRharvestGFF' in paths:
	paths['CurrentGFF'] = paths['LTRharvestGFF']

# These functions modify the global var paths
# They run various procedures and generate files
# They'll be skipped if they appear to have been done already and the requisite files for subsequent steps exist
# They'll also be skipped if the are not supposed to run for the requested procedure
# If they run they will append to the log


sys.setrecursionlimit(50000) # For WickerFam() recursive subroutine

#  I. Identification and classification of elements
#
#	1. Run LTRharvest
#
ltrharvest()    # Predict LTR retrotransposons using structural criteria
#		# Input: Sequences (FASTA)
#		# Output: LTRharvest GFF3
#
#	2. Run LTRdigest
#
ltrdigest()	# Identify parts of element internal regions with evidence of homology to LTR RT protein coding domains
#		# Input: Sequences (FASTA), LTRharvest GFF3, pHMMs
#		# Output: LTRdigest GFF3
#
AnnotateORFs(minLen=min_orf_len)
#
#	3. Classify elements to superfamily using homology evidence in Dfam and/or Repbase
#
classify_by_homology(KEEPCONFLICTS=KEEPCONFLICTS, KEEPNOCLASSIFICATION=KEEPNOCLASSIFICATION, repbase_tblastx_evalue=repbase_tblastx_evalue, nhmmer_reporting_evalue=nhmmer_reporting_evalue, nhmmer_inclusion_evalue=nhmmer_inclusion_evalue)  # Extract LTR_retrotransposon sequences for classification using homology
#			# Find evidence of homology to repeats in Dfam using HMMER. NEED TO CHANGE THIS SO REVERSE COMPLEMENT IS ALSO SEARCHED (nhmmsearch I think)
#			# Find evidence of homology to repeats in Repbase using tblastx
#			# Remove false positives from LTRdigest GFF3
#			# Input: Sequences (FASTA), DBs (Dfam & Repbase), LTRdigest or LTRharvest GFF3
#			# Output: LTR RT GFF3 with classifications and false positives removed (FPs have homology to non-LTR RTs in DBs or no homology)
#
#		Global variables containg superfamily assignment for each element
#
clusters_by_classif = shortClassif()
classifs_by_element = shortClassif(ElNames=True)
classifs = set(list(clusters_by_classif.keys()))
#
#		# Input: LTR RT structures (GFF3) and Sequences (FASTA)
#		# Output: List of elements with evidence of intra element LTR gene conversion (text table)
#
#Circos(window='1000000', plots='clusters', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
#if WICKER:
#	Circos(window='1000000', plots='clusters', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
#if USEMCL:
#	Circos(window='1000000', plots='clusters', I=MCL_I, clustering_method='MCL', WickerParams=None)
#sys.exit()
#
#  II. Clustering, divergence, gene conversion, and phylogenetic analysis
#
if WICKER:
	# 1. Perform clustering
	WickerFam(pId=wicker_pId, percAln=wicker_pAln, minLen=wicker_minLen, use_ltrs=wicker_use_ltrs, use_internal=wicker_use_internal)

	if GENECONVCLUSTERS or LTRDIVERGENCE:
		# 2. MSA for each cluster
		AutoAlign(I=None, part='entire', rmgeneconv=False, minClustSize=MinClustSize, align='clusters', rmhomologflank=False, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=SMALLS, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)

		# 3. GENECONV for each MSA
		if GENECONV_G0:
			geneconvClusters(trimal=True, g='/g0', force=False, clust=None, I=None, minClustSize=MinClustSize, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, combine_and_do_small_clusters=SMALLS)
		if GENECONV_G1:
			geneconvClusters(trimal=True, g='/g1', force=False, clust=None, I=None, minClustSize=MinClustSize, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, combine_and_do_small_clusters=SMALLS)
		if GENECONV_G2:
			geneconvClusters(trimal=True, g='/g2', force=False, clust=None, I=None, minClustSize=MinClustSize, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, combine_and_do_small_clusters=SMALLS)

		# 4. Modeltesting  for each cluster
		modeltest(iters=1, I=None, removegeneconv=remove_GC_from_modeltest_aln, part='entire', clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, minClustSize=MinClustSize, bpflank=bpflank, combine_and_do_small_clusters=SMALLS)

if USEMCL:
	# 1. Perform clustering
	MCL(I=MCL_I, minClustSize=MinClustSize)	# Run MCL. I is the inflation paramater that controls granularity. Default is 6. MCL docs recommend 1.4, 2, 4, and 6, and between 1.1 and 10 for most cases.

	if GENECONVCLUSTERS or LTRDIVERGENCE:
		# 2. MSA for each cluster
		AutoAlign(I=MCL_I, part='entire', rmgeneconv=False, minClustSize=MinClustSize, align='clusters', rmhomologflank=False, clustering_method='MCL', WickerParams=None, auto_outgroup=False, bpflank=bpflank, combine_and_do_small_clusters=SMALLS, flank_pId=flank_pId, flank_evalue=flank_evalue, flank_plencutoff=flank_plencutoff)

		# 3. GENECONV for each MSA
		if GENECONV_G0:
			geneconvClusters(trimal=True, g='/g0', force=False, clust=None, I=MCL_I, minClustSize=MinClustSize, clustering_method='MCL', WickerParams=None, combine_and_do_small_clusters=SMALLS)
		if GENECONV_G1:
			geneconvClusters(trimal=True, g='/g1', force=False, clust=None, I=MCL_I, minClustSize=MinClustSize, clustering_method='MCL', WickerParams=None, combine_and_do_small_clusters=SMALLS)
		if GENECONV_G2:
			geneconvClusters(trimal=True, g='/g2', force=False, clust=None, I=MCL_I, minClustSize=MinClustSize, clustering_method='MCL', WickerParams=None, combine_and_do_small_clusters=SMALLS)

		# 4. Modeltesting  for each cluster
		modeltest(iters=1, I=MCL_I, removegeneconv=remove_GC_from_modeltest_aln, part='entire', clustering_method='MCL', WickerParams=None, minClustSize=MinClustSize, bpflank=bpflank, combine_and_do_small_clusters=SMALLS)
#
#
#
#	GENECONV stringency settings
#
#	g = '/g0'	Most stringent - no mismatches in fragments
#	g = '/g1'	Most relaxed - mismatches in fragments
#	g = '/g2'	Medium relaxed
#
#
#  
if WICKER:
	Circos(window='1000000', plots='clusters', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
if USEMCL:
	Circos(window='1000000', plots='clusters', I=MCL_I, clustering_method='MCL', WickerParams=None)
#
#
if WICKER:
	align_ltrs(I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})	# Runs if need to use geneconvLTRs or estimate divergences
#
	if GENECONV_G0:
		geneconvLTRs(g='/g0', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})	# Identify LTR RTs with gene conversion. Most stringent - no mismatches in fragments
	if GENECONV_G1:
		geneconvLTRs(g='/g1', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})	# Identify LTR RTs with gene conversion. Most relaxed - mismatches in fragments
	if GENECONV_G2:
		geneconvLTRs(g='/g2', I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})	# Identify LTR RTs with gene conversion. Medium relaxed
if USEMCL:
	align_ltrs(I=MCL_I, clustering_method='MCL', WickerParams=None)	# Runs if need to use geneconvLTRs or estimate divergences
	#
	if GENECONV_G0:
		geneconvLTRs(g='/g0', I=MCL_I, clustering_method='MCL', WickerParams=None)	# Identify LTR RTs with gene conversion. Most stringent - no mismatches in fragments
	if GENECONV_G1:
		geneconvLTRs(g='/g1', I=MCL_I, clustering_method='MCL', WickerParams=None)	# Identify LTR RTs with gene conversion. Most relaxed - mismatches in fragments
	if GENECONV_G2:
		geneconvLTRs(g='/g2', I=MCL_I, clustering_method='MCL', WickerParams=None)	# Identify LTR RTs with gene conversion. Medium relaxed
#
#
if WICKER:
	ltr_divergence(I=None, clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
	phylo(removegeneconv=False, BOOTSTRAP=True, I=None, align='cluster', removehomologouspair=RMHOMOFLANK, part='entire', clustering_method='WickerFam', WickerParams={'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, auto_outgroup=AUTO_OUTGROUP, bootstrap_reps=bootstrap_reps, minClustSize=MinClustSize, convert_to_ultrametric=ULTRAMETRIC, bpflank=bpflank, combine_and_do_small_clusters=SMALLS)
#
#
if USEMCL:
	ltr_divergence(I=MCL_I, clustering_method='MCL', WickerParams=None)
	phylo(removegeneconv=False, BOOTSTRAP=True, I=MCL_I, align='cluster', removehomologouspair=RMHOMOFLANK, part='entire', clustering_method='MCL', WickerParams=None, auto_outgroup=AUTO_OUTGROUP,  bootstrap_reps=bootstrap_reps, minClustSize=MinClustSize, convert_to_ultrametric=ULTRAMETRIC, bpflank=bpflank, combine_and_do_small_clusters=SMALLS)
#
print('Fin!')
sys.exit()
div2Rplots(I=MCL_I)
#
#print('Done with estimate_divergence()')
#sys.exit()
#
#ts_tv() # Calculate transition/transversion for elements
#
#phylo(removegeneconv=False, BOOTSTRAP=True, I=MCL_I, align='clusters', removehomologouspair=True, part='entire') # Run FastTree for selections of elements/clusters
#
#xgmml() # for Cytoscape import
#
#hive()
#
#jbrowse() # track for jBrowse
#
#Rdensity() # ggplot2 density plot or histogram
#
#printPhylo() # ETE tree
#
print('made it')

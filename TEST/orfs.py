#!/usr/bin/env python3

import sys
import subprocess
from Bio import SeqIO

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



def bestORFs(fasta, strandsDct, outdir):
	'''
	Finds all ORFs in fasta using EMBOSS getorf and writes the best ones to a GFF3 and protein FASTA.
	The best ORFs are the set of non-overlapping ORFs containing the longest ORF out of sets on that
	strand,
	given by strandsDct, or, if strand is unknown, the set containing

	The coordinates output by getorf are 1-based
	'''
	outseq = '{0}/{1}.orfs'.format(outdir, fasta.split('/')[-1])
	#getorf_call = [ executables['getorf'], '-sequence', fasta, '-outseq', outseq ]
	#makecall(getorf_call)
	subprocess.call('/home/derstudent/software/EMBOSS-6.6.0/bin/getorf -sequence {0} -outseq {1}'.format(fasta, outseq), shell=True)
	orfs = list(SeqIO.parse(outseq, 'fasta'))
	print(orfs[0].description.split()[1][1:])
	print(orfs[-1].description.split()[3][:-1])
	print(orfs[-1].id)
	print()
	print(orfs[0])
	# Store needed information
	orfs_ordered_length_dct = {}
	orfs_ordered_coord_dct = {}
	for orf in orfs:
		element = '_'.join(orf.id.split('_')[:-1])
		orfnum = orf.id.split('_')[-1]
		desc = orf.description
		if 'REVERSE' in desc:
			strand = '-'
			end = int(desc.split()[1][1:])
			start = int(desc.split()[3][:-1])
		else:
			strand = '+'
			start = int(desc.split()[1][1:])
			end = int(desc.split()[3][:-1])

		if element in orfs_ordered_length_dct
			orfs_ordered_length_dct[element].append((desc, start-end+1)) # name, length
		else:
			orfs_ordered_length_dct[element] = (desc, start-end+1) # name, length
		
		if element in orfs_ordered_coord_dct:

	# For each element, find ORFs for its strand or both strands if strand = ?

	# Order ORFs names in list from element/strand longest to shortest.

	# Order another ORF name list with coordinate position of starts, smallest to largest

	# Make a dict with the sequences to use after selecting which orfs to keep

	i = 0
	while len(orfs_ordered_length) > i: # orfs_ordered_length is a list that is modified. i gets incremented
		el = orfs_ordered_length[i]
		j = orfs_ordered_coord.index(el)
		k = j+1
		m = j-1
		if k < len(orfs_ordered_coord)-1:
			while Overlaps( orfs_ordered_coord[j], orfs_ordered_coord[k] ): # write this function
				orfs_ordered_coord.pop(orfs_ordered_coord.index(orfs_ordered_coord[k]))
				orfs_ordered_length.pop(orfs_ordered_length.index(orfs_ordered_length[k]))
				k += 1
		if m < len(orfs_ordered_coord)+1:
			while Overlaps( orfs_ordered_coord[j], orfs_ordered_coord[k] ): # write this function
				orfs_ordered_coord.pop(orfs_ordered_coord.index(orfs_ordered_coord[k]))
				orfs_ordered_length.pop(orfs_ordered_length.index(orfs_ordered_length[k]))
				m -= 1

		i += 1
	
	print(orfs_ordered_length)

	# ORFs in dict with el:coords





args = sys.argv
outdir = '.'
strandsDct = {}
fasta = 'LTRharvest_LTR_retrotransposons.fasta'
bestORFs(fasta, strandsDct, outdir)



#>LTR_retrotransposon1_1 [19 - 63] 
#QVVMMKQFPGLLLQW
#>LTR_retrotransposon1_2 [35 - 73] 
#SSFLACCCSGNRW
#>LTR_retrotransposon1_3 [107 - 145] 
#FRRDRRGRRDRTS
#>LTR_retrotransposon1_4 [154 - 195] 
#FRSQKNQRWNYWLE
#>LTR_retrotransposon1_5 [149 - 202] 
#TSFALRRTRGGTTGWNDG
#>LTR_retrotransposon1_6 [215 - 271] 
#FREYRSGGCCRRNDRGRNE
#>LTR_retrotransposon1_7 [331 - 390] 
#LEKKCLRKKLLRRKSQKFHC
#>LTR_retrotransposon1_8 [3 - 395] 
#VVVLVTGGNDEAVSWLVVAVVIGGNDEDEYVDFGDSEEIEEEEETELAELVSLSEEPEVE
#LLVGMTDDLNDSENIEVEVVVEEMTGEEMNEEEMTGEEVTREEVTEEEMTREEVSEEEIT
#EEEITEVPLLG
#>LTR_retrotransposon1_9 [365 - 448] 
#GGNHRSSIARIEWGNNWMLRRDRAGYLF
#>LTR_retrotransposon78_576 [360 - 253] (REVERSE SENSE) 
#IYKSRFLTKPGPRRRNRGRGWVIITRILNIYRAHIE
#>LTR_retrotransposon78_577 [254 - 213] (REVERSE SENSE) 
#NKFYNIINKYIIRI
#>LTR_retrotransposon78_578 [187 - 155] (REVERSE SENSE) 
#SPRPGGDPEGV
#>LTR_retrotransposon78_579 [228 - 142] (REVERSE SENSE) 
#IYYKNIIYYIIKYKVPAPAGTPKEYKRYN
#>LTR_retrotransposon78_580 [194 - 129] (REVERSE SENSE) 
#NIKSPPRRGPRRSINDIINYII
#>LTR_retrotransposon78_581 [130 - 92] (REVERSE SENSE) 
#YKYKLKIIINLIK
#>LTR_retrotransposon78_582 [85 - 56] (REVERSE SENSE) 
#MINKKISGSQ
#>LTR_retrotransposon78_583 [78 - 43] (REVERSE SENSE) 
#TRRYPGPNNNYY
#>LTR_retrotransposon78_584 [49 - 11] (REVERSE SENSE) 
#LLLKIIIGTPTIE
#>LTR_retrotransposon78_585 [119 - 3] (REVERSE SENSE) 
#IKNNNKFNKIINDKQEDIRVPIIIIIENNNWDPHNRIEN
#

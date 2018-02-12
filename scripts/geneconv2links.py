#!/usr/bin/env python3

import sys

args = sys.argv

def geneconv2circoslinks(geneconvfile, ltrharvestgff):
	'''
	Converts GI tract pairs from geneconvClusters() output and writes a links file for Circos
	The GFF3 is needed to get the scaffold name.
	seqlengths needs to be a dictionary with the lengths of the sequences whose names correspond
	to the sequence names in the gff for the features with gene conversion tracts.
	'''
	global paths

	#lengths = { l[0]:int(l[1]) for  l in open(seqlengths).read().strip().split('\n')}

	seqs = {}
	with open(ltrharvestgff, 'r') as inFl:
		for line in inFl:
			if line.startswith('#'):
				continue
			if '\tLTR_retrotransposon\t' in line:
				gffLine = GFF3_line(line)
				element = gffLine.attributes['ID']
				scaf = gffLine.seqid
				seqs[element] = scaf
	
	with open(geneconvfile, 'r') as inFl:
		for line in inFl:
			if line.startswith('GI'):
				rec = line.strip().split('\t')
				el1, el2 = [ 'LTR_retrotransposon{0}'.format(e[1:]) for e in rec[1].split(';') ]
				el1start  = int(rec[7])
				el1end  = int(rec[8])
				el2start  = int(rec[10])
				el2end  = int(rec[11])
				el1seq = seqs[el1]
				el2seq = seqs[el2]



# Frag	Sequence	Sim	BC KA		Aligned Offsets			In Seq1			In Seq2		Num	Num	Total	Mismatch	
#   type	Names	Pvalue	Pvalue	Begin	End	Length	Begin	End	Length	Begin	End	Length	Poly	Dif	Diffs	Penalty	
GI	S554;S877	0.0016	> 1.0	be	en	le	be	en	len	be	en	le	pol	nd	ndd	pen	Other	0
GI	S1356;S999	0.0016	> 1.0	43	62	20	31	50	20	41	60	20	20	0	10	None	Other	0
GI	S600;S860	0.0018	> 1.0	45	60	16	18	33	16	43	58	16	16	0	12	None	Other	0
GI	S1106;S999	0.0018	> 1.0	41	62	22	39	60	22	39	60	22	22	0	9	None	Other	0
GI	S1380;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
GI	S1433;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
GI	S1515;S877	0.0018	> 1.0	41	62	22	39	60	22	26	47	22	22	0	9	None	Other	0
GI	S521;S999	0.0018	> 1.0	41	62	22	27	48	22	39	60	22	22	0	9	None	Other	0

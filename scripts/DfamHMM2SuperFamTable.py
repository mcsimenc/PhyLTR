#!/usr/bin/env python3

import sys

args = sys.argv

if '-h' in args:
	print('''
		usage:
			DfamHMM2SuperFamTable.py < <Dfam.hmm>

		description:
			Takes Dfam.hmm on stdin and returns a 2-column list
			on stdout with (1) element name and (2) superfamily

		''')
	sys.exit()


name = None
superfamily = None
for line in sys.stdin:
	if line.startswith('NAME'):
		name = line.strip().split(' ')[-1]
	elif line.startswith('CT'):
		if 'Superfamily' in line:
			superfamily = line.strip().split(';')[-2]
			if name == None or superfamily == None:
				sys.exit('Problem: a record does not have superfamily or name entry.\nNAME found: {0},Superfamily found:{1}'.format(name,superfamily))
			else:
				print(name, superfamily, sep='\t')

# Example Dfam.hmm record portion:

#NAME  ACCORD2_I
#ACC   DF0001530.0
#DESC  ACCORD2_I is an internal portion of ACCORD2 endogenous retrovirus.
#LENG  7212
#MAXL  7906
#ALPH  DNA
#RF    yes
#MM    no
#CONS  yes
#CS    no
#MAP   yes
#DATE  Mon Aug 17 21:19:27 2015
#NSEQ  237
#EFFN  2.038221
#CKSUM 2330758380
#GA    19.31;
#TC    19.31;
#NC    19.31;
#TH    TaxId:7227; TaxName:Drosophila melanogaster; GA:1.47; TC:19.31; NC:1.40; fdr:0.002;
#BM    hmmbuild --eentexp --ere 0.62 --maxinsertlen 10 --hand  --dna --fragthresh 1HMM.ann SEED.ann
#SM    nhmmer -Z 143 --dfamtblout dm6-full_hits  -E 100 --noali HMM.ann dfamseq
#CT    Type; Retrotransposon;
#CT    Class; LTR;
#CT    Superfamily; Gypsy;
